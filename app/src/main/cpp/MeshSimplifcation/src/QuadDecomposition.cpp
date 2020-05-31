/*
** Copyright (c) 2019 Qualcomm Technologies, Inc.
** All Rights Reserved.
** Confidential and Proprietary - Qualcomm Technologies, Inc.
**
** Code from https://github.com/rawrunprotected/rasterizer/blob/master/SoftwareRasterizer/QuadDecomposition.cpp under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication.
*/

#include "QuadDecomposition.h"
#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>
#include <random>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "MathUtil.h"

using namespace common;

typedef int Vertex;
typedef std::vector<Vertex> Path;

namespace
{
    struct PairHash
    {
    public:
        template <typename T, typename U>
        std::size_t operator()(const std::pair<T, U> &x) const
        {
            auto hashT = std::hash<T>{}(x.first);
            auto hashU = std::hash<U>{}(x.second);
            return hashT ^ (hashU + 0x9e3779b9 + (hashT << 6) + (hashT >> 2));
        }
    };

    bool canMergeTrianglesToQuad(__m128 v0, __m128 v1, __m128 v2, __m128 v3)
    {
        // Maximum distance of vertices from original plane in world space units
        float maximumDepthError = 0.5f;

        __m128 n0 = normalize(normal(v0, v1, v2));
        __m128 planeDistA = _mm_andnot_ps(_mm_set1_ps(-0.0f), _mm_dp_ps(n0, _mm_sub_ps(v1, v3), 0x7F));
        if (_mm_comigt_ss(planeDistA, _mm_set1_ps(maximumDepthError)))
        {
            return false;
        }

        __m128 n2 = normalize(normal(v2, v3, v0));
        __m128 planeDistB = _mm_andnot_ps(_mm_set1_ps(-0.0f), _mm_dp_ps(n2, _mm_sub_ps(v1, v3), 0x7F));
        if (_mm_comigt_ss(planeDistB, _mm_set1_ps(maximumDepthError)))
        {
            return false;
        }
        return true;
    }
} // namespace

namespace common
{
    struct MeshVertex;
    struct MeshFaceEdge;
    struct MeshQuad;
    struct MeshFace;

    //this edge is owned by vertex
    struct MeshFaceEdge
    {
        int vertIdx; //owned by vertIdx in the face
        MeshFace * face = nullptr;//the face it belong

        MeshVertex *start;
        MeshVertex *end;
        MeshVertex *apex;

        MeshFaceEdge*next = nullptr; //a linked list to store other Edge owned by the same vertex
    };

    struct MeshVertex
    {
        MeshVertex() {}
        void init(int idx)
        {
            vIdx = idx;
            edgeList = nullptr;
        }
        int vIdx;

        MeshFaceEdge * edgeList= nullptr; // edge linked list with start point this vertex
        void addVertexEdge(MeshFaceEdge * e)
        {
            if (edgeList == nullptr) edgeList = e;
            else
            {
                e->next = edgeList;
                edgeList = e;
            }
        }
    };

    struct Node
    {
        Node()
        {
            m_clearToken = 0;
        }

        int m_depth;
        Vertex m_parent;
        Vertex m_blossom;

        int m_clearToken;

        Vertex bridgeV;
        Vertex bridgeW;

        void reset() {
            m_depth = 0;
            m_parent = m_blossom = 0;
            bridgeV = bridgeW = 0;
            m_clearToken = 0;
        }
    };

    struct MeshFace
    {
        MeshFace() {}
        void init(int id, MeshVertex* v0, MeshVertex* v1, MeshVertex* v2)
        {
            faceId = id;

            Vertices[0] = v0;
            Vertices[1] = v1;
            Vertices[2] = v2;
            for (int i = 0; i < 3; i++)
            {
                auto e = &Edges[i];
                e->face = this;
                e->start = Vertices[i];
                e->vertIdx = e->start->vIdx;
                e->end = Vertices[(i+1)%3];
                e->apex = Vertices[(i + 2) % 3];
                e->next = nullptr;
                Vertices[i]->addVertexEdge(e);
            }
            pairFace = -1;
            neighbors.clear();
            GraphNode.reset();
        }
        int faceId = 0;
        MeshFaceEdge Edges[3];

        MeshVertex* Vertices[3] = { nullptr };

        __m128 unitNormal;// 0 = normalize(normal(v0, v1, v2));

        Node GraphNode;
        std::vector<Vertex> neighbors;
        int pairFace = -1;
    };

    struct MeshQuad
    {
        MeshFaceEdge * faceEdge1;
        MeshFaceEdge * faceEdge2;

        MeshQuad() {}
        void init(MeshFaceEdge*e1, MeshFaceEdge*e2)
        {
            this->faceEdge1 = e1;
            this->faceEdge1 = e2;
        }

        __m128 faceNormalSum;
        //unsigned int bestCentroid;
        float* pNormal;

        static inline bool compareX(MeshQuad* const &rhs, MeshQuad* const &lhs)
        {
            return (rhs->pNormal[0] < lhs->pNormal[0]);
        }
        static inline bool compareY(MeshQuad* const &rhs, MeshQuad* const &lhs)
        {
            return (rhs->pNormal[1] < lhs->pNormal[1]);
        }
        static inline bool compareZ(MeshQuad* const &rhs, MeshQuad* const &lhs)
        {
            return (rhs->pNormal[2] < lhs->pNormal[2]);
        }

        int getFirstVertex()
        {
            if (faceEdge1 != faceEdge2) return faceEdge1->start->vIdx;
            return faceEdge1->face->Vertices[0]->vIdx;
        }
    };

    static	std::vector<MeshFace*> Faces;
    static	std::vector<MeshVertex*> Vertices;
    static	std::vector<MeshQuad*> Quads;
    static	std::vector<MeshQuad*> QuadsKeys;

    class Matching
    {
    public:
        Matching(int faceNum) : m_clearToken(0)
        {
            // Start with a greedy maximal matching
            for (Vertex v = 0; v < faceNum; ++v)
            {
                if (Faces[v]->pairFace == -1)
                {
                    for (auto w : Faces[v]->neighbors)
                    {
                        if (Faces[w]->pairFace == -1)
                        {
                            match(v, w);
                            break;
                        }
                    }
                }
            }
            std::vector<Vertex> path;
            for (Vertex v = 0; v < faceNum; ++v)
            {
                if (Faces[v]->pairFace == -1)
                {
                    if (findAugmentingPath(v, path))
                    {
                        augment(path);
                        path.clear();
                    }
                }
            }
        }

        Vertex getMatchedVertex(Vertex v)
        {
            return Faces[v]->pairFace;
        }

    private:
        void match(Vertex v, Vertex w)
        {
            Faces[v]->pairFace = w;
            Faces[w]->pairFace = v;
        }

        void augment(std::vector<Vertex> &path)
        {
            for (int i = 0; i < path.size(); i += 2)
            {
                match(path[i], path[i + 1]);
            }
        }

        bool findAugmentingPath(Vertex root, std::vector<Vertex> &path)
        {
            // Clear out the forest
            m_clearToken++;

            // Start our tree root
            Faces[root]->GraphNode.m_depth = 0;
            Faces[root]->GraphNode.m_parent = -1;
            Faces[root]->GraphNode.m_clearToken = m_clearToken;
            Faces[root]->GraphNode.m_blossom = root;

            m_queue.push(root);

            while (!m_queue.empty())
            {
                Vertex v = m_queue.front();
                m_queue.pop();

                for (auto w : Faces[v]->neighbors)
                {
                    if (examineEdge(root, v, w, path))
                    {
                        while (!m_queue.empty())
                        {
                            //std::cout << "clearing " << m_queue.size() << std::endl;
                            std::queue<Vertex> empty;
                            std::swap(m_queue, empty);
                        }

                        return true;
                    }
                }
            }

            return false;
        }

        bool examineEdge(Vertex root, Vertex v, Vertex w, std::vector<Vertex> &path)
        {
            Vertex vBar = find(v);
            Vertex wBar = find(w);

            if (vBar != wBar)
            {
                if (Faces[wBar]->GraphNode.m_clearToken != m_clearToken)
                {
                    if (Faces[w]->pairFace == -1)
                    {
                        buildAugmentingPath(root, v, w, path);
                        return true;
                    }
                    else
                    {
                        extendTree(v, w);
                    }
                }
                else if (Faces[wBar]->GraphNode.m_depth % 2 == 0)
                {
                    shrinkBlossom(v, w);
                }
            }

            return false;
        }

        void buildAugmentingPath(Vertex root, Vertex v, Vertex w, std::vector<Vertex> &path)
        {
            path.push_back(w);
            findPath(v, root, path);
        }

        void extendTree(Vertex v, Vertex w)
        {
            Vertex u = Faces[w]->pairFace;

            Node &nodeV = Faces[v]->GraphNode;
            Node &nodeW = Faces[w]->GraphNode;
            Node &nodeU = Faces[u]->GraphNode;

            nodeW.m_depth = nodeV.m_depth + 1 + (nodeV.m_depth & 1); // Must be odd, so we add either 1 or 2
            nodeW.m_parent = v;
            nodeW.m_clearToken = m_clearToken;
            nodeW.m_blossom = w;

            nodeU.m_depth = nodeW.m_depth + 1;
            nodeU.m_parent = w;
            nodeU.m_clearToken = m_clearToken;
            nodeU.m_blossom = u;

            m_queue.push(u);
        }

        void shrinkBlossom(Vertex v, Vertex w)
        {
            Vertex b = findCommonAncestor(v, w);

            shrinkPath(b, v, w);
            shrinkPath(b, w, v);
        }

        void shrinkPath(Vertex b, Vertex v, Vertex w)
        {
            Vertex u = find(v);

            while (u != b)
            {
                makeUnion(b, u);
                assert(u != -1);
                assert(Faces[u]->pairFace != -1);
                u = Faces[u]->pairFace;
                makeUnion(b, u);
                makeRepresentative(b);
                m_queue.push(u);
                Faces[u]->GraphNode.bridgeV = v;
                Faces[u]->GraphNode.bridgeW = w;

                //m_bridges[u] = std::make_pair(v, w);
                u = find(Faces[u]->GraphNode.m_parent);
            }
        }

        Vertex findCommonAncestor(Vertex v, Vertex w)
        {
            while (w != v)
            {
                if (Faces[v]->GraphNode.m_depth > Faces[w]->GraphNode.m_depth)
                {
                    v = Faces[v]->GraphNode.m_parent;
                }
                else
                {
                    w = Faces[w]->GraphNode.m_parent;
                }
            }

            return find(v);
        }

        void findPath(Vertex s, Vertex t, Path &path)
        {
            if (s == t)
            {
                path.push_back(s);
            }
            else if (Faces[s]->GraphNode.m_depth % 2 == 0)
            {
                path.push_back(s);
                path.push_back(Faces[s]->pairFace);
                findPath(Faces[Faces[s]->pairFace]->GraphNode.m_parent, t, path);
            }
            else
            {
                Vertex v = 0, w = 0;
                v = Faces[s]->GraphNode.bridgeV;
                w = Faces[s]->GraphNode.bridgeW;
                //            std::tie(v, w) = m_bridges[s];

                path.push_back(s);

                size_t offset = path.size();
                findPath(v, Faces[s]->pairFace, path);
                std::reverse(path.begin() + offset, path.end());

                findPath(w, t, path);
            }
        }

        void makeUnion(int x, int y)
        {
            int xRoot = find(x);
            Faces[xRoot]->GraphNode.m_blossom = find(y);
        }

        void makeRepresentative(int x)
        {
            int xRoot = find(x);
            Faces[xRoot]->GraphNode.m_blossom = x;
            Faces[x]->GraphNode.m_blossom = x;
        }

        int find(int x)
        {
            if (Faces[x]->GraphNode.m_clearToken != m_clearToken)
            {
                return x;
            }

            if (x != Faces[x]->GraphNode.m_blossom)
            {
                // Path compression
                Faces[x]->GraphNode.m_blossom = find(Faces[x]->GraphNode.m_blossom);
            }

            return Faces[x]->GraphNode.m_blossom;
        }

    private:
        int m_clearToken;


        std::queue<Vertex> m_queue;



        //std::vector<std::pair<Vertex, Vertex>> m_bridges;
    };

    static const unsigned int PoolUnit = 10;


    static void requestSpawn(int vertCount, int triangleCount)
    {
        while(Faces.size() < triangleCount)
        {
            MeshFace* f = new MeshFace[1 << PoolUnit];
            for (int idx = 0; idx < 1 << PoolUnit; idx++)
                Faces.push_back(&f[idx]);

        }

        while (Vertices.size() < vertCount)
        {
            MeshVertex* v = new MeshVertex[1 << PoolUnit];
            for (int idx = 0; idx < 1 << PoolUnit; idx++)
                Vertices.push_back(&v[idx]);

        }
    }
    static MeshQuad* spawnMeshQuad(int idx)
    {
        while (Quads.size() <= idx)
        {
            MeshQuad* q = new MeshQuad[1<< PoolUnit];
            QuadsKeys.push_back(&q[0]);

            for (int idx = 0; idx < 1 << PoolUnit; idx++)
                Quads.push_back(&q[idx]);
        }
        return Quads[idx];
    }

    void QuadDecomposition::release()
    {
        for (int idx = 0; idx < QuadsKeys.size(); idx++)
            delete[] QuadsKeys[idx];
        QuadsKeys.clear();
        Quads.clear();

        for (int idx = 0; idx < Vertices.size()>> PoolUnit; idx++)
            delete[] Vertices[idx<< PoolUnit];
        Vertices.clear();

        for (int idx = 0; idx < Faces.size() >> PoolUnit; idx++)
            delete[] Faces[idx << PoolUnit];
        Faces.clear();
    }



    static bool canMergeTrianglesToQuadV2(MeshFace* f1, MeshFace*f2,  __m128 v1, __m128 v3)
    {
        // Maximum distance of vertices from original plane in world space units
        float maximumDepthError = 0.5f;

        auto v1v3 = _mm_sub_ps(v1, v3);
        __m128 planeDistA = _mm_andnot_ps(_mm_set1_ps(-0.0f), _mm_dp_ps(f1->unitNormal, v1v3, 0x7F));
        if (_mm_comigt_ss(planeDistA, _mm_set1_ps(maximumDepthError)))
        {
            return false;
        }

        __m128 planeDistB = _mm_andnot_ps(_mm_set1_ps(-0.0f), _mm_dp_ps(f2->unitNormal, v1v3, 0x7F));
        if (_mm_comigt_ss(planeDistB, _mm_set1_ps(maximumDepthError)))
        {
            return false;
        }

        return true;
    }

    static void matchFaceEdge(MeshFace * face1, MeshFace * face2, MeshFaceEdge *& faceEdge1, MeshFaceEdge *& faceEdge2)
    {
        for (int idx = 0; idx < 3; idx++)
        {
            faceEdge1 = &face1->Edges[idx];
            for (int j = 0; j < 3; j++)
            {
                faceEdge2 = &face2->Edges[j];
                if (faceEdge1->start == faceEdge2->end &&
                    faceEdge1->end == faceEdge2->start)
                {
                    return;
                }
            }
        }
        faceEdge1 = nullptr;
        faceEdge2 = nullptr;
    }

    static void fillInVertices(__m128 * orderedVertices, MeshQuad * q, const std::vector<__m128> &vertices)
    {
        auto faceEdge1 = q->faceEdge1;
        auto faceEdge2 = q->faceEdge2;
        if (faceEdge1 != faceEdge2)
        {
            orderedVertices[0] = vertices[faceEdge1->start->vIdx];
            orderedVertices[1] = vertices[faceEdge1->apex->vIdx];
            orderedVertices[2] = vertices[faceEdge1->end->vIdx];
            orderedVertices[3] = vertices[faceEdge2->apex->vIdx];
        }
        else
        {
            auto& vs = faceEdge1->face->Vertices;
            orderedVertices[0] = vertices[vs[0]->vIdx];
            orderedVertices[1] = vertices[vs[2]->vIdx];
            orderedVertices[2] = vertices[vs[1]->vIdx];
            orderedVertices[3] = vertices[vs[0]->vIdx];
        }
    }

    std::vector<uint32_t> QuadDecomposition::decompose(const std::vector<uint32_t> &indices, const std::vector<__m128> &vertices)
    {
        std::vector<uint32_t> result;
        size_t triangleCount = indices.size() / 3;
        size_t vertCount = vertices.size();

        requestSpawn(vertCount, triangleCount);

        for (int idx = 0; idx < vertCount; idx++)
        {
            Vertices[idx]->init(idx);
        }
        int vIdx = 0;
        for (int idx = 0; idx < triangleCount; idx++)
        {
            auto v0 = Vertices[indices[vIdx++]];
            auto v1 = Vertices[indices[vIdx++]];
            auto v2 = Vertices[indices[vIdx++]];
            Faces[idx]->init(idx, v0, v1, v2);
        }

        for (uint32_t triangleIdx = 0; triangleIdx < triangleCount; ++triangleIdx)
        {
            MeshFace * f = Faces[triangleIdx];
            f->unitNormal =  normalize(normal(vertices[f->Vertices[0]->vIdx],
                                              vertices[f->Vertices[1]->vIdx],
                                              vertices[f->Vertices[2]->vIdx]));
        }
        int facePairIdx = 0;
        for (uint32_t triangleIdx = 0; triangleIdx < triangleCount; ++triangleIdx)
        {
            MeshFace * f = Faces[triangleIdx];

            for (int edgeIdx = 0; edgeIdx < 3; ++edgeIdx)
            {
                auto fv = f->Vertices[edgeIdx];
                MeshFaceEdge* fe = &f->Edges[edgeIdx];

                auto neighbors = fe->end->edgeList;
                for (auto e = neighbors; e != nullptr; e = e->next)
                {
                    if (e->end == fv && e->face->faceId <= f->faceId)
                    {
                        if (canMergeTrianglesToQuadV2(f, e->face,
                                                      vertices[fe->apex->vIdx],
                                                      vertices[e->apex->vIdx]))
                        {
                            f->neighbors.push_back(e->face->faceId);
                            e->face->neighbors.push_back(triangleIdx);
                        }
                    }
                }
            }
        }

        Matching matching(triangleCount);
        int quadIdx = 0;
        for (uint32_t triangleIdx = 0; triangleIdx < triangleCount; ++triangleIdx)
        {
            int neighbor = matching.getMatchedVertex(triangleIdx);

            // No quad found
            if (neighbor == -1)
            {
                auto q = spawnMeshQuad(quadIdx++); //0  2 1 0
                q->faceEdge1 = &Faces[triangleIdx]->Edges[0];
                q->faceEdge2 = q->faceEdge1;
                result.push_back(q->faceEdge1->start->vIdx);
                result.push_back(q->faceEdge1->apex->vIdx);
                result.push_back(q->faceEdge1->end->vIdx);
                result.push_back(q->faceEdge1->start->vIdx);
            }
            else if (triangleIdx < uint32_t(neighbor))
            {
                MeshFace * face1 = Faces[triangleIdx];
                MeshFace * face2 = Faces[neighbor];
                MeshFaceEdge * faceEdge1;
                MeshFaceEdge * faceEdge2;
                matchFaceEdge(face1, face2, faceEdge1, faceEdge2);

                if (faceEdge1 != nullptr && faceEdge2 != nullptr)
                {
                    auto q = spawnMeshQuad(quadIdx++);
                    q->faceEdge1 = faceEdge1;
                    q->faceEdge2 = faceEdge2;
                    result.push_back(q->faceEdge1->start->vIdx);
                    result.push_back(q->faceEdge1->apex->vIdx);
                    result.push_back(q->faceEdge1->end->vIdx);
                    result.push_back(q->faceEdge2->apex->vIdx);
                }
            }
        }

        return result;
    }



} // namespace common
