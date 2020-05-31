/*========================================================================*/
/**
\file MeshReducer.cpp

Copyright (c) 2019 Qualcomm Technologies, Inc.
All rights reserved.
Confidential and Proprietary - Qualcomm Technologies, Inc.
*/
/*========================================================================*/

#include "MeshReducer.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_map>
#include "MathUtil.h"

namespace
{
    using namespace common;
    static const unsigned int POOL_SIZE = 1024;
    static const unsigned int TIMEOUT = 5;
    static const unsigned int MAX_ITERATION = 100;
    static const double QUADRIC_EPSILON = 1e-15;
    static const double VALID_THRESHOLD = 0.9604;
    static const double AGGRESSIVE = 8.0;
    static const double AREA_TOLERANCE = 0.999;
    static const double NORMAL_TOLERANCE = 0.2;

    struct Edge;
    struct Face;
    struct FaceNode;

    struct Vertex
    {
        Vertex() : CurrentEdge(ToUpdateEdge) {}

        Vec3 Pos;
        uint16_t ID = 0;

        // attribute
        bool Removed = false;

        // Need to update?
        Edge *ToUpdateEdge = nullptr;

        // Aliasing name
        Edge *&CurrentEdge;

        // use LOCAL_MARK to update
        unsigned int LOCAL_MARK = 0;

        // Quadric Matrix
        SymetricMatrix Q;

        // adjacent face list
        FaceNode *Neighbors = nullptr;

        // adjacent edges
        std::list<Edge *> Edges;

        void reset()
        {
            Edges.clear();
            LOCAL_MARK = 0;
            ID = 0;
            ToUpdateEdge = nullptr;
            Neighbors = nullptr;
            Removed = false;
            Q.reset();
        }

        static inline bool comparefunc(Vertex* const &lhs, Vertex* const &rhs)
        {
            return (lhs->Pos < rhs->Pos);
        }
    };

    struct Edge
    {
        Edge() : NextEdge(ToBeReplaced) {}
        Edge(Vertex *v0, Vertex *v1) : Start(v0), End(v1), NextEdge(ToBeReplaced) {}

        // start & end vertex
        Vertex *Start = nullptr;
        Vertex *End = nullptr;

        // collapse property for edge-based simplification
        unsigned int LOCAL_MARK = 0;

        // Aliasing name
        unsigned int AdjFaces = 0;

        double Priority = 0.0;
        Vec3 OptPos;

        // heap related
        int HeapIndex = -1;
        Edge *ToBeReplaced = nullptr;

        // Aliasing name
        Edge *&NextEdge;

        void replaceVertex(Vertex *dst, Vertex *src)
        {
            if (Start == dst)
            {
                Start = src;
            }
            else if (End == dst)
            {
                End = src;
            }
        }

        bool containVertex(Vertex *v)
        {
            return (Start == v) || (End == v);
        }

        void init(Vertex *start, Vertex *end)
        {
            // update vertices
            Start = start;
            End = end;

            // adjacent faces
            AdjFaces = 0;

            // attributes
            LOCAL_MARK = 0;
            Priority = 0.0;
            HeapIndex = -1;
            ToBeReplaced = nullptr;
        }
    };

    struct FaceNode
    {
        Face *Instance = nullptr;
        FaceNode *Next = nullptr;
    };

    struct Face
    {
        Face() {}

        void reset()
        {
            for (uint32_t i = 0; i < 3; ++i)
            {
                Vertices[i] = nullptr;
                Edges[i] = nullptr;
            }

            // attributes
            Valid = true;
            LOCAL_MARK = 0;
            OptEdge = 0;
        }

        // adjacent vertices
        Vertex *Vertices[3] = {nullptr};

        // edges
        Edge *Edges[3] = {nullptr};

        // valid & dirty(need to update)
        bool Valid = true;

        // related face node
        FaceNode FaceNodes[3];

        // collapse property for fast simplification
        uint16_t OptEdge = 0;
        unsigned int LOCAL_MARK = 0; //to check whether face are updated for current iteration

        // Priority
        double priority = 0;

        // get connected edges
        void getEdges(Vertex *v, Edge *&e1, Edge *&e2)
        {
            for (uint32_t i = 0; i < 3; ++i)
            {
                if (v == Vertices[i])
                {
                    e1 = Edges[(i + 2) % 3];
                    e2 = Edges[i];
                    break;
                }
            }
        }

        void setCentral(Vertex *v0, Vertex *&v1, Vertex *&v2) const
        {
            for (uint32_t i = 0; i < 3; ++i)
            {
                if (v0 == Vertices[i])
                {
                    v1 = Vertices[(i + 1) % 3];
                    v2 = Vertices[(i + 2) % 3];
                    break;
                }
            }
        }

        void replace(Vertex *dst, Vertex *src)
        {
            for (uint32_t i = 0; i < 3; ++i)
            {
                if (dst == Vertices[i])
                {
                    Vertices[i] = src;
                    Edges[i]->replaceVertex(dst, src);
                    Edges[(i + 2) % 3]->replaceVertex(dst, src);
                    break;
                }
            }
        }

        bool containVertex(Vertex *v) const
        {
            return Vertices[0] == v || Vertices[1] == v || Vertices[2] == v;
        }

        Vec3 normal() const
        {
            const Vec3 &v0 = Vertices[0]->Pos;
            const Vec3 &v1 = Vertices[1]->Pos;
            const Vec3 &v2 = Vertices[2]->Pos;

            const Vec3 &e0 = v1 - v0;
            const Vec3 &e1 = v2 - v0;

            return e0.cross(e1);
        }

        double computeQuality(Vertex *start, const Vec3 &optPos)
        {
            Vertex *v1 = nullptr;
            Vertex *v2 = nullptr;
            this->setCentral(start, v1, v2);

            const Vec3 &e0 = v2->Pos - v1->Pos;
            const Vec3 &e1 = v1->Pos - optPos;

            const Vec3 &normal = e0.cross(e1);

            double len_norm = normal.length();
            if (len_norm == 0.0)
            {
                return 0;
            }

            const Vec3 &e2 = v2->Pos - optPos;
            double len_e0 = e0.squaredLength();
            double len_e1 = e1.squaredLength();
            double len_e2 = e2.squaredLength();
            double max_edge = std::max(len_e0, std::max(len_e1, len_e2));

            return len_norm / max_edge;
        }

        // get opposite edges
        void getEdges(Vertex *v0, Vertex *v1, Edge *&v2v1, Edge *&v2v0)
        {
            if (Vertices[0] == v0)
            {
                if (Vertices[1] == v1)
                {
                    v2v1 = this->Edges[1];
                    v2v0 = this->Edges[2];
                }
                else
                {
                    v2v0 = this->Edges[0];
                    v2v1 = this->Edges[1];
                }
            }
            else if (Vertices[1] == v0)
            {
                if (Vertices[2] == v1)
                {
                    v2v1 = this->Edges[2];
                    v2v0 = this->Edges[0];
                }
                else
                {
                    v2v1 = this->Edges[2];
                    v2v0 = this->Edges[1];
                }
            }
            else
            {
                if (Vertices[0] == v1)
                {
                    v2v1 = this->Edges[0];
                    v2v0 = this->Edges[1];
                }
                else
                {
                    v2v1 = this->Edges[0];
                    v2v0 = this->Edges[2];
                }
            }
        }

        void replaceDuplicatedEdge()
        {
            for (uint32_t i = 0; i < 3; ++i)
            {
                if (Edges[i]->ToBeReplaced)
                {
                    Edges[i] = Edges[i]->ToBeReplaced;
                }
            }
        }

        Edge* replaceEdge(Vertex *start, Vertex *end)
        {
            Edge *sEdge = nullptr;
            Edge *eEdge = nullptr;
            getEdges(start, end, sEdge, eEdge);
            sEdge->ToBeReplaced = eEdge;
            return sEdge;
        }

        Vertex* getRestVertex(Vertex *v0, Vertex *v1)
        {
            if (Vertices[0] != v0 && Vertices[0] != v1)
            {
                return Vertices[0];
            }
            if (Vertices[1] != v0 && Vertices[1] != v1)
            {
                return Vertices[1];
            }
            return Vertices[2];
        }
    };

    struct EdgeHeap
    {
        EdgeHeap(std::vector<Edge*> &Edges) : Container(Edges), Length(Edges.size())
        {
            //how to get
            for (int idx = 1; idx < Length; ++idx)
            {
                initElement(idx);
            }
            for (int i = 0; i < Length; ++i)
            {
                Container[i]->HeapIndex = i;
            }
        }

        void initElement(int idx)
        {
            int parent = (idx - 1) >> 1;
            auto e = Container[idx];
            // heap up
            while (e->Priority < Container[parent]->Priority)
            {
                // swap and keep indexes updated
                Container[idx] = Container[parent];

                // top
                if (parent == 0)
                {
                    Container[0] = e;
                    return;
                }

                idx = parent;
                parent = (idx - 1) >> 1;
            }

            Container[idx] = e;
            return;
        }

        Edge *top()
        {
            Edge *item = Container[0];
            // mark as not inside heap anymore
            item->HeapIndex = -1;
            if (Length > 1)
            {
                // move the last one to top
                Container[0] = Container[Length - 1];
                Container[0]->HeapIndex = 0;
                // heap down
                heapifyDown(Container[0]);
            }

            // Remove the last element
            --Length;
            return item;
        }

        void update(Edge *e, double prev_priority)
        {
            if (e->Priority < prev_priority)
            {
                heapifyUp(e);
            }
            else if (e->Priority > prev_priority)
            {
                heapifyDown(e);
            }
        }

        void pop(Edge *edge)
        {
            // get current index
            int cur_idx = edge->HeapIndex;
            if (cur_idx < 0 || cur_idx >= Length || edge != Container[cur_idx])
            {
                return;
            }

            // mark as removed
            edge->HeapIndex = -2;
            if (cur_idx != (Length - 1))
            {
                auto last = Container[Length - 1];
                Container[cur_idx] = last;
                last->HeapIndex = cur_idx;
                if (last->Priority < edge->Priority)
                {
                    heapifyUp(last);
                }
                else if (last->Priority > edge->Priority)
                {
                    heapifyDown(last);
                }
            }

            // remove
            --Length;
        }

        void heapifyUp(Edge *e)
        {
            if (e->HeapIndex == 0)
            {
                return;
            }

            // get current index & parent
            int idx = e->HeapIndex;
            int parent = (idx - 1) >> 1;

            // heap up
            while (e->Priority < Container[parent]->Priority)
            {
                // swap
                std::swap(Container[idx], Container[parent]);
                Container[idx]->HeapIndex = idx;
                Container[parent]->HeapIndex = parent;

                // top
                if (parent == 0)
                {
                    return;
                }

                // update index
                idx = parent;
                parent = (idx - 1) >> 1;
            }
        }

        void heapifyDown(Edge *e)
        {
            int pos = e->HeapIndex;

            // get left & right children
            int left = (pos << 1) + 1;
            if (left >= Length)
            {
                return;
            }

            int right = left + 1;
            int smallest = pos;
            while (1)
            {
                if (Container[left]->Priority < e->Priority)
                {
                    if (right < Length && Container[right]->Priority < Container[left]->Priority)
                    {
                        smallest = right;
                    }
                    else
                    {
                        smallest = left;
                    }
                }
                else if (right < Length && Container[right]->Priority < e->Priority)
                {
                    smallest = right;
                }
                else
                {
                    Container[pos] = e;
                    e->HeapIndex = pos;
                    return;
                }

                // swap
                std::swap(Container[pos], Container[smallest]);
                Container[pos]->HeapIndex = pos;
                Container[smallest]->HeapIndex = smallest;

                left = (smallest << 1) + 1;
                if (left >= Length)
                {
                    return;
                }

                pos = smallest;
                right = left + 1;
            }
        }

        std::vector<Edge *> &Container;
        int Length = 0;
    };

    class CollapseHelper {
    public:
        static unsigned int GLOBAL_MARK;
        static double ScaleFactor;

        // Face pool to reduce memory reallocation
        static std::vector<Face *> FacePool;

        // Edge pool to reduce memory reallocation
        static std::vector<Edge *> EdgePool;
        static uint32_t EdgePoolIdx;

        // Edge pool to reduce memory reallocation
        static std::vector<Vertex *> VertexPool;

        static void reset()
        {
            GLOBAL_MARK = 0;
            ScaleFactor = 1.0;
            EdgePoolIdx = 0;
        }

        static Edge *spawnEdgeFromPool(Vertex *v0, Vertex *v1)
        {
            if (EdgePoolIdx >= EdgePool.size())
            {
                Edge* pool = new Edge[POOL_SIZE];
                for (int idx = 0; idx < POOL_SIZE; idx++)
                {
                    EdgePool.push_back(&pool[idx]);
                }
            }

            Edge* e = EdgePool[EdgePoolIdx++];
            e->init(v0, v1);
            return e;
        }

        static void reservePool(int nVert, int nFace)
        {
            VertexPool.reserve(nVert);
            FacePool.reserve(nFace);

            // allocate memory for vertex
            while (VertexPool.size() < nVert)
            {
                Vertex *v = new Vertex[POOL_SIZE];
                for (unsigned int i_vert = 0; i_vert < POOL_SIZE; ++i_vert)
                {
                    VertexPool.push_back(&v[i_vert]);
                }
            }

            // reset
            for (unsigned int i_vert = 0; i_vert < nVert; ++i_vert)
            {
                VertexPool[i_vert]->reset();
            }

            // allocate memory for face
            while (FacePool.size() < nFace)
            {
                Face *f = new Face[POOL_SIZE];
                for (unsigned int i_face = 0; i_face < POOL_SIZE; ++i_face)
                {
                    FacePool.push_back(&f[i_face]);
                }
            }

            // reset
            for (unsigned int i_face = 0; i_face < nFace; ++i_face)
            {
                FacePool[i_face]->reset();
            }
        }

        // cost = VT * Q * V
        static double getQuadricCost(const Vec3 &v, const SymetricMatrix &Q)
        {
            double x = v.X;
            double y = v.Y;
            double z = v.Z;
            // Q[0] * x * x + 2 * Q[1] * x * y + 2 * Q[2] * x * z + 2 * Q[3] * x + Q[4] * y * y + 2 * Q[5] * y * z + 2 * Q[6] * y + Q[7] * z * z + 2 * Q[8] * z + Q[9];
            return (Q[0] * x + 2 * Q[1] * y + 2 * Q[2] * z + 2 * Q[3]) * x +
                   (Q[4] * y + 2 * Q[5] * z + 2 * Q[6]) * y + (Q[7] * z + 2 * Q[8]) * z + Q[9];
        }

        // priority = cost / (normal * tri_quality)
        static double computePriority(Edge &e, const double &QuadricCost)
        {
            const Vec3 &optPos = e.OptPos;
            auto start = e.Start;
            auto end = e.End;

            // for each face related to start vertex
            double minQual = 1;
            for (auto node = start->Neighbors; node != nullptr; node = node->Next)
            {
                auto face = node->Instance;
                if (!face->Valid || face->containVertex(end))
                {
                    continue;
                }

                double quality = face->computeQuality(start, optPos);
                if (quality == 0)
                {
                    return std::numeric_limits<double>::infinity();
                }
                minQual = std::min(minQual, quality);
            }

            // for each face related to end vertex
            for (auto node = end->Neighbors; node != nullptr; node = node->Next)
            {
                auto face = node->Instance;
                if (!face->Valid || face->containVertex(start))
                {
                    continue;
                }

                double quality = face->computeQuality(end, optPos);
                if (quality == 0)
                {
                    return std::numeric_limits<double>::infinity();
                }
                minQual = std::min(minQual, quality);
            }

            // cost
            double cost = ScaleFactor * QuadricCost;
            if (cost <= QUADRIC_EPSILON)
            {
                cost = -1 / (start->Pos - end->Pos).length();
            }
            cost /= minQual;
            return cost;
        }

        static void calcOptimalPosition(Edge &e, double &cost)
        {
            static const double COST_THRESHOLD = 200.0 * QUADRIC_EPSILON;
            auto start = e.Start;
            auto end = e.End;

            const SymetricMatrix &Q = start->Q + end->Q;
            Vec3 &optPos = e.OptPos;
            optPos = (start->Pos + end->Pos) / 2.0;
            if (getQuadricCost(optPos, start->Q) + getQuadricCost(optPos, end->Q) > COST_THRESHOLD)
            {
                if (!Q.solve(optPos))
                {
                    // calculate the cost
                    const Vec3 &v0 = start->Pos;
                    const Vec3 &v1 = end->Pos;

                    double cost0 = getQuadricCost(start->Pos, Q);
                    double cost1 = getQuadricCost(end->Pos, Q);
                    double costm = getQuadricCost(optPos, Q);

                    double min = std::min(cost0, std::min(cost1, costm));
                    cost = min;
                    if (min == cost0)
                    {
                        optPos = start->Pos;
                    }
                    else if (min == cost1)
                    {
                        optPos = end->Pos;
                    }
                    // else use mid point
                    return;
                }
            }
            cost = getQuadricCost(optPos, Q);
            return;
        }

        static void solve(Face &face)
        {
            double min = std::numeric_limits<double>::max();
            for (unsigned int j = 0; j < 3; ++j)
            {
                Edge *e = face.Edges[j];
                // need to update if edge LOCAL_MARK is less than vertex LOCAL_MARK
                if (e->LOCAL_MARK <= e->Start->LOCAL_MARK || e->LOCAL_MARK <= e->End->LOCAL_MARK)
                {
                    e->LOCAL_MARK = GLOBAL_MARK;
                    calcOptimalPosition(*e, e->Priority);
                }
                if (e->Priority < min)
                {
                    face.OptEdge = j;
                    min = e->Priority;
                }
            }
            face.priority = face.Edges[face.OptEdge]->Priority;
        }

        static void solve(Edge &edge)
        {
            if (edge.Start->Pos == edge.End->Pos)
            {
                edge.OptPos = edge.Start->Pos;
                edge.Priority = -std::numeric_limits<double>::infinity();
            }

            double cost = 0.0;
            calcOptimalPosition(edge, cost);
            edge.Priority = computePriority(edge, cost);
        }

        static bool flipped(Vertex *start, Vertex *end, const Vec3 &optPos)
        {
            if (optPos == start->Pos)
            {
                return false;
            }
            for (auto node = start->Neighbors; node != nullptr; node = node->Next)
            {
                auto neighbor = node->Instance;
                if (!neighbor->Valid)
                {
                    continue;
                }
                if (neighbor->containVertex(end))
                {
                    continue;
                }

                Vertex *v1 = nullptr;
                Vertex *v2 = nullptr;
                neighbor->setCentral(start, v1, v2);
                const Vec3 &d1 = (v1->Pos - optPos).normalize();
                const Vec3 &d2 = (v2->Pos - optPos).normalize();

                if (std::fabs(d1.dot(d2)) > AREA_TOLERANCE)
                {
                    return true;
                }

                const Vec3 &unitNormal = (d1.cross(d2)).normalize();
                if (unitNormal.dot(neighbor->normal().normalize()) < NORMAL_TOLERANCE)
                {
                    return true;
                }
            }
            return false;
        }

        static unsigned int update(Face &face)
        {
            // get vertices
            Vertex *v0 = face.Vertices[face.OptEdge];
            Vertex *v1 = face.Vertices[(face.OptEdge + 1) % 3];

            const Vec3 &optPos = face.Edges[face.OptEdge]->OptPos;
            if (flipped(v0, v1, optPos) || flipped(v1, v0, optPos))
            {
                return 0;
            }

            // mark as Removed
            v1->Removed = true;

            // update v1 faces
            unsigned int nDeleted = 0;
            for (auto node = v1->Neighbors; node != nullptr; node = node->Next)
            {
                auto f = node->Instance;
                if (!f->Valid)
                {
                    continue;
                }

                // try to remove the face
                if (f->containVertex(v0))
                {
                    nDeleted++;
                    f->Valid = false;
                    f->replaceEdge(v0, v1);

                    auto v2 = f->getRestVertex(v0, v1);
                    // v2->Neighbors->clear();
                    removeInvalidFace(v2);
                    if(v2->Neighbors == nullptr)
                    {
                        v2->Removed = true;
                    }
                }
            }
            // v0->Neighbors = v0->Neighbors->import(v1);
            importVertexNeighbors(v0, v1);
            if (v0->Neighbors == nullptr)
            {
                v0->Removed = true;
                return nDeleted;
            }

            GLOBAL_MARK++;
            // use v0 to store new vertex
            v0->Pos = optPos;

            // update v0
            v0->Q += v1->Q;
            v0->LOCAL_MARK = GLOBAL_MARK;
            for (auto node = v1->Neighbors; node != nullptr; node = node->Next)
            {
                auto f = node->Instance;
                f->replace(v1, v0);
                f->replaceDuplicatedEdge();
            }

            // increase GLOBAL_MARK to update the faces
            GLOBAL_MARK++;
            for (auto node = v0->Neighbors; node != nullptr; node = node->Next)
            {
                auto f = node->Instance;

                // mark face as dirty
                f->LOCAL_MARK = GLOBAL_MARK;

                // update
                CollapseHelper::solve(*f);
            }
            return nDeleted;
        }

        static unsigned int updateEdge(EdgeHeap &heap, Edge &edge)
        {
            // edge cache
            static std::vector<Edge *> EdgeCache;

            // update global mark
            GLOBAL_MARK++;

            // get vertex
            Vertex *v0 = edge.Start;
            Vertex *v1 = edge.End;

            //prefer to use the one with larger capacity
            v1->Removed = true;

            // update v1 faces
            unsigned int nDeleted = 0;
            for (auto node = v1->Neighbors; node != nullptr; node = node->Next)
            {
                auto f = node->Instance;
                if (!f->Valid)
                {
                    continue;
                }
                // remove the face
                if (f->containVertex(v0))
                {
                    // set face as invalid
                    f->Valid = false;

                    // get another edge connected to v1, then remove it.
                    auto v2v1 = f->replaceEdge(v0, v1);
                    if (v2v1->HeapIndex != -2)
                    {
                        heap.pop(v2v1);
                    }

                    auto v2 = f->getRestVertex(v0, v1);
                    // v2->Neighbors->clear();
                    removeInvalidFace(v2);
                    if (!v2->Neighbors->Instance->Valid)
                    {
                        v2->Neighbors = v2->Neighbors->Next;
                    }

                    if (v2->Neighbors == nullptr) {
                        v2->Removed = true;
                    }
                    nDeleted++;
                    continue;
                }
            }

            // v0->Neighbors = v0->Neighbors->import(v1);
            importVertexNeighbors(v0, v1);
            if (v0->Neighbors == nullptr) {
                v0->Removed = true;
                return nDeleted;
            }

            // update v0 neighbors
            for (auto node = v1->Neighbors; node != nullptr; node = node->Next)
            {
                auto f = node->Instance;
                if (!(f)->Valid)
                {
                    continue;
                }

                // replace
                f->replace(v1, v0);
                f->replaceDuplicatedEdge();
            }

            // use v0 to store new vertex
            v0->Pos = edge.OptPos;
            v0->LOCAL_MARK = GLOBAL_MARK;
            v0->Q += v1->Q;

            Edge *e1 = nullptr;
            Edge *e2 = nullptr;
            EdgeCache.clear();
            for (auto node = v0->Neighbors; node != nullptr; node = node->Next)
            {
                auto f = node->Instance;
                f->getEdges(v0, e1, e2);
                syncIntoEdgeCache(e1, v0, EdgeCache);
                syncIntoEdgeCache(e2, v0, EdgeCache);
                if (e1->ToBeReplaced != nullptr || e2->ToBeReplaced != nullptr)
                {
                    f->replaceDuplicatedEdge();
                }
            }
            for (auto e : EdgeCache)
            {
                if (e->ToBeReplaced != nullptr)
                {
                    heap.pop(e);
                }
                else
                {
                    double old_priority = e->Priority;
                    solve(*e);
                    heap.update(e, old_priority);
                    e->Start->ToUpdateEdge = nullptr;
                }
            }
            return nDeleted;
        }

        static void release()
        {
            for (auto f : FacePool)
            {
                delete f;
            }
            for (auto v : VertexPool)
            {
                delete v;
            }
            for (auto e : EdgePool)
            {
                delete e;
            }
            FacePool.clear();
            VertexPool.clear();
            EdgePool.clear();
        }

        static void syncIntoEdgeCache(Edge *e1, Vertex *v0, std::vector<Edge *> &EdgeCache)
        {
            // force to be same topology
            if (e1->Start == v0)
            {
                std::swap(e1->Start, e1->End);
            }

            if (e1->LOCAL_MARK != GLOBAL_MARK)
            {
                if (e1->Start->ToUpdateEdge != nullptr && e1->Start->ToUpdateEdge != e1)
                {
                    e1->ToBeReplaced = e1->Start->ToUpdateEdge;
                }
                else
                {
                    e1->Start->ToUpdateEdge = e1;
                }
                e1->LOCAL_MARK = GLOBAL_MARK;
                EdgeCache.push_back(e1);
            }
        }

        static void removeInvalidFace(Vertex *v)
        {
            auto lastValid = v->Neighbors;
            for (auto node = lastValid->Next; node != nullptr; node = node->Next)
            {
                if (!node->Instance->Valid)
                {
                    continue;
                }
                lastValid->Next = node;
                lastValid = node;
            }
            lastValid->Next = nullptr;
        }

        static void importVertexNeighbors(Vertex *v0, Vertex *v1)
        {
            auto lastValid = v0->Neighbors;
            for (auto node = lastValid->Next; node != nullptr; node = node->Next)
            {
                if (!node->Instance->Valid) {
                    continue;
                }
                lastValid->Next = node;
                lastValid = node;
            }

            while (v1->Neighbors != nullptr && v1->Neighbors->Instance->Valid == false)
            {
                v1->Neighbors = v1->Neighbors->Next;
            }
            lastValid->Next = v1->Neighbors;
            if (v1->Neighbors != nullptr)
            {
                for (auto node = lastValid->Next; node != nullptr; node = node->Next)
                {
                    if (!node->Instance->Valid)
                    {
                        continue;
                    }
                    lastValid->Next = node;
                    lastValid = node;
                }
                lastValid->Next = nullptr;
            }

            if (v0->Neighbors->Instance->Valid == false)
            {
                v0->Neighbors = v0->Neighbors->Next;
            }
        }
    };

    // data section
    unsigned int CollapseHelper::GLOBAL_MARK = 0;
    double CollapseHelper::ScaleFactor = 1.0;

    // pool
    std::vector<Vertex *> CollapseHelper::VertexPool;
    std::vector<Face *> CollapseHelper::FacePool;
    std::vector<Edge *> CollapseHelper::EdgePool;
    uint32_t CollapseHelper::EdgePoolIdx = 0;

    class MeshReducerPrivate
    {
    public:
        MeshReducerPrivate();
        ~MeshReducerPrivate();

        void reset();
        void reduce(unsigned int nTarget, bool bForceStrict = false);

        void load(const double *vertices, const uint16_t *indices, unsigned int nVert,
                  unsigned int nInd);
        void store(std::vector<double> &vertices, std::vector<uint16_t> &indices);

        bool isNonClosedMesh() const;
        bool isValid() const;

    private:
        void buildQuadricMatrix();
        void initCollapses();

        void doFastLoop(unsigned int nTarget);
        void doStrictLoop(unsigned int nTarget);

        void cleanUp();
    private:
        // non-manifold
        bool m_bStrictConstraint = false;

        // data section
        std::vector<Vertex *> Vertices;
        std::vector<Face *> Faces;
        std::vector<Edge *> Edges;

        // Bounding box Diagonal
        double OriginBBoxDiagonal = 0.0;
    };

    MeshReducerPrivate::MeshReducerPrivate() {}

    MeshReducerPrivate::~MeshReducerPrivate()
    {
        CollapseHelper::release();
    }

    void MeshReducerPrivate::reset()
    {
        m_bStrictConstraint = false;
        Vertices.clear();
        Faces.clear();
        Edges.clear();
        CollapseHelper::reset();
    }

    void MeshReducerPrivate::reduce(unsigned int nTarget, bool bForceStrict) {
        m_bStrictConstraint = m_bStrictConstraint || bForceStrict;

        // build Quadric matrix for each vertex
        buildQuadricMatrix();

        // compute the new position and quadric error
        initCollapses();

        // loop
        if (bForceStrict || m_bStrictConstraint)
        {
            doStrictLoop(nTarget);
        }
        else
        {
            doFastLoop(nTarget);
        }

        // clean up
        cleanUp();
    }

    // TODO: CNECO-2636 Find out whether the index type should change to uint32.
    void MeshReducerPrivate::load(const double *vertices, const uint16_t *indices, unsigned int nVert,
                                  unsigned int nInd)
    {
        if (!vertices || !indices || !nVert || !nInd)
        {
            return;
        }

        // build vertices
        Vertices.resize(nVert);
        Vec3 min(std::numeric_limits<double>::max());
        Vec3 max(std::numeric_limits<double>::min());

        // reserve first
        unsigned int nFace = nInd / 3;
        CollapseHelper::reservePool(nVert, nFace);
        unsigned int vCounter = 0;
        for (unsigned int i_vert = 0; i_vert < nVert; ++i_vert) {
            double x = vertices[vCounter++];
            double y = vertices[vCounter++];
            double z = vertices[vCounter++];

            // get vertex from pool
            auto vert = CollapseHelper::VertexPool[i_vert];
            Vertices[i_vert] = vert;
            vert->Pos = Vec3(x, y, z);
            vert->ID = i_vert;

            // get scale
            min.X = std::min(min.X, x);
            max.X = std::max(max.X, x);

            min.Y = std::min(min.Y, y);
            max.Y = std::max(max.Y, y);

            min.Z = std::min(min.Z, z);
            max.Z = std::max(max.Z, z);
        }

        //CollapseHelper::ScaleFactor = double(1e8 * std::pow(1.0 / double((max - min).length()), 6));
        CollapseHelper::ScaleFactor = double(1e8 * std::pow(1.0 / double((max - min).squaredLength()), 3));
        OriginBBoxDiagonal = (max - min).squaredLength();

        // build faces
        unsigned int fCounter = 0;
        Faces.resize(nFace);
        Edges.reserve(2 * nFace + 1);
        for (unsigned int idx = 0; idx < nFace; ++idx)
        {
            // get face from pool
            auto face = CollapseHelper::FacePool[idx];
            Faces[idx] = face;
            // create vertex neighbor list
            for (unsigned int i_vert = 0; i_vert < 3; ++i_vert)
            {
                auto vert = Vertices[indices[fCounter++]];
                face->Vertices[i_vert] = vert;
                FaceNode* node = &face->FaceNodes[i_vert];
                node->Instance = face;
                node->Next = vert->Neighbors;
                vert->Neighbors = node;
            }

            for (unsigned int i_vert = 0; i_vert < 3; ++i_vert)
            {
                Vertex *v0 = face->Vertices[i_vert];
                Vertex *v1 = face->Vertices[(i_vert + 1) % 3];

                if (v0->ID > v1->ID)
                {
                    std::swap(v0, v1);
                }

                Edge* e = nullptr;
                // ToUpdateEdge and ToBeReplaced are used as Current & Next Edge
                // as a list structure
                Edge * possibleE = v0->CurrentEdge;
                while (possibleE != nullptr)
                {
                    if (possibleE->containVertex(v1))
                    {
                        e = possibleE;
                        break;
                    }
                    if (possibleE->NextEdge != nullptr)
                    {
                        possibleE = possibleE->NextEdge;
                    }
                    else
                    {
                        break;
                    }
                }

                if (e == nullptr)
                {
                    e = CollapseHelper::spawnEdgeFromPool(v0, v1);
                    Edges.push_back(e);
                    if (possibleE != nullptr)
                    {
                        e->NextEdge = v0->CurrentEdge;
                        v0->CurrentEdge = e;
                    }
                    else
                    {
                        v0->CurrentEdge = e;
                    }

                }
                e->AdjFaces++;
                face->Edges[i_vert] = e;
            }
        }
        for (auto v : Vertices)
        {
            v->CurrentEdge = nullptr;
        }
        for (auto e : Edges)
        {
            e->NextEdge = nullptr;
        }

        // manifold or non-manifold
        for (auto edge : Edges)
        {
            if (edge->AdjFaces == 1)
            {
                m_bStrictConstraint = true;
                break;
            }
        }
    }

    void MeshReducerPrivate::store(std::vector<double> &vertices, std::vector<uint16_t> &indices)
    {
        unsigned int nValid = 0;
        for (auto &vert : Vertices)
        {
            if (vert->Removed)
            {
                continue;
            }
            nValid++;
        }

        vertices.resize(nValid * 3);
        unsigned int i_valid = 0;
        unsigned int vCounter = 0;
        for (auto &vert : Vertices)
        {
            if (vert->Removed)
            {
                continue;
            }

            vertices[vCounter++] = vert->Pos.X;
            vertices[vCounter++] = vert->Pos.Y;
            vertices[vCounter++] = vert->Pos.Z;
            vert->ID = i_valid++;
        }

        indices.resize(Faces.size() * 3);
        int faceIndex = 0;
        for (auto face : Faces)
        {
            if (!face->Valid)
            {
                continue;
            }
            indices[faceIndex++] = Vertices[face->Vertices[0]->LOCAL_MARK]->ID;
            indices[faceIndex++] = Vertices[face->Vertices[1]->LOCAL_MARK]->ID;
            indices[faceIndex++] = Vertices[face->Vertices[2]->LOCAL_MARK]->ID;
        }
    }

    bool MeshReducerPrivate::isNonClosedMesh() const
    {
        return m_bStrictConstraint;
    }

    void MeshReducerPrivate::buildQuadricMatrix()
    {
        for (auto &face : Faces)
        {
            // ax + by + cz + d = 0
            Vec3 normal = face->normal();
            if (!m_bStrictConstraint)
            {
                normal = normal.normalize();
            }

            double d = -normal.dot(face->Vertices[0]->Pos);
            // assemble quadric matrix
            SymetricMatrix Q(normal.X, normal.Y, normal.Z, d);
            for (int i = 0; i < 3; ++i)
            {
                face->Vertices[i]->Q += Q;
            }
        }

        if (m_bStrictConstraint)
        {
            for (auto &face : Faces)
            {
                if (face->Edges[0]->AdjFaces != 1 && face->Edges[1]->AdjFaces != 1 && face->Edges[2]->AdjFaces != 1)
                {
                    continue;
                }

                Vec3 normal = face->normal();
                for (unsigned int j = 0; j < 3; ++j)
                {
                    if (face->Edges[j]->AdjFaces == 1)
                    {
                        Vertex *start = face->Vertices[j];
                        Vertex *end = face->Vertices[(j + 1) % 3];

                        const Vec3 &pStart = start->Pos;
                        const Vec3 &pEnd = end->Pos;

                        const Vec3 &edgePlane = normal.cross((pEnd - pStart).normalize());
                        double offset = -edgePlane.dot(pStart);
                        SymetricMatrix EQ(edgePlane.X, edgePlane.Y, edgePlane.Z, offset);

                        // add to related vertices
                        start->Q += EQ;
                        end->Q += EQ;
                    }
                }
            }
        }
    }

    void MeshReducerPrivate::initCollapses()
    {
        if (m_bStrictConstraint)
        {
            for (auto &edge : Edges)
            {
                CollapseHelper::solve(*edge);
            }
        }
        else
        {
            CollapseHelper::GLOBAL_MARK++;
            for (auto &face : Faces)
            {
                CollapseHelper::solve(*face);
            }
        }
    }

    void MeshReducerPrivate::doFastLoop(unsigned int nTarget)
    {
        unsigned int faceCount = static_cast < unsigned int> (Faces.size());
        unsigned int nowCount = faceCount;

        double minPriority = std::numeric_limits<double>::max();
        for (auto &face : Faces)
        {
            double priority = face->priority;
            if (priority < minPriority)
            {
                minPriority = priority;
            }
        }

        for (unsigned int iter = 0; iter < MAX_ITERATION; ++iter)
        {
            if (nowCount <= nTarget)
            {
                break;
            }

            double threshold = 1e-9 * static_cast<double>(std::pow(iter + 3, AGGRESSIVE));
            if (threshold <= minPriority)
            {
                continue;
            }

            int deletePerRound = 0;
            int currGlobalMark = CollapseHelper::GLOBAL_MARK;
            for (auto &face : Faces)
            {
                // get the top heap element
                if (!face->Valid || face->LOCAL_MARK > currGlobalMark || face->priority > threshold)
                {
                    continue;
                }

                // update
                unsigned int nDeleted = CollapseHelper::update(*face);
                deletePerRound += nDeleted;

                nowCount -= nDeleted;
                if (nowCount <= nTarget)
                {
                    break;
                }
            }

            if (nowCount > nTarget && deletePerRound > 0)
            {
                // update minimal priority
                minPriority = std::numeric_limits<double>::max();
                for (auto &face : Faces)
                {
                    if (face->Valid)
                    {
                        double d = face->priority;
                        if (d < minPriority)
                        {
                            minPriority = d;
                        }
                    }
                }
            }
        }
    }

    void MeshReducerPrivate::doStrictLoop(unsigned int nTarget)
    {
        unsigned int faceCount = static_cast<unsigned int>(Faces.size());
        unsigned int nowCount = faceCount;
        EdgeHeap heap(Edges);

        // clock
        auto start = std::chrono::steady_clock::now();
        auto end = start;

        int loopRound = 0;
        // collapse in loop
        while ((nowCount > nTarget) && heap.Length != 0 &&
               std::chrono::duration_cast<std::chrono::seconds>(end - start).count() < TIMEOUT)
        {
            // get top element
            Edge &top = *heap.top();

            if (top.Start->Removed || top.End->Removed)
            {
                continue;
            }

            // update time
            loopRound++;
            if (loopRound & 1)
            {
                end = std::chrono::steady_clock::now();
            }
            // update
            unsigned int nDeleted = CollapseHelper::updateEdge(heap, top);
            nowCount -= nDeleted;
        }
    }

    void MeshReducerPrivate::cleanUp()
    {
        int nValidVert = 0;
        for (auto &vert : Vertices)
        {
            if (vert->Removed)
            {
                continue;
            }

            Vertices[nValidVert++] = vert;
        }

        if (nValidVert == 0)
        {
            Vertices.clear();
            Faces.clear();
            return;
        }

        Vertices.resize(nValidVert);
        std::sort(Vertices.begin(), Vertices.end(), Vertex::comparefunc);

        // remove duplicated vertex
        unsigned int j = 0;
        unsigned int i = 1;
        unsigned int nDeleted = 0;
        for (; i != nValidVert;)
        {
            Vertices[i]->Removed = true;
            if (Vertices[i]->Pos == Vertices[j]->Pos)
            {
                Vertices[i]->LOCAL_MARK = j;
                i++;
                nDeleted++;
            }
            else
            {
                Vertices[i]->LOCAL_MARK = i;
                j = i;
                ++i;
            }
        }

        Vertices[0]->LOCAL_MARK = 0;
        Vertices[0]->Removed = true;

        // remove face out of range area
        // if share same vertex id, means face invalid
        unsigned int validFaceNum = 0;
        unsigned int faceCount = Faces.size();
        for (auto &face : Faces)
        {
            if (!face->Valid)
            {
                continue;
            }

            auto v0 = face->Vertices[0];
            auto v1 = face->Vertices[1];
            auto v2 = face->Vertices[2];
            if (v0->LOCAL_MARK == v1->LOCAL_MARK || v0->LOCAL_MARK == v2->LOCAL_MARK || v1->LOCAL_MARK == v2->LOCAL_MARK)
            {
                face->Valid = false;
                continue;
            }
            Faces[validFaceNum++] = face;
            Vertices[v0->LOCAL_MARK]->Removed = false;
            Vertices[v1->LOCAL_MARK]->Removed = false;
            Vertices[v2->LOCAL_MARK]->Removed = false;
        }
        Faces.resize(validFaceNum);
    }

    bool MeshReducerPrivate::isValid() const
    {
        if (Vertices.empty())
        {
            return false;
        }

        Vec3 min(std::numeric_limits<double>::max());
        Vec3 max(std::numeric_limits<double>::min());

        min.X = Vertices[0]->Pos.X;
        max.X = Vertices.back()->Pos.X;

        for (auto vert : Vertices)
        {
            if (vert->Removed)
            {
                continue;
            }

            min.Y = std::min(min.Y, vert->Pos.Y);
            max.Y = std::max(max.Y, vert->Pos.Y);

            min.Z = std::min(min.Z, vert->Pos.Z);
            max.Z = std::max(max.Z, vert->Pos.Z);
        }

        double len_diag = (max - min).squaredLength();
        double minDiag = std::min(len_diag, OriginBBoxDiagonal);
        double maxDiag = std::max(len_diag, OriginBBoxDiagonal);
        return minDiag > VALID_THRESHOLD * maxDiag;
    }
} // namespace

namespace common
{
    static MeshReducerPrivate reducer;
    bool MeshReducer::reduce(const double *vertices, const uint16_t *indices, unsigned int nVert,
                             unsigned int nIdx,
                             std::vector<double> &reducedVertices,
                             std::vector<uint16_t> &reducedIndices,
                             unsigned int nTarget)
    {
        if (vertices == nullptr || indices == nullptr || nVert == 0 || nIdx == 0 || nTarget == 0)
        {
            return false;
        }

        reducer.reset();
        reducer.load(vertices, indices, nVert, nIdx);
        reducer.reduce(nTarget);
        bool bValid = reducer.isValid();
        if (!bValid && !reducer.isNonClosedMesh())
        {
            // do again if the result is invalid and the mesh is non-manifold
            // force to go through the strict route
            reducer.reset();
            reducer.load(vertices, indices, nVert, nIdx);
            reducer.reduce(nTarget, true);
            bValid = reducer.isValid();
        }

        reducer.store(reducedVertices, reducedIndices);
        return bValid;
    }
} // namespace common