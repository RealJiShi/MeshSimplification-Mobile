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
static const unsigned int TIMEOUT = 5;
static const unsigned int MAX_ITERATION = 100;
static const double QUADRIC_EPSILON = 1e-15;
static const double VALID_THRESHOLD = 0.98;
static const double AGGRESSIVE = 8.0;
static const double AREA_TOLERANCE = 0.999;
static const double NORMAL_TOLERANCE = 0.2;

struct Vertex;
struct Face;
struct Edge;

struct Vertex
{
    Vertex() {}

    Vec3 Pos;
    unsigned int ID = 0;

    // attribute
    bool Removed = false;

    // Need to update?
    bool NeedToUpdate = false;

    // use LOCAL_MARK to update
    unsigned int LOCAL_MARK = 0;

    // inter variable
    Vec3 Pivot;
    double PivotSquare = 0.0;

    // Quadric Matrix
    SymetricMatrix Q;

    // adjecent faces
    std::list<Face *> Neighbors;
};

struct Edge
{
    Edge(Vertex *v0, Vertex *v1) : Start(v0), End(v1) {}

    // start & end vertex
    Vertex *Start = nullptr;
    Vertex *End = nullptr;

    // Number of adjecent faces
    unsigned int AdjFaces = 0;

    // collapse property for edge-based simplification
    unsigned int LOCAL_MARK = 0;
    double Priority = 0.0;
    Vec3 OptPos;

    // for little heap
    bool operator<(const Edge &rhs) const
    {
        return this->Priority > rhs.Priority;
    }
};

struct Face
{
    Face() {}
    Face(Vertex *v0, Vertex *v1, Vertex *v2)
    {
        Vertices[0] = v0;
        Vertices[1] = v1;
        Vertices[2] = v2;
    }

    // adjecent vertices
    Vertex *Vertices[3] = {nullptr};

    // edges
    unsigned int Edges[3] = {0};

    // valid & dirty(need to update)
    bool Valid = true;
    bool Dirty = false;

    // record v0 index
    unsigned int Central = 0;

    // collapse property for fast simplification
    Vec3 OptPos;
    unsigned int OptEdge = 0;
    double Priority = 0.0;

    void setCentral(Vertex *vert)
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            if (vert == Vertices[i])
            {
                Central = i;
                break;
            }
        }
    }

    void replace(Vertex *dst, Vertex *src)
    {
        for (unsigned int i = 0; i < 3; ++i)
        {
            if (dst == Vertices[i])
            {
                Vertices[i] = src;
                break;
            }
        }
    }

    Vertex *V0()
    {
        return Vertices[Central];
    }

    Vertex *V1()
    {
        return Vertices[(Central + 1) % 3];
    }

    Vertex *V2()
    {
        return Vertices[(Central + 2) % 3];
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

    double computeQuality() const
    {
        Vertex *v0 = Vertices[Central];
        Vertex *v1 = Vertices[(Central + 1) % 3];
        Vertex *v2 = Vertices[(Central + 2) % 3];

        const Vec3 &e0 = v2 - v1;
        const Vec3 &e1 = v1->Pivot;
        const Vec3 &e2 = v2->Pivot;

        Vec3 normal = e0.cross(e1);
        double len_norm = normal.length();
        if (len_norm == 0.0)
        {
            return 0;
        }

        double len_e0 = e0.squaredLength();
        double len_e1 = v1->PivotSquare;
        double len_e2 = v2->PivotSquare;
        double max_edge = std::max(len_e0, std::max(len_e1, len_e2));
        if (max_edge == 0.0)
        {
            return 0;
        }
        return len_norm / max_edge;
    }
};

class CollapseHelper
{
public:
    static unsigned int GLOBAL_MARK;
    static double ScaleFactor;

    static void reset()
    {
        GLOBAL_MARK = 0;
        ScaleFactor = 1.0;
    }

    // cost = VT * Q * V
    static double getQuadricCost(const Vec3 &v, const SymetricMatrix &Q)
    {
        double x = v.X;
        double y = v.Y;
        double z = v.Z;
        // Q[0] * x * x + 2 * Q[1] * x * y + 2 * Q[2] * x * z + 2 * Q[3] * x + Q[4] * y * y + 2 * Q[5] * y * z + 2 * Q[6] * y + Q[7] * z * z + 2 * Q[8] * z + Q[9];
        return (Q[0] * x + 2 * Q[1] * y + 2 * Q[2] * z + 2 * Q[3]) * x + (Q[4] * y + 2 * Q[5] * z + 2 * Q[6]) * y + (Q[7] * z + 2 * Q[8]) * z + Q[9];
    }

    // priority = cost / (normal * tri_quality)
    static double computePriority(const Vec3 &optPos, const double &QuadricCost, Vertex *start, Vertex *end)
    {
        // replace the position with optPos
        Vec3 start_old = start->Pos;
        Vec3 end_old = end->Pos;
        start->Pos = optPos;
        end->Pos = optPos;

        std::vector<Face *> FaceCollections;
        for (auto &face : start->Neighbors)
        {
            if (!face->Valid)
            {
                continue;
            }
            face->setCentral(start);
            if (face->V1() == end || face->V2() == end)
            {
                continue;
            }
            FaceCollections.push_back(face);
        }

        for (auto &face : end->Neighbors)
        {
            if (!face->Valid)
            {
                continue;
            }
            face->setCentral(end);
            if (face->V1() == start || face->V2() == start)
            {
                continue;
            }
            FaceCollections.push_back(face);
        }

        for (auto &face : FaceCollections)
        {
            if (!face->V1()->NeedToUpdate)
            {
                face->V1()->NeedToUpdate = true;
                face->V1()->Pivot = face->V1()->Pos - optPos;
                face->V1()->PivotSquare = face->V1()->Pivot.squaredLength();
            }
            if (!face->V2()->NeedToUpdate)
            {
                face->V2()->NeedToUpdate = true;
                face->V2()->Pivot = face->V2()->Pos - optPos;
                face->V2()->PivotSquare = face->V2()->Pivot.squaredLength();
            }
        }

        for (auto &face : FaceCollections)
        {
            for (auto &vert : face->Vertices)
            {
                vert->NeedToUpdate = false;
            }
        }

        // for each face related to start vertex
        double minQual = std::numeric_limits<double>::max();
        for (auto &face : FaceCollections)
        {
            double quality = face->computeQuality();
            minQual = std::min(minQual, quality);
        }

        // clamp
        minQual = std::min(minQual, 1.0);

        // cost
        double cost = ScaleFactor * QuadricCost;
        if (cost <= QUADRIC_EPSILON)
        {
            cost = -1 / (start_old - end_old).length();
        }
        cost /= minQual;

        // restore
        start->Pos = start_old;
        end->Pos = end_old;
        return cost;
    }

    static Vec3 calcOptimalPosition(const SymetricMatrix &Q, Vertex *start, Vertex *end, double &cost)
    {
        static const double COST_THRESHOLD = 200.0 * QUADRIC_EPSILON;
        Vec3 optPos = (start->Pos + end->Pos) / 2.0;
        if (getQuadricCost(optPos, start->Q) + getQuadricCost(optPos, end->Q) > COST_THRESHOLD)
        {
            if (!Q.solve(optPos))
            {
                // calculate the cost
                const Vec3 &v0 = start->Pos;
                const Vec3 &v1 = end->Pos;
                // const Vec3 &mid = (v0 + v1) / 2.0;
                const Vec3 &mid = optPos;

                double cost0 = getQuadricCost(v0, Q);
                double cost1 = getQuadricCost(v1, Q);
                double costm = getQuadricCost(mid, Q);

                double min = std::min(cost0, std::min(cost1, costm));
                cost = min;
                if (min == cost0)
                {
                    optPos = v0;
                }
                else if (min == cost1)
                {
                    optPos = v1;
                }
                else
                {
                    optPos = mid;
                }
            }
            else
            {
                cost = getQuadricCost(optPos, Q);
            }
        }
        return optPos;
    }

    static void solve(Face &face)
    {
        double min = std::numeric_limits<double>::max();
        for (unsigned int j = 0; j < 3; ++j)
        {
            Vertex *start = face.Vertices[j];
            Vertex *end = face.Vertices[(j + 1) % 3];
            SymetricMatrix Q = start->Q + end->Q;
            double cost = 0.0;
            Vec3 pos = calcOptimalPosition(Q, start, end, cost);
            if (cost < min)
            {
                face.OptPos = pos;
                face.Priority = cost;
                face.OptEdge = j;
                min = cost;
            }
        }
    }

    static void solve(Edge &edge)
    {
        SymetricMatrix Q = edge.Start->Q + edge.End->Q;
        double cost = 0.0;
        edge.OptPos = calcOptimalPosition(Q, edge.Start, edge.End, cost);
        edge.Priority = computePriority(edge.OptPos, cost, edge.Start, edge.End);
    }

    static bool flipped(Vertex *start, Vertex *end, const Vec3 &optPos)
    {
        if (optPos == start->Pos)
        {
            return false;
        }

        for (auto neighbor : start->Neighbors)
        {
            if (!neighbor->Valid)
            {
                continue;
            }

            neighbor->setCentral(start);
            if (neighbor->V1() == end || neighbor->V2() == end)
            {
                continue;
            }

            Vec3 d1 = (neighbor->V1()->Pos - optPos).normalize();
            Vec3 d2 = (neighbor->V2()->Pos - optPos).normalize();
            if (std::fabs(d1.dot(d2)) > AREA_TOLERANCE)
            {
                return true;
            }

            Vec3 unitNormal = (d1.cross(d2)).normalize();
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

        if (flipped(v0, v1, face.OptPos) || flipped(v1, v0, face.OptPos))
        {
            return 0;
        }

        // mark as Removed
        v1->Removed = true;

        // use v0 to store new vertex
        v0->Pos = face.OptPos;

        // update v0
        v0->Q += v1->Q;

        // update v1 faces
        unsigned int nDeleted = 0;
        for (auto iter = v1->Neighbors.begin(); iter != v1->Neighbors.end(); ++iter)
        {
            if (!(*iter)->Valid)
            {
                continue;
            }

            // try to remove the face
            if (v0 == (*iter)->Vertices[0] || v0 == (*iter)->Vertices[1] || v0 == (*iter)->Vertices[2])
            {
                nDeleted++;
                (*iter)->Valid = false;
                v1->Neighbors.erase(iter);
                continue;
            }

            // mark face as dirty
            (*iter)->Dirty = true;

            // replace
            (*iter)->replace(v1, v0);

            // update
            CollapseHelper::solve(**iter);

            // add to v0
            v0->Neighbors.push_back(*iter);
        }

        for (auto iter = v0->Neighbors.begin(); iter != v0->Neighbors.end(); ++iter)
        {
            if (!(*iter)->Valid)
            {
                continue;
            }

            // try to remove the face
            if (v1 == (*iter)->Vertices[0] || v1 == (*iter)->Vertices[1] || v1 == (*iter)->Vertices[2])
            {
                (*iter)->Valid = false;
                v0->Neighbors.erase(iter);
                continue;
            }

            // mark face as dirty
            (*iter)->Dirty = true;

            // update
            CollapseHelper::solve(**iter);
        }

        return nDeleted;
    }

    static unsigned int update(Edge &edge, std::vector<Edge> &Edges)
    {
        // update global mark
        GLOBAL_MARK++;

        // get vertex
        Vertex *v0 = edge.Start;
        Vertex *v1 = edge.End;
        v1->Removed = true;

        // use v0 to store new vertex
        v0->Pos = edge.OptPos;
        v0->LOCAL_MARK = GLOBAL_MARK;

        // update v1 faces
        unsigned int nDeleted = 0;
        for (auto iter = v1->Neighbors.begin(); iter != v1->Neighbors.end(); ++iter)
        {
            if (!(*iter)->Valid)
            {
                continue;
            }
            // remove the face
            if (v0 == (*iter)->Vertices[0] || v0 == (*iter)->Vertices[1] || v0 == (*iter)->Vertices[2])
            {
                // set face as invalid
                (*iter)->Valid = false;
                v1->Neighbors.erase(iter);
                nDeleted++;
                continue;
            }
            // replace
            (*iter)->replace(v1, v0);
            v0->Neighbors.push_back(*iter);
        }

        std::vector<Vertex *> VertexVec;
        for (auto &neighbor : v0->Neighbors)
        {
            if (!neighbor->Valid)
            {
                continue;
            }
            neighbor->setCentral(v0);
            if (!neighbor->V1()->NeedToUpdate)
            {
                neighbor->V1()->NeedToUpdate = true;
                VertexVec.push_back(neighbor->V1());
            }
            if (!neighbor->V2()->NeedToUpdate)
            {
                neighbor->V2()->NeedToUpdate = true;
                VertexVec.push_back(neighbor->V2());
            }
        }

        v0->Q += v1->Q;
        for (auto &v : VertexVec)
        {
            // reset first
            v->NeedToUpdate = false;
            Edges.emplace_back(v0, v);

            solve(Edges.back());
            // update mark
            Edges.back().LOCAL_MARK = GLOBAL_MARK;
            std::push_heap(Edges.begin(), Edges.end());
        }
        return nDeleted;
    }
};

unsigned int CollapseHelper::GLOBAL_MARK = 0;
double CollapseHelper::ScaleFactor = 1.0;

class MeshReducerPrivate
{
public:
    MeshReducerPrivate();
    ~MeshReducerPrivate();
    void reset();

    void reduce(unsigned int nTarget, bool bForceStrict = false);
    void load(const double *vertices, const uint16_t *indices, unsigned int nVert, unsigned int nInd);
    void store(std::vector<double> &vertices, std::vector<uint16_t> &indices);
    bool isManifoldMesh() const;
    bool isValid() const;

private:
    class DuplicateVertexCmp
    {
    public:
        inline bool operator()(Vertex *const &rhs, Vertex *const &lhs)
        {
            return (rhs->Pos == lhs->Pos ? (rhs < lhs) : (rhs->Pos < lhs->Pos));
        }
    };

    unsigned int buildEdge(Vertex *v0, Vertex *v1);
    void buildQuadricMatrix();
    void initCollapses();
    void doFastLoop(unsigned int nTarget);
    void doStrictLoop(unsigned int nTarget);
    void cleanUp();
    unsigned int removeFaceOutOfRangeArea();
    unsigned int removeDuplicateVertex();
    unsigned int removeUnreferenceVertex();

private:
    // non-manifold ratio
    bool m_bStrictConstraint = false;

    // data section
    std::vector<Vertex> Vertices;
    std::vector<Face> Faces;

    // build edge for border setup
    std::unordered_map<uint64_t, unsigned int> EdgeMap;
    std::vector<Edge> Edges;

    // Bounding box
    Vec3 OriginBBoxMin;
    Vec3 OriginBBoxMax;
};

MeshReducerPrivate::MeshReducerPrivate() {}

MeshReducerPrivate::~MeshReducerPrivate() {}

void MeshReducerPrivate::reset()
{
    m_bStrictConstraint = false;
    Vertices.clear();
    Faces.clear();
    EdgeMap.clear();
    Edges.clear();
    CollapseHelper::reset();
}

void MeshReducerPrivate::reduce(unsigned int nTarget, bool bForceStrict)
{
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
void MeshReducerPrivate::load(const double *vertices, const uint16_t *indices, unsigned int nVert, unsigned int nInd)
{
    if (!vertices || !indices || !nVert || !nInd)
    {
        return;
    }

    // build vertices
    Vertices.resize(nVert);
    Vec3 min(std::numeric_limits<double>::max());
    Vec3 max(std::numeric_limits<double>::min());

    for (unsigned int i_vert = 0; i_vert < nVert; ++i_vert)
    {
        double x = vertices[i_vert * 3 + 0];
        double y = vertices[i_vert * 3 + 1];
        double z = vertices[i_vert * 3 + 2];
        Vertices[i_vert].Pos = Vec3(x, y, z);
        Vertices[i_vert].ID = i_vert;

        // get scale
        min.X = std::min(min.X, x);
        max.X = std::max(max.X, x);

        min.Y = std::min(min.Y, y);
        max.Y = std::max(max.Y, y);

        min.Z = std::min(min.Z, z);
        max.Z = std::max(max.Z, z);
    }

    CollapseHelper::ScaleFactor = double(1e8 * std::pow(1.0 / double((max - min).length()), 6));
    OriginBBoxMax = max;
    OriginBBoxMin = min;

    // build faces
    Faces.resize(nInd / 3);
    for (unsigned int i_face = 0; i_face < nInd / 3; ++i_face)
    {
        auto &curr = Faces[i_face];
        for (unsigned int j = 0; j < 3; ++j)
        {
            curr.Vertices[j] = &Vertices[indices[i_face * 3 + j]];
            curr.Vertices[j]->Neighbors.push_back(&curr);
        }

        for (unsigned int j = 0; j < 3; ++j)
        {
            Vertex *v0 = curr.Vertices[j];
            Vertex *v1 = curr.Vertices[(j + 1) % 3];

            unsigned int i_edge = buildEdge(v0, v1);
            Edge &edge = Edges[i_edge];
            edge.AdjFaces++;
            curr.Edges[j] = i_edge;
        }
    }

    // manifold or non-manifold
    unsigned int nonManiEdge = 0;
    for (auto &edge : Edges)
    {
        if (edge.AdjFaces == 1)
        {
            nonManiEdge++;
        }
    }

    if (nonManiEdge > 0)
    {
        m_bStrictConstraint = true;
    }
}

void MeshReducerPrivate::store(std::vector<double> &vertices, std::vector<uint16_t> &indices)
{
    unsigned int i_valid = 0;
    for (unsigned int i = 0; i < Vertices.size(); ++i)
    {
        auto &vert = Vertices[i];
        if (vert.Removed)
        {
            continue;
        }

        vert.LOCAL_MARK = i_valid++;
        vertices.push_back(vert.Pos.X);
        vertices.push_back(vert.Pos.Y);
        vertices.push_back(vert.Pos.Z);
    }

    // faces
    for (auto &face : Faces)
    {
        if (!face.Valid)
        {
            continue;
        }

        for (auto &vert : face.Vertices)
        {
            indices.push_back((unsigned short)vert->LOCAL_MARK);
        }
    }
}

bool MeshReducerPrivate::isManifoldMesh() const
{
    return m_bStrictConstraint;
}

unsigned int MeshReducerPrivate::buildEdge(Vertex *v0, Vertex *v1)
{
    // find first
    uint64_t key = v0->ID < v1->ID ? (uint64_t(v0->ID) << 32) | v1->ID : (uint64_t(v1->ID) << 32) | v0->ID;
    if (EdgeMap.find(key) != EdgeMap.end())
    {
        return EdgeMap[key];
    }

    Edges.emplace_back(v0, v1);
    unsigned int index = Edges.size() - 1;
    EdgeMap[key] = index;
    return index;
}

void MeshReducerPrivate::buildQuadricMatrix()
{
    for (auto &face : Faces)
    {
        // ax + by + cz + d = 0
        Vec3 normal = face.normal();
        if (!m_bStrictConstraint)
        {
            normal = normal.normalize();
        }
        double d = -normal.dot(face.Vertices[0]->Pos);

        // assemble quadric matrix
        SymetricMatrix Q(normal.X, normal.Y, normal.Z, d);
        for (auto &vert : face.Vertices)
        {
            vert->Q += Q;
        }
    }

    if (m_bStrictConstraint)
    {
        for (auto &face : Faces)
        {
            Vec3 normal = face.normal();
            for (unsigned int j = 0; j < 3; ++j)
            {
                if (Edges[face.Edges[j]].AdjFaces == 1)
                {
                    const Vec3 &start = face.Vertices[j]->Pos;
                    const Vec3 &end = face.Vertices[(j + 1) % 3]->Pos;
                    const Vec3 &edgePlane = normal.cross((end - start).normalize());
                    double offset = -edgePlane.dot(start);
                    SymetricMatrix EQ(edgePlane.X, edgePlane.Y, edgePlane.Z, offset);

                    // add to related vertices
                    face.Vertices[j]->Q += EQ;
                    face.Vertices[(j + 1) % 3]->Q += EQ;
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
            CollapseHelper::solve(edge);
        }
    }
    else
    {
        for (auto &face : Faces)
        {
            CollapseHelper::solve(face);
        }
    }
}

void MeshReducerPrivate::doFastLoop(unsigned int nTarget)
{
    unsigned int faceCount = static_cast<unsigned int>(Faces.size());
    unsigned int nowCount = faceCount;
    for (unsigned int iter = 0; iter < MAX_ITERATION; ++iter)
    {
        if (nowCount <= nTarget)
        {
            break;
        }

        double threshold = 1e-9 * static_cast<double>(std::pow(iter + 3, AGGRESSIVE));
        for (auto &face : Faces)
        {
            // get the top heap element
            if (!face.Valid || face.Dirty || face.Priority > threshold)
            {
                continue;
            }

            // update
            unsigned int nDeleted = CollapseHelper::update(face);
            nowCount -= nDeleted;
            if (nowCount <= nTarget)
            {
                break;
            }
        }

        // clear dirty flag
        for (auto &face : Faces)
        {
            face.Dirty = false;
        }
    }
}

void MeshReducerPrivate::doStrictLoop(unsigned int nTarget)
{
    unsigned int faceCount = static_cast<unsigned int>(Faces.size());
    unsigned int nowCount = faceCount;

    // sort first
    std::make_heap(Edges.begin(), Edges.end());

    // clock
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();

    // collapse in loop
    while ((nowCount > nTarget) && !Edges.empty() &&
           std::chrono::duration_cast<std::chrono::seconds>(end - start).count() < TIMEOUT)
    {
        // update time
        end = std::chrono::steady_clock::now();

        Edge &top = Edges.front();
        std::pop_heap(Edges.begin(), Edges.end());
        Edges.pop_back();

        if (top.Start->Removed || top.End->Removed ||
            top.LOCAL_MARK < top.Start->LOCAL_MARK || top.LOCAL_MARK < top.End->LOCAL_MARK)
        {
            continue;
        }

        // update
        unsigned int nDeleted = CollapseHelper::update(top, Edges);
        nowCount -= nDeleted;
    }
}

void MeshReducerPrivate::cleanUp()
{
    removeFaceOutOfRangeArea();
    removeDuplicateVertex();
    removeUnreferenceVertex();
}

unsigned int MeshReducerPrivate::removeFaceOutOfRangeArea()
{
    unsigned int nDeleted = 0;
    for (auto &face : Faces)
    {
        if (!face.Valid)
        {
            continue;
        }

        const Vec3 &normal = face.normal();
        // double area = normal.length();
        // if (area == 0)
        if (normal.X == 0 && normal.Y == 0 && normal.Z == 0)
        {
            face.Valid = false;
            nDeleted++;
        }
    }
    return nDeleted;
}

unsigned int MeshReducerPrivate::removeDuplicateVertex()
{
    std::size_t nVert = Vertices.size();
    std::vector<Vertex *> pVec(nVert);
    unsigned int i_valid = 0;
    for (auto &vert : Vertices)
    {
        if (!vert.Removed)
        {
            pVec[i_valid] = &vert;
            i_valid++;
        }
    }
    pVec.resize(i_valid);

    DuplicateVertexCmp cmp;
    std::sort(pVec.begin(), pVec.end(), cmp);

    unsigned int j = 0;
    unsigned int i = 1;
    unsigned int nDeleted = 0;
    unsigned int INVALID_MARK = Vertices.size();
    for (auto &vert : pVec)
    {
        vert->LOCAL_MARK = vert->ID;
    }
    for (; i != i_valid;)
    {
        if (!pVec[i]->Removed && !pVec[j]->Removed && pVec[i]->Pos == pVec[j]->Pos)
        {
            pVec[i]->LOCAL_MARK = pVec[j]->ID;
            pVec[i]->Removed = true;
            i++;
            nDeleted++;
        }
        else
        {
            j = i;
            ++i;
        }
    }

    for (auto &face : Faces)
    {
        if (!face.Valid)
        {
            continue;
        }
        for (unsigned int k = 0; k < 3; ++k)
        {
            auto &vert = face.Vertices[k];
            if (vert->LOCAL_MARK != vert->ID)
            {
                face.Vertices[k] = &Vertices[vert->LOCAL_MARK];
            }
        }
    }
    return nDeleted;
}

unsigned int MeshReducerPrivate::removeUnreferenceVertex()
{
    std::size_t nVert = Vertices.size();
    std::vector<bool> bVec(nVert, false);

    // clear LOCAL_MARK
    for (auto &vert : Vertices)
    {
        vert.LOCAL_MARK = 0;
    }
    unsigned int VALID_MARK = 1;
    for (auto &face : Faces)
    {
        if (!face.Valid)
        {
            continue;
        }
        for (auto &vert : face.Vertices)
        {
            vert->LOCAL_MARK = VALID_MARK;
        }
    }

    unsigned int nDeleted = 0;
    for (auto &vert : Vertices)
    {
        if (vert.Removed)
        {
            continue;
        }
        if (vert.LOCAL_MARK ^ VALID_MARK)
        {
            vert.Removed = true;
            nDeleted++;
        }
    }
    return nDeleted;
}

bool MeshReducerPrivate::isValid() const
{
    Vec3 min(std::numeric_limits<double>::max());
    Vec3 max(std::numeric_limits<double>::min());
    for (auto vert : Vertices)
    {
        if (vert.Removed)
        {
            continue;
        }

        // get scale
        min.X = std::min(min.X, vert.Pos.X);
        max.X = std::max(max.X, vert.Pos.X);

        min.Y = std::min(min.Y, vert.Pos.Y);
        max.Y = std::max(max.Y, vert.Pos.Y);

        min.Z = std::min(min.Z, vert.Pos.Z);
        max.Z = std::max(max.Z, vert.Pos.Z);
    }

    double len_diag = (max - min).length();
    double len_diag_old = (OriginBBoxMax - OriginBBoxMin).length();
    double ratio = std::min(len_diag, len_diag_old) / std::max(len_diag, len_diag_old);
    return ratio > VALID_THRESHOLD;
}
} // namespace

namespace common
{

static MeshReducerPrivate reducer;
bool MeshReducer::reduce(const double *vertices, const uint16_t *indices, unsigned int nVert, unsigned int nIdx,
                         std::vector<double> &reducedVertices, std::vector<uint16_t> &reducedIndices,
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
    if (!bValid && !reducer.isManifoldMesh())
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
