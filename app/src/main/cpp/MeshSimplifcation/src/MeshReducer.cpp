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

namespace {
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

    struct Vertex {
        Vertex() {}

        Vec3 Pos;
        uint16_t ID = 0;

        // attribute
        bool Removed = false;

        // Need to update?
        Edge *toUpdateEdge = nullptr;
        //bool NeedToUpdate = false;

        // use LOCAL_MARK to update
        unsigned int LOCAL_MARK = 0;

        // Quadric Matrix
        SymetricMatrix Q;

        // adjacent faces
        std::vector<Face *> Neighbors;

        void reset() {
            Neighbors.clear();
            LOCAL_MARK = 0;
            toUpdateEdge = nullptr;
            ID = 0;
            Removed = false;
            Q.reset();
            //NeedToUpdate = false;
        }
    };

    struct Edge {
        Edge(Vertex *v0, Vertex *v1) : Start(v0), End(v1) {}

        // start & end vertex
        Vertex *Start = nullptr;
        Vertex *End = nullptr;

        // Number of adjacent faces
        unsigned int AdjFaces = 0;

        // collapse property for edge-based simplification
        unsigned int LOCAL_MARK = 0;
        double Priority = 0.0;
        Vec3 OptPos;

        void replaceVertex(Vertex *dst, Vertex *src) {
            if (Start == dst) {
                Start = src;
            } else if (End == dst) {
                End = src;
            }
        }

        bool containVertex(Vertex *v1) {
            return (Start == v1) || (End == v1);
        }

        bool containVertex(Vertex *v1, Vertex *v2) {
            return (Start == v1 && End == v2) || (Start == v2 && End == v1);
        }


        void updateEdge(Vertex *v0, Vertex *v1) {
            Start = v0;
            End = v1;

            // Number of adjacent faces
            AdjFaces = 0;

            // collapse property for edge-based simplification
            LOCAL_MARK = 0;
            Priority = 0.0;
            HeapIndex = -1;
            toBeReplaced = nullptr;
        }


        int HeapIndex = -1;
        Edge *toBeReplaced = nullptr;;

        bool sameEndWith(Edge *b) {
            if (Start == b->Start && End == b->End) return true;
            if (Start == b->End && End == b->Start) return true;
            return false;
        }

        uint64_t getKey() {
            auto v0 = Start;
            auto v1 = End;
            return v0->ID < v1->ID ? (uint64_t(v0->ID) << 32) | v1->ID : (uint64_t(v1->ID) << 32) |
                                                                         v0->ID;
        }

        bool smallerPriority(Edge *p) {
            return this->Priority < p->Priority;
        }

    };


    struct Face {
        Face() {}

        void reset() {
            Vertices[0] = Vertices[1] = Vertices[2] = nullptr;
            ConsistentEdges[0] = ConsistentEdges[1] = ConsistentEdges[2] = nullptr;

            Valid = true;
            LOCAL_MARK = 0;

            OptEdge = 0;
        }

        // adjacent vertices
        Vertex *Vertices[3] = {nullptr};

        // edges
        Edge *ConsistentEdges[3] = {0};

        // valid & dirty(need to update)
        bool Valid = true;

        // collapse property for fast simplification
        uint16_t OptEdge = 0;
        unsigned int LOCAL_MARK = 0; //to check whether face are updated for current iteration

        double priority = 0;

        void GetEdgesByVertex(Vertex *v0, Edge *&e1, Edge *&e2) {
            if (v0 == Vertices[0]) {
                e1 = ConsistentEdges[0];
                e2 = ConsistentEdges[2];
            } else if (v0 == Vertices[1]) {
                e1 = ConsistentEdges[0];
                e2 = ConsistentEdges[1];
            } else //if (v0 == Vertices[2])
            {
                e1 = ConsistentEdges[1];
                e2 = ConsistentEdges[2];
            }
        }

        void setCentral(Vertex *v0, Vertex *&v1, Vertex *&v2) const {
            if (v0 == Vertices[0]) {
                v1 = Vertices[1];
                v2 = Vertices[2];
            } else if (v0 == Vertices[1]) {
                v1 = Vertices[2];
                v2 = Vertices[0];
            } else //if (v0 == Vertices[2])
            {
                v1 = Vertices[0];
                v2 = Vertices[1];
            }
        }


        void replaceKeepTopology(Vertex *dst, Vertex *src) {
            if (dst == Vertices[0]) {
                Vertices[0] = src;
                ConsistentEdges[0]->replaceVertex(dst, src);
                ConsistentEdges[2]->replaceVertex(dst, src);
            } else if (dst == Vertices[1]) {
                Vertices[1] = src;
                ConsistentEdges[0]->replaceVertex(dst, src);
                ConsistentEdges[1]->replaceVertex(dst, src);
            } else if (dst == Vertices[2]) {
                Vertices[2] = src;
                ConsistentEdges[1]->replaceVertex(dst, src);
                ConsistentEdges[2]->replaceVertex(dst, src);
            }
        }

        bool containVertex(Vertex *p) {
            return Vertices[0] == p || Vertices[1] == p || Vertices[2] == p;
        }

        Vec3 normal() const {
            const Vec3 &v0 = Vertices[0]->Pos;
            const Vec3 &v1 = Vertices[1]->Pos;
            const Vec3 &v2 = Vertices[2]->Pos;

            const Vec3 &e0 = v1 - v0;
            const Vec3 &e1 = v2 - v0;

            return e0.cross(e1);
        }

        double computeQuality(Vertex *start, const Vec3 &optPos) {
            Vertex *v1 = nullptr;
            Vertex *v2 = nullptr;
            this->setCentral(start, v1, v2);

            const Vec3 &e0 = v2->Pos - v1->Pos;
            const Vec3 &e1 = v1->Pos - optPos;

            const Vec3 &normal = e0.cross(e1);

            double len_norm = normal.length();
            if (len_norm == 0.0) {
                return 0; //the opts pos on the same line
            }

            const Vec3 &e2 = v2->Pos - optPos;
            double len_e0 = e0.squaredLength();
            double len_e1 = e1.squaredLength();
            double len_e2 = e2.squaredLength();
            double max_edge = std::max(len_e0, std::max(len_e1, len_e2));

            return len_norm / max_edge;

        }

        void getEdgesByVertex(Vertex *v0, Vertex *v1, Edge *&v2v1, Edge *&v2v0) {
            if (Vertices[0] == v0) {
                if (Vertices[1] == v1)  //v0 0 v1 1 v2 2
                {
                    v2v1 = this->ConsistentEdges[1];
                    v2v0 = this->ConsistentEdges[2];
                } else //v1 -> 2        //v0 0  v2 1 v1 2
                {
                    v2v0 = this->ConsistentEdges[0];
                    v2v1 = this->ConsistentEdges[1];
                }
            } else if (Vertices[1] == v0) {
                if (Vertices[2] == v1)  //v2 0  v0 1 v1 2
                {
                    v2v1 = this->ConsistentEdges[2];
                    v2v0 = this->ConsistentEdges[0];
                } else  //v1-> 0			//v1 0  v0 1 v2 2
                {
                    v2v1 = this->ConsistentEdges[2];
                    v2v0 = this->ConsistentEdges[1];
                }
            } else {// if (Vertices[2] == v0) {
                if (Vertices[0] == v1)  //v1 0  v2 1 v0 2
                {
                    v2v1 = this->ConsistentEdges[0];
                    v2v0 = this->ConsistentEdges[1];
                } else                    //v2 0  v1 1 v0 2
                {
                    v2v1 = this->ConsistentEdges[0];
                    v2v0 = this->ConsistentEdges[2];
                }
            }
        }

        void replaceDuplicateEdge() {
            if (ConsistentEdges[0]->toBeReplaced != nullptr) {
                ConsistentEdges[0] = ConsistentEdges[0]->toBeReplaced;
            }
            if (ConsistentEdges[1]->toBeReplaced != nullptr) {
                ConsistentEdges[1] = ConsistentEdges[1]->toBeReplaced;
            }
            if (ConsistentEdges[2]->toBeReplaced != nullptr) {
                ConsistentEdges[2] = ConsistentEdges[2]->toBeReplaced;
            }
        }


        void markV2V1ReplacedByV2V0(Vertex *v0, Vertex *v1) {
            Edge *v2v1 = nullptr;
            Edge *v2v0 = nullptr;
            this->getEdgesByVertex(v0, v1, v2v1, v2v0);
            v2v1->toBeReplaced = v2v0;
            //return v2v1;
        }

        Edge *GetV2V1ReplacedByV2V0(Vertex *v0, Vertex *v1) {
            Edge *v2v1 = nullptr;
            Edge *v2v0 = nullptr;
            this->getEdgesByVertex(v0, v1, v2v1, v2v0);
            v2v1->toBeReplaced = v2v0;
            return v2v1;
        }
    };

    struct MinEdgeHeap {
        Edge **arr;
        int size;

        MinEdgeHeap(std::vector<Edge *> &edges) {
            arr = &edges[0];
            size = edges.size();

            //how to get
            for (int idx = 1; idx < size; ++idx) {
                heapifyUpInit(idx);
            }
            for (int i = size - 1; i >= 0; --i) {
                arr[i]->HeapIndex = i;
            }
        }

        Edge *extractMin() {
            Edge *item = arr[0];
            item->HeapIndex = -1; //mark as not inside heap anymore

            //Remove the last element
            size--;
            if (size != 0) //current item is not the previous last item. Swap needed
            {
                //This messes up the index field
                arr[0] = arr[size]; //get previous last element
                //So reset it
                arr[0]->HeapIndex = 0;

                //////Heapify from the position we just messed with
                heapifyDown(arr[0]);
            }
            return item;
        }

        void updateHeapEdgePriority(Edge *e, double priority) {
            if (arr[e->HeapIndex] != e) {
                std::cout << "ERROR>  heapifyDown  MUST FIX............";
                return;
            }

            if (priority < e->Priority) {
                e->Priority = priority;
                heapifyUp(e);
            } else if (priority > e->Priority) {
                e->Priority = priority;
                heapifyDown(e);
            }
        }

        //No performance gain if we do the delete edge..
        void deleteEdge(Edge *deletedEdge) {
            int pos = deletedEdge->HeapIndex;
            if (pos >= 0 && pos < size && deletedEdge == arr[pos]) {
            } else {
                if (deletedEdge->HeapIndex !=
                    -2) //already deleted, as same face with diff orientation exist
                {
                    std::cout << "ERROR...DELETE EDGE should not be here" << std::endl;
                }
                return;
            }
            deletedEdge->HeapIndex = -2;//-2 means deleteEdge
            //Remove the last element
            size--;
            if (pos != size)//if pos is not the previous last element, swap needed
            {
                //This messes up the index field
                auto e = arr[size];// get the previous last element
                //So reset it
                arr[pos] = e;
                e->HeapIndex = pos;

                if (e->Priority < deletedEdge->Priority) {
                    heapifyUp(e);
                } else if (e->Priority > deletedEdge->Priority) {
                    heapifyDown(e);
                }
            }
        }

    private:// function to heapify the tree

        //heapifyUpInit expect pos >= 1
        void heapifyUpInit(int pos) {
            int parent = (pos - 1) >> 1;// / 2;

            auto e = arr[pos];
            while (e->Priority < arr[parent]->Priority) {
                // swap and keep indexes updated
                arr[pos] = arr[parent];
                //heap[parent] = e;

                if (parent == 0) {
                    arr[0] = e;
                    return;
                }

                pos = parent;
                parent = (pos - 1) >> 1;
            }

            arr[pos] = e;
            return;
        }

        void heapifyUp(Edge *e) {
            if (e->HeapIndex == 0) return;
            int pos = e->HeapIndex;

            int parent = (pos - 1) >> 1;// / 2;
            while (e->smallerPriority(arr[parent])) {
                // swap and keep indexes updated
                arr[pos] = arr[parent];
                arr[pos]->HeapIndex = pos; //update latest pos element's index
                //heap[parent] = e;

                if (parent == 0) //already at the top of heap
                {
                    e->HeapIndex = 0;
                    arr[0] = e;
                    return;
                }

                pos = parent;
                parent = (pos - 1) >> 1;
            }

            //pos has been updated. need to sync now
            arr[pos] = e;
            e->HeapIndex = pos; //update target edge's index
            return;
        }

        void heapifyDown(Edge *e) {
            int pos = e->HeapIndex;

            int l = (pos << 1) + 1;// 2 * pos + 1;
            if (l >= size) return;//already the leaf node. no need to do anything

            int r = l + 1;// 2 * pos + 2;
            int smallest = pos; // root is the Smallest element
            // If left child is smaller than root
            while (true) {

                //if (l < n && heap[l]->Priority < e->Priority)
                if (arr[l]->smallerPriority(e))
                    smallest = l;

                // If right child is smaller than largest so far
                if (r < size && arr[r]->smallerPriority(arr[smallest]))
                    smallest = r;

                if (smallest == pos) //already the smallest
                {
                    e->HeapIndex = pos; //no need to swap down
                    return;
                }

                // swap and keep indexes updated
                arr[pos] = arr[smallest];
                arr[pos]->HeapIndex = pos; //update latest pos element's index

                arr[smallest] = e; //heap[smallest] might be used

                l = (smallest << 1) + 1;// 2 * pos + 1;
                if (l >= size) {
                    e->HeapIndex = smallest; //no need to swap down, currently e located at smallest
                    return;//already the leaf node. no need to do anything
                }

                pos = smallest;
                r = l + 1;// (pos << 1) + 2;// 2 * pos + 2;
            }
        }
    };

    class CollapseHelper {
    public:
        static unsigned int GLOBAL_MARK;
        static double ScaleFactor;
        static SymetricMatrix tempQ;


        static std::vector<Face *> FacesPool;
        static size_t facePoolIdx;

        static std::vector<Edge *> EdgesPool;
        static size_t edgePoolIdx;
        static std::vector<Vertex *> VertexPool;
        static size_t vertexPoolIdx;

        static void reset() {
            GLOBAL_MARK = 0;
            ScaleFactor = 1.0;
            edgePoolIdx = 0;
            vertexPoolIdx = 0;
            facePoolIdx = 0;
        }

        static Vertex *spawnVertexFromPool() {
            if (vertexPoolIdx >= VertexPool.size()) {
                Vertex *v = new Vertex();
                VertexPool.push_back(v);
                vertexPoolIdx++;
                return v;
            }

            Vertex *v = VertexPool[vertexPoolIdx++];
            v->reset();
            return v;
        }

        static Face *spawnFaceFromPool() {
            if (facePoolIdx >= FacesPool.size()) {
                FacesPool.push_back(new Face());
                return FacesPool[facePoolIdx++];
            }
            auto f = FacesPool[facePoolIdx++];
            f->reset();
            return f;
        }

        static void reserveFacePool(int nFace) {
            if (FacesPool.capacity() < nFace)
                FacesPool.reserve(nFace);
        }

        static Edge *spawnEdgeFromPool(Vertex *v0, Vertex *v1) {
            if (EdgesPool.size() > edgePoolIdx) {
                Edge *e = EdgesPool[edgePoolIdx++];
                e->updateEdge(v0, v1);
                return e;
            }
            auto e = new Edge(v0, v1);
            EdgesPool.push_back(e);
            edgePoolIdx++;
            return e;
        }

        // cost = VT * Q * V
        static double getQuadricCost(const Vec3 &v, const SymetricMatrix &Q) {
            double x = v.X;
            double y = v.Y;
            double z = v.Z;
            // Q[0] * x * x + 2 * Q[1] * x * y + 2 * Q[2] * x * z + 2 * Q[3] * x + Q[4] * y * y + 2 * Q[5] * y * z + 2 * Q[6] * y + Q[7] * z * z + 2 * Q[8] * z + Q[9];
            return (Q[0] * x + 2 * Q[1] * y + 2 * Q[2] * z + 2 * Q[3]) * x +
                   (Q[4] * y + 2 * Q[5] * z + 2 * Q[6]) * y + (Q[7] * z + 2 * Q[8]) * z + Q[9];
        }

        // priority = cost / (normal * tri_quality)
        static double computePriority(Edge &e, const double &QuadricCost) {
            const Vec3 &optPos = e.OptPos;
            auto start = e.Start;
            auto end = e.End;

            // for each face related to start vertex
            double minQual = 1;
            for (auto face : start->Neighbors) {
                if (!face->Valid || face->containVertex(end)) {
                    continue;
                }

                double quality = face->computeQuality(start, optPos);
                if (quality == 0) {
                    minQual = 0;
                    break;
                }
                minQual = std::min(minQual, quality);
            }

            // for each face related to end vertex
            if (minQual != 0)
                for (auto face : end->Neighbors) {
                    if (!face->Valid || face->containVertex(start)) {
                        continue;
                    }

                    double quality = face->computeQuality(end, optPos);
                    if (quality == 0) {
                        minQual = 0;
                        break;
                    }
                    minQual = std::min(minQual, quality);
                }

            if (minQual == 0) {
                //avoid selection of edge which would cause regression of faces
                return std::numeric_limits<double>::infinity(); //- std::numeric_limits < double >::max();
            }

            // cost
            double cost = ScaleFactor * QuadricCost;
            if (cost <= QUADRIC_EPSILON) {
                cost = -1 / (start->Pos - end->Pos).length();
            }
            cost /= minQual;
            return cost;
        }

        static void calcOptimalPosition(Edge &e, double &cost) {

            static const double COST_THRESHOLD = 200.0 * QUADRIC_EPSILON;
            auto start = e.Start;
            auto end = e.End;

            const SymetricMatrix &Q = start->Q + end->Q;
            Vec3 &optPos = e.OptPos;
            optPos = (start->Pos + end->Pos) / 2.0;
            if (getQuadricCost(optPos, start->Q) + getQuadricCost(optPos, end->Q) >
                COST_THRESHOLD) {
                if (!Q.solve(optPos)) {
                    // calculate the cost
                    const Vec3 &v0 = start->Pos;
                    const Vec3 &v1 = end->Pos;

                    double cost0 = getQuadricCost(start->Pos, Q);
                    double cost1 = getQuadricCost(end->Pos, Q);
                    double costm = getQuadricCost(optPos, Q);

                    double min = std::min(cost0, std::min(cost1, costm));
                    cost = min;
                    if (min == cost0) {
                        optPos = start->Pos;
                    } else if (min == cost1) {
                        optPos = end->Pos;
                    } else {
                    }
                    return;
                }
            }
            cost = getQuadricCost(optPos, Q);
            return;

        }

        static void solveFaceOptimized(Face &face) {
            double min = std::numeric_limits<double>::max();
            for (unsigned int j = 0; j < 3; ++j) {
                Edge *e = face.ConsistentEdges[j];
                if (e->LOCAL_MARK < GLOBAL_MARK) {
                    if (e->LOCAL_MARK <= e->Start->LOCAL_MARK ||
                        e->LOCAL_MARK <= e->End->LOCAL_MARK) {
                        e->LOCAL_MARK = GLOBAL_MARK;
                        calcOptimalPosition(*e, e->Priority);
                    }
                }
                if (e->Priority < min) {
                    face.OptEdge = j;
                    min = e->Priority;
                }
            }
            face.priority = face.ConsistentEdges[face.OptEdge]->Priority;
        }

        static void solve(Edge &edge) {
            edge.Priority = solveEdgePriority(edge);

        }

        inline static double solveEdgePriority(Edge &edge) {
            if (edge.Start->Pos == edge.End->Pos) {
                edge.OptPos = edge.Start->Pos;
                return -std::numeric_limits<double>::infinity();
            }

            double cost = 0.0;
            calcOptimalPosition(edge, cost);
            return computePriority(edge, cost);
        }

        static bool flipped(Vertex *start, Vertex *end, const Vec3 &optPos) {
            if (optPos == start->Pos) {
                return false;
            }
            for (auto neighbor : start->Neighbors) {
                if (!neighbor->Valid) {
                    continue;
                }
                if (neighbor->containVertex(end)) {
                    continue;
                }

                Vertex *v1 = nullptr;
                Vertex *v2 = nullptr;
                neighbor->setCentral(start, v1, v2);
                const Vec3 &d1 = (v1->Pos - optPos).normalize();
                const Vec3 &d2 = (v2->Pos - optPos).normalize();

                if (std::fabs(d1.dot(d2)) > AREA_TOLERANCE) {
                    return true;
                }

                const Vec3 &unitNormal = (d1.cross(d2)).normalize();
                if (unitNormal.dot(neighbor->normal().normalize()) < NORMAL_TOLERANCE) {
                    return true;
                }
            }
            return false;
        }

        static int SwapToBack(std::vector<Face *> &Neighbors, int idx, int v0Size) {
            int bIdx = v0Size - 1;
            for (; bIdx > idx; --bIdx) {
                if (Neighbors[bIdx]->Valid) {
                    Neighbors[idx] = Neighbors[bIdx];
                    return bIdx;
                }
            }
            return bIdx;
        }

        static void UnorderedSwapInvalidFace(std::vector<Face *> &InputFaces) {

            if (InputFaces.size() == 0) return;
            int frontIdx = 0;
            int backIdx = InputFaces.size() - 1;
            for (; frontIdx <= backIdx; ++frontIdx) {
                auto e = InputFaces[frontIdx];
                if (e->Valid == false) {
                    bool swap = false;
                    for (; backIdx > frontIdx; --backIdx) {
                        if (InputFaces[backIdx]->Valid) {
                            InputFaces[frontIdx] = InputFaces[backIdx];
                            InputFaces[backIdx] = e;
                            backIdx--;
                            swap = true;
                            break;
                        }
                    }
                    if (swap == false) break;
                }
            }
            frontIdx = std::min(frontIdx, backIdx);
            if (InputFaces[frontIdx]->Valid) {

                InputFaces.resize(frontIdx + 1);
            } else
                InputFaces.resize(frontIdx);

        }

        static unsigned int update(Face &face) {

            // get vertices
            Vertex *v0 = face.Vertices[face.OptEdge];
            Vertex *v1 = face.Vertices[(face.OptEdge + 1) % 3];

            Vec3 &optPos = face.ConsistentEdges[face.OptEdge]->OptPos;
            face.Valid = false; //assume the remove of this face would be successfully

            if (flipped(v0, v1, optPos) || flipped(v1, v0, optPos)) {
                face.Valid = true; //due to flip, the face remove is reverted..
                // this would means the face would be flipped. no point to try flip again
                //face.LOCAL_MARK = GLOBAL_MARK + 100000;
                return 0;
            }
            //prefer to use the one with larger capacity
            if (v0->Neighbors.capacity() < v1->Neighbors.capacity()) {
                //swap
                auto temp = v0;
                v0 = v1;
                v1 = temp;
            }

            GLOBAL_MARK++;
            // mark as Removed
            v1->Removed = true;

            // use v0 to store new vertex
            v0->Pos = optPos;

            // update v0
            v0->Q += v1->Q;
            v0->LOCAL_MARK = GLOBAL_MARK;

            // update v1 faces
            unsigned int nDeleted = 1; //surely face would be invalid. already set valid to false.
            face.markV2V1ReplacedByV2V0(v0, v1);

            for (auto f : v1->Neighbors) {
                if (!f->Valid) {
                    continue;
                }
                // try to remove the face
                if (f->containVertex(v0)) {
                    nDeleted++;
                    f->Valid = false;
                    f->markV2V1ReplacedByV2V0(v0, v1);
                }
            }

            UnorderedSwapInvalidFace(v0->Neighbors);
            for (auto f : v1->Neighbors) {
                if (!f->Valid) {
                    continue;
                }

                f->replaceKeepTopology(v1, v0);
                f->replaceDuplicateEdge();

                // add to v0
                v0->Neighbors.push_back(f);
            }
            v1->Neighbors.clear();

            GLOBAL_MARK++;
            //it has been validated that cache would take effect.
            //sometime there would be same triangle with diff index
            //which would make the num of cal edge work < v0->neighbors.size();
            for (auto f : v0->Neighbors) {
                // mark face as dirty
                f->LOCAL_MARK = GLOBAL_MARK;

                // update
                CollapseHelper::solveFaceOptimized(*f);
            }
            return nDeleted;
        }

        static unsigned int
        updateEdge(MinEdgeHeap &heap, Edge &edge, std::vector<Edge *> &EdgeCache) {
            // update global mark
            GLOBAL_MARK++;

            // get vertex
            Vertex *v0 = edge.Start;
            Vertex *v1 = edge.End;

            //prefer to use the one with larger capacity
            if (v0->Neighbors.capacity() < v1->Neighbors.capacity()) {
                //swap
                auto temp = v0;
                v0 = v1;
                v1 = temp;
            }
            v1->Removed = true;
            // update v1 faces
            unsigned int nDeleted = 0;
            //v1->Neighbors.erase(iter);
            for (auto f : v1->Neighbors) {
                if (!f->Valid) {
                    continue;
                }
                // remove the face
                if (f->containVertex(v0)) {
                    // set face as invalid
                    f->Valid = false;

                    auto v2v1 = f->GetV2V1ReplacedByV2V0(v0, v1);
                    if (v2v1->HeapIndex != -2) {
                        heap.deleteEdge(v2v1); //if not deleted
                    }
                    nDeleted++;
                    continue;
                }
            }
            UnorderedSwapInvalidFace(v0->Neighbors);
            for (auto f : v1->Neighbors) {
                if (!(f)->Valid) {
                    continue;
                }

                // replace
                f->replaceKeepTopology(v1, v0);
                f->replaceDuplicateEdge();
                v0->Neighbors.push_back(f);

            }
            v1->Neighbors.clear();
            UnorderedSwapInvalidFace(v0->Neighbors);
            if (v0->Neighbors.size() == 0) {
                //std::cout << "remove useless vertices" << std::endl;
                v0->Removed = true;
                return nDeleted;
            }
            // use v0 to store new vertex
            v0->Pos = edge.OptPos;
            v0->LOCAL_MARK = GLOBAL_MARK;
            EdgeCache.clear();

            v0->Q += v1->Q;
            Edge *e1 = nullptr;
            Edge *e2 = nullptr;

            for (auto f : v0->Neighbors) {
                f->GetEdgesByVertex(v0, e1, e2);
                syncIntoEdgeCache(e1, v0, EdgeCache);
                syncIntoEdgeCache(e2, v0, EdgeCache);
                if (e1->toBeReplaced != nullptr || e2->toBeReplaced != nullptr) {
                    f->replaceDuplicateEdge();
                }
            }
            for (auto e : EdgeCache) {
                if (e->toBeReplaced != nullptr) {
                    heap.deleteEdge(e);
                } else {
                    double priority = solveEdgePriority(*e);
                    heap.updateHeapEdgePriority(e, priority);
                    e->Start->toUpdateEdge = nullptr;
                }
            }
            v1->Neighbors.clear();
            return nDeleted;
        }


    private:

        static void syncIntoEdgeCache(Edge *e1, Vertex *v0, std::vector<Edge *> &EdgeCache) {
            if (e1->Start == v0) {
                e1->Start = e1->End;
                e1->End = v0;
            }


            if (e1->LOCAL_MARK != GLOBAL_MARK) {
                if (e1->Start->toUpdateEdge != nullptr && e1->Start->toUpdateEdge != e1) {
                    e1->toBeReplaced = e1->Start->toUpdateEdge;
                } else {
                    e1->Start->toUpdateEdge = e1;
                }
                e1->LOCAL_MARK = GLOBAL_MARK;
                EdgeCache.push_back(e1);
            }
        }
    };

    unsigned int CollapseHelper::GLOBAL_MARK = 0;
    double CollapseHelper::ScaleFactor = 1.0;
    SymetricMatrix CollapseHelper::tempQ;

    std::vector<Face *> CollapseHelper::FacesPool;
    size_t CollapseHelper::facePoolIdx = 0;

    std::vector<Edge *> CollapseHelper::EdgesPool;
    size_t CollapseHelper::edgePoolIdx = 0;


    std::vector<Vertex *> CollapseHelper::VertexPool;
    size_t CollapseHelper::vertexPoolIdx = 0;

    class MeshReducerPrivate {
    public:
        MeshReducerPrivate();

        ~MeshReducerPrivate();

        void reset();

        void reduce(unsigned int nTarget, bool bForceStrict = false);

        void load(const double *vertices, const uint16_t *indices, unsigned int nVert,
                  unsigned int nInd);

        void store(std::vector<double> &vertices, std::vector<uint16_t> &indices);

        bool isManifoldMesh() const;

        bool isValid() const;

        //Edge * spawnEdgeFromPool(Vertex * v0, Vertex * v1);
    private:
        class DuplicateVertexCmp {
        public:
            inline bool operator()(Vertex *const &rhs, Vertex *const &lhs) {
                {
                    return (rhs->Pos < lhs->Pos);
                }
                //            return (rhs->Pos == lhs->Pos ? (rhs < lhs) : (rhs->Pos < lhs->Pos));
            }
        };


        Edge *buildEdge(Vertex *v0, Vertex *v1);

        void buildQuadricMatrix();

        void initCollapses();

        void doFastLoop(unsigned int nTarget);

        void doStrictLoop(unsigned int nTarget);

        void cleanUp();

        std::vector<Vertex *> pVecResult;
        int validFaceNum = 0;
    private:
        // non-manifold ratio
        bool m_bStrictConstraint = false;

        // data section
        std::vector<Vertex *> Vertices;
        std::vector<Face *> Faces;

        std::vector<Vec3> EdgeOptList;
        std::vector<Vec3> EdgeDefaultList;


        // build edge for border setup
        std::unordered_map<uint32_t, Edge *> EdgeMap;


        std::vector<Edge *> Edges;
        std::vector<Edge *> EdgeCache;

        // Bounding box
        Vec3 OriginBBoxMin;
        Vec3 OriginBBoxMax;
    };

    MeshReducerPrivate::MeshReducerPrivate() {}

    MeshReducerPrivate::~MeshReducerPrivate() {}

    void MeshReducerPrivate::reset() {
        m_bStrictConstraint = false;
        Vertices.clear();
        //Faces.clear();
        Faces.clear();

        EdgeMap.clear();
        Edges.clear();
        pVecResult.clear();
        CollapseHelper::reset();
    }

    void MeshReducerPrivate::reduce(unsigned int nTarget, bool bForceStrict) {
        m_bStrictConstraint = m_bStrictConstraint || bForceStrict;

        // build Quadric matrix for each vertex
        buildQuadricMatrix();

        // compute the new position and quadric error
        initCollapses();

        // loop
        if (bForceStrict || m_bStrictConstraint) {
            doStrictLoop(nTarget);
        } else {
            doFastLoop(nTarget);
        }

        // clean up
        cleanUp();
    }

    // TODO: CNECO-2636 Find out whether the index type should change to uint32.
    void
    MeshReducerPrivate::load(const double *vertices, const uint16_t *indices, unsigned int nVert,
                             unsigned int nInd) {
        if (!vertices || !indices || !nVert || !nInd) {
            return;
        }

        // build vertices
        Vertices.resize(nVert);
        Vec3 min(std::numeric_limits<double>::max());
        Vec3 max(std::numeric_limits<double>::min());

        int i_vert_idx = 0;
        for (unsigned int i_vert = 0; i_vert < nVert; ++i_vert) {
            double x = vertices[i_vert_idx++];
            double y = vertices[i_vert_idx++];
            double z = vertices[i_vert_idx++];

            auto vertPtrs = CollapseHelper::spawnVertexFromPool();
            Vertices[i_vert] = vertPtrs;
            vertPtrs->Pos = Vec3(x, y, z);
            vertPtrs->ID = i_vert;

            // get scale
            min.X = std::min(min.X, x);
            max.X = std::max(max.X, x);

            min.Y = std::min(min.Y, y);
            max.Y = std::max(max.Y, y);

            min.Z = std::min(min.Z, z);
            max.Z = std::max(max.Z, z);
        }

        //CollapseHelper::ScaleFactor = double(1e8 * std::pow(1.0 / double((max - min).length()), 6));
        CollapseHelper::ScaleFactor = double(
                1e8 * std::pow(1.0 / double((max - min).squaredLength()), 3));


        OriginBBoxMax = max;
        OriginBBoxMin = min;

        // build faces
        unsigned int nFace = nInd / 3;

        Faces.resize(nFace);
        CollapseHelper::reserveFacePool(nFace);
        for (unsigned int idx = 0; idx < nFace; ++idx) {
            Faces[idx] = CollapseHelper::spawnFaceFromPool();
        }

        int i_face_vertIdx = 0;


        for (unsigned int i_face = 0; i_face < nFace; ++i_face) {
            auto &curr = *Faces[i_face];
            for (unsigned int j = 0; j < 3; ++j) {
                curr.Vertices[j] = Vertices[indices[i_face_vertIdx++]];
                curr.Vertices[j]->Neighbors.push_back(&curr);
            }

            for (unsigned int j = 0; j < 3; ++j) {
                Vertex *v0 = curr.Vertices[j];
                Vertex *v1 = curr.Vertices[(j + 1) % 3];

                Edge *e = buildEdge(v0, v1);
                if (e != nullptr)
                    e->AdjFaces++;
                curr.ConsistentEdges[j] = e;
            }
        }

        // manifold or non-manifold
        unsigned int nonManiEdge = 0;
        for (auto edge : Edges) {
            if (edge->AdjFaces == 1) {
                m_bStrictConstraint = true;
                break;
            }
        }
    }

    void MeshReducerPrivate::store(std::vector<double> &vertices, std::vector<uint16_t> &indices) {
        int valid_Vertex_Num = 0;
        for (auto &vertPtr : pVecResult) //scan pVec only to get rid of removed vertex
        {
            if (vertPtr->Removed) continue;
            valid_Vertex_Num++;
        }

        vertices.resize(valid_Vertex_Num * 3);
        unsigned int i_valid = 0;

        int vIdx = 0;
        for (auto &vert : pVecResult) {
            if (vert->Removed) {
                continue;
            }

            vert->ID = i_valid++;
            vertices[vIdx++] = vert->Pos.X;
            vertices[vIdx++] = vert->Pos.Y;
            vertices[vIdx++] = vert->Pos.Z;
        }


        indices.resize(this->validFaceNum * 3);
        int faceIndex = 0;
        for (auto face : Faces) {
            if (face->Valid == false) continue;
            indices[faceIndex++] = pVecResult[face->Vertices[0]->LOCAL_MARK]->ID;
            indices[faceIndex++] = pVecResult[face->Vertices[1]->LOCAL_MARK]->ID;
            indices[faceIndex++] = pVecResult[face->Vertices[2]->LOCAL_MARK]->ID;
        }

    }

    bool MeshReducerPrivate::isManifoldMesh() const {
        return m_bStrictConstraint;
    }


    Edge *MeshReducerPrivate::buildEdge(Vertex *v0, Vertex *v1) {
        // find first
        uint32_t key =
                v0->ID < v1->ID ? (uint32_t(v0->ID) << 16) | v1->ID : (uint32_t(v1->ID) << 16) |
                                                                      v0->ID;
        if (EdgeMap.find(key) != EdgeMap.end()) {
            auto e = EdgeMap[key];
            if (e != nullptr) {
                return e;
            } else {
                std::cout << "h" << std::endl;
            }
        }

        auto e = CollapseHelper::spawnEdgeFromPool(v0, v1);
        Edges.push_back(e);// .emplace_back(v0, v1);
        EdgeMap[key] = e;
        return e;
    }

    void MeshReducerPrivate::buildQuadricMatrix() {

        for (auto &face : Faces) {
            // ax + by + cz + d = 0
            Vec3 normal = face->normal();
            if (!m_bStrictConstraint) {
                normal = normal.normalize();
            }

            double d = -normal.dot(face->Vertices[0]->Pos);

            // assemble quadric matrix
            SymetricMatrix Q(normal.X, normal.Y, normal.Z, d);
            for (int i = 0; i < 3; ++i) {
                face->Vertices[i]->Q += Q;
            }
        }

        if (m_bStrictConstraint) {
            //int edgeSize = Edges.size();
            int faceSize = Faces.size();
            for (int fidx = 0; fidx < faceSize; fidx++) {
                auto &face = *Faces[fidx];
                if (face.ConsistentEdges[0]->AdjFaces != 1 &&
                    face.ConsistentEdges[1]->AdjFaces != 1 &&
                    face.ConsistentEdges[2]->AdjFaces != 1)
                    continue;;
                Vec3 normal = face.normal();
                for (unsigned int j = 0; j < 3; ++j) {
                    if (face.ConsistentEdges[j]->AdjFaces == 1) {
                        Vertex *start = face.Vertices[j];
                        Vertex *end = face.Vertices[(j + 1) % 3];

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

    void MeshReducerPrivate::initCollapses() {
        if (m_bStrictConstraint) {
            for (auto &edge : Edges) {
                CollapseHelper::solve(*edge);
            }
        } else {
            CollapseHelper::GLOBAL_MARK++;
            for (auto &face : Faces) {
                CollapseHelper::solveFaceOptimized(*face);
            }
        }
    }

    void MeshReducerPrivate::doFastLoop(unsigned int nTarget) {
        unsigned int faceCount = static_cast < unsigned int> (Faces.size());
        unsigned int nowCount = faceCount;

        double minPriority = std::numeric_limits<double>::max();
        {
            for (auto &face : Faces) {
                double d = face->priority;
                if (d < minPriority) {
                    minPriority = d;
                }
            }
        }
        int removedCount = 0;
        for (unsigned int iter = 0; iter < MAX_ITERATION; ++iter) {
            if (nowCount <= nTarget) {
                //std::cout << "Operation end at iter " << iter << std::endl;
                break;
            }

            //std::cout << "Operation at iter " << iter << std::endl;
            double threshold = 1e-9 * static_cast<double>(std::pow(iter + 3, AGGRESSIVE));

            {
                if (threshold <= minPriority) {
                    continue;
                }
            }
            //std::cout << "threshold " << threshold << std::endl;
            int start_mark = CollapseHelper::GLOBAL_MARK;
            int deletePerRound = 0;
            int bubleCount = 0;
            int startMark = CollapseHelper::GLOBAL_MARK;
            for (auto &face : Faces) {
                // get the top heap element
                if (!face->Valid || face->LOCAL_MARK > startMark || face->priority > threshold) {
                    continue;
                }

                // update
                unsigned int nDeleted = CollapseHelper::update(*face);
                deletePerRound += nDeleted;

                nowCount -= nDeleted;
                if (nowCount <= nTarget) {
                    break;
                }
            }


            if (nowCount > nTarget && deletePerRound > 0) {
                bubleCount += deletePerRound;
                //here I would just bubble out invalid faces
                if (bubleCount > Faces.size() >> 2 && bubleCount > 100) { //1/4 bubble which is bad
                    //CollapseHelper::UnorderedSwapInvalidFace(Faces);
                    bubleCount = 0;
                }
                minPriority = std::numeric_limits<double>::max();
                // clear dirty flag
                for (auto &face : Faces) {
                    //face->Dirty = false;
                    if (face->Valid) {
                        double d = face->priority;
                        if (d < minPriority) {
                            minPriority = d;
                        }
                    }
                }
            }
        }
    }

    void MeshReducerPrivate::doStrictLoop(unsigned int nTarget) {
        unsigned int faceCount = static_cast < unsigned int> (Faces.size());
        unsigned int nowCount = faceCount;
        MinEdgeHeap heap(Edges);

        // clock
        auto start = std::chrono::steady_clock::now();
        auto end = start;// std::chrono::steady_clock::now();

        int loopRound = 0;
        ////int first10 = 10000;
        // collapse in loop
        while ((nowCount > nTarget) && heap.size != 0 &&
               std::chrono::duration_cast<std::chrono::seconds>(end - start).count() < TIMEOUT) {
            Edge &top = *heap.extractMin();


            if (top.Start->Removed || top.End->Removed ||
                top.LOCAL_MARK < top.Start->LOCAL_MARK || top.LOCAL_MARK < top.End->LOCAL_MARK) {
                //here means that isolated edges exist, all faces are not valid
                //this is the interest part. as we could find many corner cases
                continue;
            }

            // update time
            loopRound++;
            if (loopRound & 1) {
                end = std::chrono::steady_clock::now();
            }
            // update
            unsigned int nDeleted = CollapseHelper::updateEdge(heap, top, EdgeCache);
            nowCount -= nDeleted;
        }
    }

    void MeshReducerPrivate::cleanUp() {

        std::vector<Vertex *> &pVec = pVecResult;
        pVec.clear();
        validFaceNum = 0;

        int valid_Vertex_Num = 0;
        for (auto &vertPtr : Vertices) {
            if (vertPtr->Removed) continue;
            valid_Vertex_Num++;
        }
        if (valid_Vertex_Num == 0) {
            return;
        }
        pVec.resize(valid_Vertex_Num);
        unsigned int i_valid = 0;
        for (auto &vertPtr : Vertices) {
            auto &vert = *vertPtr;
            if (vert.Removed) continue;
            pVec[i_valid++] = &vert;
        }
        DuplicateVertexCmp cmp;
        std::sort(pVec.begin(), pVec.end(), cmp);

        unsigned int j = 0;
        unsigned int i = 1;
        unsigned int nDeleted = 0;
        for (; i != i_valid;) {
            pVec[i]->Removed = true;
            if (pVec[i]->Pos == pVec[j]->Pos) {
                pVec[i]->LOCAL_MARK = j; //the index in pVec
                i++;
                nDeleted++;
            } else {
                pVec[i]->LOCAL_MARK = i;
                j = i;
                ++i;
            }
        }
        pVec[0]->LOCAL_MARK = 0;
        pVec[0]->Removed = true;
        //std::cout << "V2 remove duplicate " << nDeleted << std::endl;
        //****************************************************************8

        for (auto &face : Faces) //all face valid here
        {
            if (face->Valid == false) continue;

            auto v0 = face->Vertices[0];// ->Pos;
            auto v1 = face->Vertices[1];//->Pos;
            auto v2 = face->Vertices[2];//->Pos;

            if (v0->LOCAL_MARK == v1->LOCAL_MARK ||
                v0->LOCAL_MARK == v2->LOCAL_MARK ||
                v1->LOCAL_MARK == v2->LOCAL_MARK) {
                face->Valid = false; //if share same vertex id, means face invalid.. Equivalent to original removeFaceOutOfRangeArea
                continue;
            }

            {
                validFaceNum++;

                //following equivalent to removeUnreferenceVertex by set the correct Remove..
                pVec[v0->LOCAL_MARK]->Removed = false;
                pVec[v1->LOCAL_MARK]->Removed = false;
                pVec[v2->LOCAL_MARK]->Removed = false;
            }
        }
    }

    bool MeshReducerPrivate::isValid() const {
        Vec3 min(std::numeric_limits<double>::max());
        Vec3 max(std::numeric_limits<double>::min());
        for (auto vertPtr : Vertices) {
            auto &vert = *vertPtr;
            if (vert.Removed) {
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

        double len_diag = (max - min).squaredLength();
        double len_diag_old = (OriginBBoxMax - OriginBBoxMin).squaredLength();

        static const double VALID_THRESHOLD_SQUARE = VALID_THRESHOLD * VALID_THRESHOLD;

        double minDiag = std::min(len_diag, len_diag_old);
        double maxDiag = std::max(len_diag, len_diag_old);
        return minDiag > VALID_THRESHOLD_SQUARE * maxDiag;
    }
} // namespace

namespace common {
    static MeshReducerPrivate reducer;
    bool MeshReducer::reduce(const double *vertices, const uint16_t *indices, unsigned int nVert,
                             unsigned int nIdx,
                             std::vector<double> &reducedVertices,
                             std::vector<uint16_t> &reducedIndices,
                             unsigned int nTarget) {
        if (vertices == nullptr || indices == nullptr || nVert == 0 || nIdx == 0 || nTarget == 0) {
            return false;
        }
        reducer.reset();

        reducer.load(vertices, indices, nVert, nIdx); //load take 1/3 of total time
        reducer.reduce(nTarget);        //take 2/3 of total time
        bool bValid = reducer.isValid();
        if (!bValid && !reducer.isManifoldMesh()) {
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
