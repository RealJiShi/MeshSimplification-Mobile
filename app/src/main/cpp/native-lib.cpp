#include <jni.h>
#include <chrono>
#include <map>
#include <string>
#include <android/log.h>
#include "MeshSimplifcation/MeshReducer.h"
#include "MeshSimplifcation/OffFileHelper.h"
#include "MeshSimplifcation/QuadDecomposition.h"
#include "MeshSimplifcation/Occluder.h"

#define  LOG_TAG    "MeshSimplication"
#define  LOG(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)

struct MeshInfo
{
    std::vector<double> Vertices;
    std::vector<uint16_t> Indices;
    float Time = 0.0f;
};

class AABB
{
public:
    AABB(const std::vector<__m128> &Vertices)
    {
        if (Vertices.empty())
        {
            Min = _mm_setzero_ps();
            Max = _mm_setzero_ps();
        }
        else
        {
            Min = _mm_set1_ps(+std::numeric_limits<float>::infinity());
            Max = _mm_set1_ps(-std::numeric_limits<float>::infinity());

            for (auto v : Vertices)
            {
                Min = _mm_min_ps(Min, v);
                Max = _mm_max_ps(Max, v);
            }
        }
    }

    __m128 getCenter() const
    {
        return _mm_add_ps(Min, Max);
    }

    __m128 getExtents() const
    {
        return _mm_sub_ps(Min, Max);
    }

    __m128 surfaceArea()
    {
        __m128 extents = getExtents();
        __m128 extents2 = _mm_shuffle_ps(extents, extents, _MM_SHUFFLE(3, 0, 2, 1));
        return _mm_dp_ps(extents, extents2, 0x7F);
    }

    // max, min
    __m128 Min, Max;
};

bool load(const std::string &file, std::vector<double> &Vertices, std::vector<uint16_t> &Indices)
{
    static const std::string PREFIX = "/sdcard/Original/";
    static const std::string SUFFIX = ".off";
    if (!OffFileHelper::load(PREFIX + file + SUFFIX, Vertices, Indices))
    {
        return false;
    }
    return true;
}

bool process(const std::string &file, std::vector<double> &Vertices, std::vector<uint16_t> &Indices,
             float &runtime, bool bIsSave = false, float ratio = 0.25f)
{
    static std::string OUPTUT = "/sdcard/LOD/";
    if (Vertices.empty() || Indices.empty())
    {
        return false;
    }

    // simplification
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> reducedVertices;
    std::vector<uint16_t > reducedIndices;

    unsigned int nVert = Vertices.size() / 3;
    unsigned int nIdx = Indices.size();
    unsigned int nTarget = static_cast<unsigned int>(float(nIdx) / 3.0f * 0.25f);
    bool bValid = common::MeshReducer::reduce(&Vertices[0], &Indices[0], nVert, nIdx,
            reducedVertices, reducedIndices, nTarget);
//
//    // keep consistent with QOC code
//    std::vector<__m128> VerticesVec;
//    std::vector<unsigned int> IndicesVec;
//    if (!bValid)
//    {
//        VerticesVec.resize(Vertices.size() / 3);
//        IndicesVec.resize(Indices.size());
//        for (unsigned int i_vert = 0; i_vert < Vertices.size() / 3; ++i_vert)
//        {
//            VerticesVec[i_vert] = _mm_setr_ps(Vertices[3 * i_vert + 0],
//                                              Vertices[3 * i_vert + 1],
//                                              Vertices[3 * i_vert + 2],
//                                              1.0f);
//        }
//
//        for (unsigned int i_index = 0; i_index < Indices.size(); ++i_index)
//        {
//            IndicesVec[i_index] = Indices[i_index];
//        }
//    }
//    else
//    {
//        VerticesVec.resize(reducedVertices.size() / 3);
//        IndicesVec.resize(reducedIndices.size());
//        for (unsigned int i_vert = 0; i_vert < reducedVertices.size() / 3; ++i_vert)
//        {
//            VerticesVec[i_vert] = _mm_setr_ps(reducedVertices[3 * i_vert + 0],
//                                              reducedVertices[3 * i_vert + 1],
//                                              reducedVertices[3 * i_vert + 2],
//                                              1.0f);
//        }
//
//        for (unsigned int i_index = 0; i_index < reducedIndices.size(); ++i_index)
//        {
//            IndicesVec[i_index] = reducedIndices[i_index];
//        }
//    }
//
//    // Quad Decomposition
//    IndicesVec = common::QuadDecomposition::decompose(IndicesVec, VerticesVec);
//
//    // padding to a multiple of 4 quads
//    // while ((IndicesVec.size() % 16) != 0)
//    while ((IndicesVec.size() & 15) != 0)
//    {
//        IndicesVec.push_back(IndicesVec[0]);
//    }
//
//    // Get AABB info
//    AABB aabb(VerticesVec);
//
//    // expand the data
//    std::vector<__m128> VerticesVecExt;
//    VerticesVecExt.reserve(IndicesVec.size());
//    for (auto index : IndicesVec)
//    {
//        VerticesVecExt.push_back(VerticesVec[index]);
//    }
//    auto occluderBaked = common::Occluder::bake(VerticesVecExt, aabb.Min, aabb.Max);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    runtime = duration.count() / 1000.0f;

    // only save at the last loop
    if (bIsSave)
    {
        std::string output = file + "LOD.off";
        if (!OffFileHelper::store(OUPTUT + output, reducedVertices, reducedIndices))
        {
            return false;
        }
    }
    return true;
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_qualcomm_meshsimplication_MainActivity_stringFromJNI(
        JNIEnv* env,
        jobject /* this */) {
    static const int NUM_MESH = 291;
    static const int NUM_LOOP = 500;
    std::map<std::string, MeshInfo> MeshMap;

    // load first
    for (int index = 0; index < NUM_MESH; ++index)
    {
        std::string filename = std::string("mesh") + std::to_string(index);
        MeshMap.insert(std::make_pair(filename, MeshInfo()));
        auto &curr = MeshMap[filename];
        load(filename, curr.Vertices, curr.Indices);
    }

    // process
    for (int loop = 0; loop < NUM_LOOP; ++loop)
    {
        LOG("Current Loop: %d", loop);
        for (auto &mesh : MeshMap)
        {
            float runtime = 0.0f;
            process(mesh.first, mesh.second.Vertices, mesh.second.Indices, runtime, (loop == (NUM_LOOP - 1)));
            mesh.second.Time += runtime;
        }
    }

    for (auto &mesh : MeshMap)
    {
        mesh.second.Time /= NUM_LOOP;
    }

    // output
    std::ofstream stat;
    stat.open("/sdcard/LOD/LodStat.csv");
    if (stat.is_open())
    {
        // header
        stat << "Mesh Name, Process Avg Time\n";
        for (auto &rec : MeshMap)
        {
            stat << rec.first << "," << rec.second.Time << "\n";
        }
        stat.close();
    }

    std::string hello = "Finished";
    return env->NewStringUTF(hello.c_str());
}
