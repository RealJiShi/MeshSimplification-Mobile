#include <jni.h>
#include <chrono>
#include <map>
#include <string>
#include <android/log.h>
#include "MeshSimplifcation/MeshReducer.h"
#include "MeshSimplifcation/OffFileHelper.h"

#define  LOG_TAG    "MeshSimplication"
#define  LOG(...)  __android_log_print(ANDROID_LOG_INFO,LOG_TAG,__VA_ARGS__)

struct MeshInfo
{
    std::vector<double> Vertices;
    std::vector<uint16_t> Indices;
    float Time = 0.0f;
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
    common::MeshReducer::reduce(&Vertices[0], &Indices[0], Vertices.size() / 3, Indices.size(),
                                reducedVertices, reducedIndices, nTarget);
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
