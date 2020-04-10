#include <jni.h>
#include <chrono>
#include <map>
#include <string>
#include "MeshSimplifcation/MeshReducer.h"
#include "MeshSimplifcation/OffFileHelper.h"

bool process(const std::string &file, float &runtime, float ratio = 0.25f)
{
    static std::string PREFIX = "/sdcard/Original/";
    static std::string OUPTUT = "/sdcard/LOD/";
    static std::string SUFFIX = ".off";
    std::vector<double> Vertices;
    std::vector<uint16_t> Indices;
    if (!OffFileHelper::load(PREFIX + file + SUFFIX, Vertices, Indices))
    {
        return false;
    }

    // simplification
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> reducedVertices;
    std::vector<uint16_t > reducedIndices;
    common::MeshReducer::reduce(&Vertices[0], &Indices[0], Vertices.size(), Indices.size(),
            reducedVertices, reducedIndices, Indices.size() / 3 * 0.25f);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    runtime = duration.count() / 1000.0f;

    //
    std::string output = file + "LOD.off";
    if (!OffFileHelper::store(OUPTUT + output, reducedVertices, reducedIndices))
    {
        return false;
    }
    return true;
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_qualcomm_meshsimplication_MainActivity_stringFromJNI(
        JNIEnv* env,
        jobject /* this */) {

    static int NUM_MESH = 299;
    std::map<std::string, float> RecorderMap;
    for (int index = 0; index < NUM_MESH; ++index)
    {
        std::string filename = std::string("mesh") + std::to_string(index);
        RecorderMap.insert(std::make_pair(filename, 0.0f));
    }

    static const int NUM_LOOP = 100;
    for (int loop = 0; loop < NUM_LOOP; ++loop)
    {
        for (auto &mesh : RecorderMap)
        {
            float runtime = 0.0f;
            process(mesh.first, runtime);
            mesh.second += runtime;
        }
    }

    for (auto &mesh : RecorderMap)
    {
        mesh.second /= NUM_LOOP;
    }

    // output
    std::ofstream stat;
    stat.open("/sdcard/LOD/LodStat.csv");
    if (stat.is_open())
    {
        // header
        stat << "Mesh Name, Process Avg Time\n";
        for (auto &rec : RecorderMap)
        {
            stat << rec.first << "," << rec.second << "\n";
        }
        stat.close();
    }

    std::string hello = "Finished";
    return env->NewStringUTF(hello.c_str());
}
