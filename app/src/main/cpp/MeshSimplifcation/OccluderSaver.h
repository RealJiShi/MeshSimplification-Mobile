#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include "Occluder.h"
#include "SSE2NEON.h"

class OccluderSaver
{
public:
    static bool store(const std::string &file_path, const common::OccluderPtr &occluder)
    {
        static __m128i maskY = _mm_set1_epi32(2047);
        static __m128i maskZ = _mm_set1_epi32(1023);
        std::ofstream fout(file_path);
        if (!fout)
        {
            return false;
        }

        float x[4], y[4], z[4];
        // int32_t data[4];
        for (auto vert: occluder->Vertices)
        {
             // _mm_storeu_si128((__m128i *)data, vert);
             // fout << data[0] << " " << data[1] << " " << data[2] << " " << data[3] << std::endl;

            // unpack & store
            __m128 X = _mm_cvtepi32_ps(_mm_srli_epi32(vert, 21));
            __m128 Y = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(vert, 10), maskY));
            __m128 Z = _mm_cvtepi32_ps(_mm_and_si128(vert, maskZ));

            _mm_storeu_ps(x, X);
            _mm_storeu_ps(y, Y);
            _mm_storeu_ps(z, Z);

            // store
            fout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
            fout << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << std::endl;
            fout << z[0] << " " << z[1] << " " << z[2] << " " << z[3] << std::endl;
        }
        fout.close();
        return true;
    }

    static bool load(const std::string &file_path, std::vector<int> &data)
    {
        std::ifstream fin(file_path);
        if (!fin)
        {
            return false;
        }

        std::string line;
        int tmp[4];
        while (std::getline(fin, line))
        {
            std::stringstream info;
            info.str(line);
            for (int i = 0; i < 4; ++i)
            {
                int tmp;
                info >> tmp;
                data.push_back(tmp);
            }
        }
        fin.close();
        return true;
    }
};