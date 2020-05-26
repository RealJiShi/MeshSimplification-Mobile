/*
** Copyright (c) 2019 Qualcomm Technologies, Inc.
** All Rights Reserved.
** Confidential and Proprietary - Qualcomm Technologies, Inc.
**
** Code from https://github.com/rawrunprotected/rasterizer/blob/master/SoftwareRasterizer/Occluder.cpp under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
*/

#include "Occluder.h"

#include <cassert>
#include <limits>

#include "MathUtil.h"
#include "SSE2NEON.h"

namespace
{
using namespace common;
__m128 normalizeQuadNormal(const __m128 &v0, const __m128 &v1, const __m128 &v2, const __m128 &v3)
{
    const __m128 &edge0 = _mm_sub_ps(v1, v0);
    const __m128 &edge1 = _mm_sub_ps(v2, v0);
    const __m128 &edge2 = _mm_sub_ps(v3, v0);
    return normalize(_mm_add_ps(cross(edge0, edge1), cross(edge1, edge2)));
}
} // namespace

namespace common
{

std::shared_ptr<Occluder> Occluder::bake(const std::vector<__m128> &vertices, __m128 refMin, __m128 refMax)
{
    assert((vertices.size() & 15) == 0);

    // Simple k-means clustering by normal direction to improve backface culling efficiency
    std::vector<__m128> quadNormals;
    const int vert_num = vertices.size();
    quadNormals.resize(vert_num / 4);
    for (unsigned int i = 0; i < vert_num; i += 4)
    {
        auto v0 = vertices[i + 0];
        auto v1 = vertices[i + 1];
        auto v2 = vertices[i + 2];
        auto v3 = vertices[i + 3];

        quadNormals[i / 4] = normalizeQuadNormal(v0, v1, v2, v3);
    }

    std::vector<__m128> centroids(6);
    std::vector<uint32_t> centroidAssignment;
    centroids[0] = _mm_setr_ps(+1.0f, 0.0f, 0.0f, 0.0f);
    centroids[1] = _mm_setr_ps(0.0f, +1.0f, 0.0f, 0.0f);
    centroids[2] = _mm_setr_ps(0.0f, 0.0f, +1.0f, 0.0f);
    centroids[3] = _mm_setr_ps(0.0f, -1.0f, 0.0f, 0.0f);
    centroids[4] = _mm_setr_ps(0.0f, 0.0f, -1.0f, 0.0f);
    centroids[5] = _mm_setr_ps(-1.0f, 0.0f, 0.0f, 0.0f);

    centroidAssignment.resize(vert_num / 4);

    bool anyChanged = true;
    const unsigned int norm_num = quadNormals.size();
    for (unsigned int iter = 0; iter < 10 && anyChanged; ++iter)
    {
        anyChanged = false;
        for (unsigned int j = 0; j < norm_num; ++j)
        {
            __m128 normal = quadNormals[j];

            __m128 bestDistance = _mm_set1_ps(-std::numeric_limits<float>::infinity());
            unsigned int bestCentroid = 10000;
            for (unsigned int k = 0; k < centroids.size(); ++k)
            {
                __m128 distance = _mm_dp_ps(centroids[k], normal, 0x7F);
                if (_mm_comige_ss(distance, bestDistance))
                {
                    bestDistance = distance;
                    bestCentroid = k;
                }
            }

            if (centroidAssignment[j] != bestCentroid)
            {
                centroidAssignment[j] = bestCentroid;
                anyChanged = true;
            }
        }
        std::fill(centroids.begin(), centroids.end(), _mm_setzero_ps());

        for (unsigned int j = 0; j < norm_num; ++j)
        {
            unsigned int k = centroidAssignment[j];

            if (k < 6)
            {
                centroids[k] = _mm_add_ps(centroids[k], quadNormals[j]);
            }
        }

        for (unsigned int k = 0; k < centroids.size(); ++k)
        {
            centroids[k] = normalize(centroids[k]);
        }
    }

    std::vector<__m128> orderedVertices;
    for (unsigned int j = 0; j < norm_num; ++j)
    {
        unsigned int val_j = centroidAssignment[j];
        if (val_j < centroids.size())
        {
            orderedVertices.push_back(vertices[4 * j + 0]);
            orderedVertices.push_back(vertices[4 * j + 1]);
            orderedVertices.push_back(vertices[4 * j + 2]);
            orderedVertices.push_back(vertices[4 * j + 3]);
        }
    }

    auto occluder = std::make_shared<Occluder>();

    __m128 invExtents = _mm_div_ps(_mm_set1_ps(1.0f), _mm_sub_ps(refMax, refMin));

    __m128 scalingX = _mm_set1_ps(2047.0f);
    __m128 scalingY = _mm_set1_ps(2047.0f);
    __m128 scalingZ = _mm_set1_ps(1023.0f);

    __m128 half = _mm_set1_ps(0.5f);

    const unsigned int NUM_ORDERED_VERT = orderedVertices.size();
    for (unsigned int i = 0; i < NUM_ORDERED_VERT; i += 16)
    {
        for (unsigned int j = 0; j < 4; ++j)
        {
            // Transform into [0,1] space relative to bounding box
            __m128 v0 = _mm_mul_ps(_mm_sub_ps(orderedVertices[i + j + 0], refMin), invExtents);
            __m128 v1 = _mm_mul_ps(_mm_sub_ps(orderedVertices[i + j + 4], refMin), invExtents);
            __m128 v2 = _mm_mul_ps(_mm_sub_ps(orderedVertices[i + j + 8], refMin), invExtents);
            __m128 v3 = _mm_mul_ps(_mm_sub_ps(orderedVertices[i + j + 12], refMin), invExtents);

            // Transpose into [xxxx][yyyy][zzzz][wwww]
            _MM_TRANSPOSE4_PS(v0, v1, v2, v3);

            // Scale and truncate to int
            v0 = _mm_fmadd_ps(v0, scalingX, half);
            v1 = _mm_fmadd_ps(v1, scalingY, half);
            v2 = _mm_fmadd_ps(v2, scalingZ, half);

            __m128i X = _mm_cvttps_epi32(v0);
            __m128i Y = _mm_cvttps_epi32(v1);
            __m128i Z = _mm_cvttps_epi32(v2);

            // Pack to 11/11/10 format
            __m128i XYZ = _mm_or_si128(_mm_slli_epi32(X, 21), _mm_or_si128(_mm_slli_epi32(Y, 10), Z));

            occluder->Vertices.push_back(XYZ);
        }
    }

    // set bbox as homogeneous coordinates
    occluder->BoundsMin = _mm_set_w(refMin, 1.0f);
    occluder->BoundsMax = _mm_set_w(refMax, 1.0f);
    occluder->Center = _mm_mul_ps(_mm_add_ps(refMax, refMin), _mm_set1_ps(0.5f));

    return occluder;
}

} // namespace common
