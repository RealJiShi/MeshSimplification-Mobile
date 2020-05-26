/*
** Copyright (c) 2019 Qualcomm Technologies, Inc.
** All Rights Reserved.
** Confidential and Proprietary - Qualcomm Technologies, Inc.
**
** Code from https://github.com/rawrunprotected/rasterizer/blob/master/SoftwareRasterizer/Occluder.h under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
*/

#pragma once

#include <memory>
#include <map>
#include <vector>
#include "SSE2NEON.h"

namespace common
{

struct Occluder
{
    static std::shared_ptr<Occluder> bake(const std::vector<__m128> &vertices, __m128 refMin, __m128 refMax);

    // center
    __m128 Center;

    // Bounding Box min/max
    __m128 BoundsMin;
    __m128 BoundsMax;

    std::vector<__m128i> Vertices;
};

typedef std::shared_ptr<Occluder> OccluderPtr;

} // namespace common
