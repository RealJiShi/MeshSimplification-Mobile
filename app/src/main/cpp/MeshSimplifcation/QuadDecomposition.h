/*
** Copyright (c) 2019 Qualcomm Technologies, Inc.
** All Rights Reserved.
** Confidential and Proprietary - Qualcomm Technologies, Inc.
**
** Code from https://github.com/rawrunprotected/rasterizer/blob/master/SoftwareRasterizer/QuadDecomposition.h under CC0 1.0 Universal (CC0 1.0) Public Domain Dedication.
*/

#pragma once

#include <vector>
#include "SSE2NEON.h"

namespace common
{
class QuadDecomposition
{
public:
    static std::vector<uint32_t> decompose(const std::vector<uint32_t> &indices, const std::vector<__m128> &vertices);
    static void release();
};
} // namespace common
