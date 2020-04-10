/*========================================================================*/
/**
\file MeshReducer.h

    Copyright (c) 2019 Qualcomm Technologies, Inc.
    All rights reserved.
    Confidential and Proprietary - Qualcomm Technologies, Inc.
*/
/*========================================================================*/

#pragma once

#include <vector>
#include <stdint.h>
#include "SSE2NEON.h"

namespace common
{
class MeshReducer
{
public:
    static bool reduce(const double *vertices, const uint16_t *indices, unsigned int nVert, unsigned int nIdx,
                       std::vector<double> &reducedVertices, std::vector<uint16_t> &reducedIndices,
                       unsigned int nTarget);

    MeshReducer(const MeshReducer &) = delete;
    MeshReducer(MeshReducer &&) = delete;
    MeshReducer &operator=(const MeshReducer &) = delete;
    MeshReducer &operator=(MeshReducer &&) = delete;

private:
    MeshReducer();
    ~MeshReducer();
};
} // namespace common
