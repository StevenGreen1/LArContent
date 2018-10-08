/**
 *  @file   larpandoracontent/LArHelpers/LArDeepLearningHelper.cc
 *
 *  @brief  Implementation of the deep learning helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArDeepLearningHelper.h"

#include <limits>

using namespace pandora;

namespace lar_content
{

void LArDeepLearningHelper::HistogramToDataBlock(const TwoDHistogram &twoDHistogram, KerasModel::DataBlock2D &dataBlock2D)
{
    Data3D data3D;
    Data2D data2D;
    for (int yBin = 0; yBin < twoDHistogram.GetNBinsY(); yBin++)
    {
        Data1D data1D;
        for (int xBin = 0; xBin < twoDHistogram.GetNBinsX(); xBin++)
        {
            data1D.Append(twoDHistogram.GetBinContent(xBin, yBin) * 256.f / 10000.f );
        }
        data2D.Append(data1D);
    }
    data3D.Append(data2D);
    dataBlock2D.SetData(data3D);
    return;
}

} // namespace lar_content
