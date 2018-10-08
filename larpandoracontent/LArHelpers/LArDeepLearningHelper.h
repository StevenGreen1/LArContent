/**
 *  @file   larpandoracontent/LArHelpers/LArDeepLearningHelper.h
 *
 *  @brief  Header file for the deep learning helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DEEP_LEARNING_HELPER_H
#define LAR_DEEP_LEARNING_HELPER_H 1

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArDeepLearning/LArKerasModel.h"

namespace lar_content
{

/**
 *  @brief  LArDeepLearningHelper class
 */
class LArDeepLearningHelper
{
public:
    /**
     *  @brief  Function to convert Pandora histogram to Keras data block
     *
     *  @param  twoDHistogram histogram to convert
     *  @param  dataBlock2D data block to populate
     */
    static void HistogramToDataBlock(const pandora::TwoDHistogram &twoDHistogram, KerasModel::DataBlock2D &dataBlock2D);
};

} // namespace lar_content

#endif // #ifndef LAR_DEEP_LEARNING_HELPER_H
