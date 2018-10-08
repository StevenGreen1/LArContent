/**
 *  @file   larpandoracontent/LArTwoDReco/ClusterSplitting/DeepLearningSplittingAlgorithm.h
 *
 *  @brief  Header file for the layer splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DEEP_LEARNING_SPLITTING_ALGORITHM_H
#define LAR_DEEP_LEARNING_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeepLearningSplittingAlgorithm class
 */
class DeepLearningSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DeepLearningSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
        pandora::CaloHitList &secondCaloHitList) const;

    KerasModel    m_kerasModel;               ///< Keras model
    std::string   m_cnnModelXml;              ///< Model xml
    std::string   m_cnnModelName;             ///< Model name
    int           m_gridSize;                 ///< Number of bins in grid
    float         m_gridDimensions;           ///< Physical dimension of grid
};

} // namespace lar_content

#endif // #ifndef LAR_DEEP_LEARNING_SPLITTING_ALGORITHM_H
