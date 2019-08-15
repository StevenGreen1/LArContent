/**
 *  @file   larpandoracontent/LArTrackShowerId/DLClusterCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the deep learning based cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_DL_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DLClusterCharacterisationAlgorithm class
 */
class DLClusterCharacterisationAlgorithm : public ClusterCharacterisationBaseAlgorithm
{
private:
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_DL_CLUSTER_CHARACTERISATION_ALGORITHM_H
