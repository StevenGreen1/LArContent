/**
 *  @file   larpandoracontent/LArTrackShowerId/DLPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the deep learning based pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DLPfoCharacterisationAlgorithm class
 */
class DLPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
private:
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
