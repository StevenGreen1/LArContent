/**
 *  @file   larpandoracontent/LArTrackShowerId/DLPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTrackShowerId/DLPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

bool DLPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    int nShowerHits(0), nTrackHits(0);

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;
            const CaloHit *pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));
            const PropertiesMap &properties(pParentCaloHit->GetPropertiesMap());
            const bool isTrack(properties.find("IsTrack") == properties.end() ? false : true);

            if (isTrack)
            {
                nTrackHits++;
            }
            else
            {
                nShowerHits++;
            }
        }
    }

    if (nTrackHits >= nShowerHits)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const /*pPfo*/) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
