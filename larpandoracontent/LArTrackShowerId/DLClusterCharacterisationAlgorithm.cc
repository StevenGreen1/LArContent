/**
 *  @file   larpandoracontent/LArTrackShowerId/DLClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning based cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTrackShowerId/DLClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

bool DLClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    int nShowerHits(0), nTrackHits(0);

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;
            // ATTN: Parent calo hit is owned by master algorithm, which has correct metadata
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

StatusCode DLClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
