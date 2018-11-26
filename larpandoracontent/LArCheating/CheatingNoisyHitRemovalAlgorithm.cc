/**
 *  @file   larpandoracontent/LArCheating/CheatingNoisyHitRemovalAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating noisy hit removal algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArCheating/CheatingNoisyHitRemovalAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingNoisyHitRemovalAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    CaloHitList outputCaloHitList;

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        const HitType hitType(pCaloHit->GetHitType());
        const float z(pCaloHit->GetPositionVector().GetZ());
        bool isInGap(false);

        for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
        {
            const LineGap *const pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

            if (!pLineGap)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            const LineGapType lineGapType(pLineGap->GetLineGapType());

            if (!(((TPC_VIEW_U == hitType) && (TPC_WIRE_GAP_VIEW_U == lineGapType)) ||
                 ((TPC_VIEW_V == hitType) && (TPC_WIRE_GAP_VIEW_V == lineGapType)) ||
                  ((TPC_VIEW_W == hitType) && (TPC_WIRE_GAP_VIEW_W == lineGapType))))
            {
                continue;
            }

            if ((pLineGap->GetLineStartZ() < z) && (pLineGap->GetLineEndZ() > z))
                isInGap = true;
        }

        if (!isInGap)
            outputCaloHitList.push_back(pCaloHit);
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, outputCaloHitList, m_outputCaloHitListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNoisyHitRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
