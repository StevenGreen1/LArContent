/**
 *  @file   larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CaloHitEventDumpAlgorithm::CaloHitEventDumpAlgorithm() :
    m_gridSize(16),
    m_gridDimensions(50),
    m_useTrainingMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitEventDumpAlgorithm::Run()
{
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    CaloHitList wCaloHitList;

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        if (TPC_VIEW_W == pCaloHit->GetHitType())
            wCaloHitList.push_back(pCaloHit);
    }

    if (m_useTrainingMode)
    {
        for (const CaloHit *pTargetCaloHit : wCaloHitList)
        {
            LArMvaHelper::MvaFeatureVector featureVector;
            int targetParticleId(0);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pTargetCaloHit));
                targetParticleId = pMCParticle->GetParticleId();
            }
            catch(...)
            {
                continue;
            }

            featureVector.push_back(static_cast<double>(targetParticleId));

            TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);

            for (const CaloHit *pNeighbourCaloHit : wCaloHitList)
            {
                CartesianVector relativePosition(pNeighbourCaloHit->GetPositionVector() - pTargetCaloHit->GetPositionVector());
                twoDHistogram.Fill(relativePosition.GetX(), relativePosition.GetZ(), pNeighbourCaloHit->GetInputEnergy());
            }

            for (int xBin = 0; xBin < twoDHistogram.GetNBinsX(); xBin++)
            {
                for (int yBin = 0; yBin < twoDHistogram.GetNBinsY(); yBin++)
                {
                    featureVector.push_back(static_cast<double>(twoDHistogram.GetBinContent(xBin, yBin)));
                }
            }

            // ATTN: The hard coded in true is a redundant variable for the deep learning cases
            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, true, featureVector);
        }
        return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitEventDumpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridSize", m_gridSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridDimensions", m_gridDimensions));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
