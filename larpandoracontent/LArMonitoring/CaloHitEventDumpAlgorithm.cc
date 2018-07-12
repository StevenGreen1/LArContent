/**
 *  @file   larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CaloHitEventDumpAlgorithm::CaloHitEventDumpAlgorithm() :
    m_textFileName("CaloHitEventDump.txt"),
    m_gridSize(28),
    m_gridDimensions(50)
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

    std::ofstream myfile;
    myfile.open (m_textFileName.c_str(), std::ios_base::app);

    for (const CaloHit *pTargetCaloHit : wCaloHitList)
    {
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pTargetCaloHit));
            myfile << "CaloHit," << pMCParticle->GetParticleId() << std::endl;
        }
        catch(...)
        {
            continue;
        }

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
                const float globalBinNumber(xBin*twoDHistogram.GetNBinsX() + yBin);
                const float content(twoDHistogram.GetBinContent(xBin, yBin));
                myfile << globalBinNumber << "," << content << std::endl;
            }
        }
    }

    myfile.close();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitEventDumpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TextFileName", m_textFileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "GridSize", m_gridSize));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "GridDimensions", m_gridDimensions));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
