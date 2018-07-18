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
            const CartesianVector &position(pTargetCaloHit->GetPositionVector());
            myfile << pMCParticle->GetParticleId() << "," << position.GetX() << "," << position.GetY() << "," << position.GetZ() << "," << pTargetCaloHit->GetInputEnergy() << std::endl;
        }
        catch(...)
        {
            continue;
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
