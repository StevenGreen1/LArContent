/**
 *  @file   larpandoracontent/LArMonitoring/MetricsAlgorithm.cc
 * 
 *  @brief  Implementation of the metrics algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

//#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArMonitoring/MetricsAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MetricsAlgorithm::MetricsAlgorithm() : 
    m_fileName("PerformanceAnalysis.root"),
    m_treeName("PerformanceAnalysisTree"),
    m_mcParticleListName("Input"),
    m_pfoListName("RecreatedPfos")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MetricsAlgorithm::~MetricsAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MetricsAlgorithm::Run()
{
//    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoEigenvalueBeamSlice", pSlicePropertiesBeam->m_minorAxisTwoEigenvalue));

    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    std::cout << "pPfoList->size() " << pPfoList->size() << std::endl;
    std::cout << "pMCParticleList->size() " << pMCParticleList->size() << std::endl;

    for (const MCParticle *pMCParticle : *pMCParticleList)
    {
        const LArMCParticle *pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));

        if (2001 == pLArMCParticle->GetNuanceCode())
        {
            std::cout << "Found target beam particle " << pLArMCParticle->GetParticleId() << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MetricsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_fileName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
