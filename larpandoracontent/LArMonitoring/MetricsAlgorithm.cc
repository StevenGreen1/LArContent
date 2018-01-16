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

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    IntVector triggeredParticleID;

    for (const MCParticle *pMCParticle : *pMCParticleList)
    {
        const LArMCParticle *pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));

        if (2001 == pLArMCParticle->GetNuanceCode())
        {
            triggeredParticleID.push_back(pLArMCParticle->GetParticleId());
            std::cout << "Found target beam particle " << pLArMCParticle->GetParticleId() << std::endl;
        }
    }

    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    int nTrackLike(0);
    int nShowerLike(0);

    for (const Pfo *const pPfo : *pPfoList)
    {
        if (12 == std::abs(pPfo->GetParticleId()))
        {
            nShowerLike++;
            std::cout << "Shower Beam Candidate Identified" << std::endl;
        }
        else if (14 == std::abs(pPfo->GetParticleId()))
        {
            nTrackLike++;
            std::cout << "Track Beam Candidate Identified" << std::endl;
        }
    }

    int nTriggeredParticles(triggeredParticleID.size());

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TriggeredParticleID", &triggeredParticleID));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTriggeredParticles", nTriggeredParticles));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nShowerLikePFOs", nShowerLike));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrackLikePFOs", nTrackLike));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

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
