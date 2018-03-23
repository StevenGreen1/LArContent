/**
 *  @file   larpandoracontent/LArMonitoring/BeamTargetsAlgorithm.cc
 * 
 *  @brief  Implementation of the beam target algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArMonitoring/BeamTargetsAlgorithm.h"

namespace lar_content
{

using namespace pandora;

StatusCode BeamTargetsAlgorithm::Run()
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);
/*
    CaloHitList uCaloHitList, vCaloHitList, wCaloHitList;

    for (const auto caloHit : *pCaloHitList)
    {
        if (caloHit->GetHitType() == TPC_VIEW_U) uCaloHitList.push_back(caloHit);
        else if (caloHit->GetHitType() == TPC_VIEW_V) vCaloHitList.push_back(caloHit);
        else if (caloHit->GetHitType() == TPC_VIEW_W) wCaloHitList.push_back(caloHit);
    }
*/
    MCParticleList beamMCParticleList;
    this->GetBeamMCParticles(pMCParticleList, beamMCParticleList);

    // Obtain vector: primary mc particles
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);

    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToPrimaryMCMap, selectedCaloHitList, true, 1000);

    // Obtain maps: [hit -> primary mc particle], [primary mc particle -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap hitToPrimaryMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

    for (const auto iter : mcToTrueHitListMap)
    {
        const MCParticle *pMCPrimary(iter.first);
        if (!LArMCParticleHelper::IsBeamParticle(pMCPrimary))
            continue;

        const CaloHitList caloHitListMC(iter.second);

        CaloHitList uCaloHitList, vCaloHitList, wCaloHitList;
        for (const auto caloHit : caloHitListMC)
        {
            if (caloHit->GetHitType() == TPC_VIEW_U) uCaloHitList.push_back(caloHit);
            else if (caloHit->GetHitType() == TPC_VIEW_V) vCaloHitList.push_back(caloHit);
            else if (caloHit->GetHitType() == TPC_VIEW_W) wCaloHitList.push_back(caloHit);
        }

        MCParticleList mcParticle;
        mcParticle.push_back(pMCPrimary);

        LArMCParticleHelper::MCContributionMap newTargetToCaloHitMap;
        this->GetCaloHitToTargetMap(wCaloHitList, pMCPrimary, newTargetToCaloHitMap);

        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        int counter(0);

        for (const auto it : newTargetToCaloHitMap)
        {
            counter++;
            std::string name("CaloHits_Target_" + std::to_string(counter));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &it.second, name.c_str(), AUTO));
        }

        PANDORA_MONITORING_API(VisualizeMCParticles(this->GetPandora(), &mcParticle, "PrimaryMCParticle", AUTO, nullptr));
        PANDORA_MONITORING_API(VisualizeMCParticles(this->GetPandora(), pMCParticleList, "AllMCParticles", AUTO, nullptr));
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uCaloHitList, "UMCParticleCaloHits", RED));
//        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vCaloHitList, "VMCParticleCaloHits", GREEN));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wCaloHitList, "AllWPrimaryMCParticleCaloHits", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
/*
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    PANDORA_MONITORING_API(VisualizeMCParticles(this->GetPandora(), pMCParticleList, "CurrentMCParticles", AUTO, nullptr));
    PANDORA_MONITORING_API(VisualizeMCParticles(this->GetPandora(), &beamMCParticleList, "BeamMCParticles", AUTO, nullptr));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uCaloHitList, "UCurrentCaloHits", RED));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vCaloHitList, "VCurrentCaloHits", GREEN));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wCaloHitList, "WCurrentCaloHits", BLUE));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, "CurrentPfos", AUTO, true, true));
*/
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::GetBeamMCParticles(const MCParticleList *pMCParticleList, MCParticleList &beamMCParticleList)
{
    for (const auto *pMCParticle : *pMCParticleList)
    {
        const LArMCParticle *const pLArMCParticleParent(dynamic_cast<const LArMCParticle*>(LArMCParticleHelper::GetParentMCParticle(pMCParticle)));
        const int nuance(pLArMCParticleParent->GetNuanceCode());

        if (nuance == 2000 || nuance == 2001)
            beamMCParticleList.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::GetCaloHitToTargetMap(const CaloHitList &caloHitList, const MCParticle *pMCPrimary, LArMCParticleHelper::MCContributionMap &mcToCaloHitMap)
{
    MCParticleList allMCParticles;
    this->GetAllDownstreamMCParticles(pMCPrimary, allMCParticles);

    MCParticleList filteredMCParticles;
    this->FilterMCParticles(allMCParticles, filteredMCParticles, caloHitList);

    for (const auto *pCaloHit : caloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_W)
            continue;

        const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

        if (pHitParticle->GetParentList().empty())
        {
            // Target is primary particle
            this->AddCaloHitToMap(pCaloHit, pHitParticle, mcToCaloHitMap);
        }
        else
        {
            // Target not primary particles
            const MCParticle *pTarget(nullptr);
            this->IsVisibleTarget(pHitParticle, filteredMCParticles, pTarget);
            if (pTarget == nullptr)
                pTarget = pMCPrimary;
            this->AddCaloHitToMap(pCaloHit, pTarget, mcToCaloHitMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::GetAllDownstreamMCParticles(const MCParticleList &inputMCParticleList, MCParticleList &outputMCParticleList)
{
    for (const MCParticle *const pMCParticle : inputMCParticleList)
        this->GetAllDownstreamMCParticles(pMCParticle, outputMCParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::GetAllDownstreamMCParticles(const MCParticle *const pMCParticle, MCParticleList &outputMCParticleList)
{
    if (outputMCParticleList.end() != std::find(outputMCParticleList.begin(), outputMCParticleList.end(), pMCParticle))
        return;

    outputMCParticleList.push_back(pMCParticle);
    this->GetAllDownstreamMCParticles(pMCParticle->GetDaughterList(), outputMCParticleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::FilterMCParticles(const MCParticleList &allMCParticles, MCParticleList &outputMCParticles, const CaloHitList &caloHitList)
{
    MCParticleList filteredMCParticles;
    for (const MCParticle *pMCParticle : allMCParticles)
    {
        if (pMCParticle->GetParentList().empty())
            continue;

        int nCaloHits(0);

        for (const CaloHit *pCaloHit : caloHitList)
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            if (pHitParticle == pMCParticle)
                nCaloHits++;
        }

        if (LArMCParticleHelper::IsVisible(pMCParticle) && std::find(filteredMCParticles.begin(), filteredMCParticles.end(), pMCParticle) == filteredMCParticles.end() && nCaloHits > 10)
            filteredMCParticles.push_back(pMCParticle);
    }

    for (const MCParticle *pMCParticle : filteredMCParticles)
    {
//        if (!this->CheckParentInList(pMCParticle, filteredMCParticles))
        outputMCParticles.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BeamTargetsAlgorithm::CheckParentInList(const MCParticle *pMCParticle, const MCParticleList &filteredMCParticles)
{
    const MCParticleList parentList(pMCParticle->GetParentList());

    if (parentList.size() == 0)
    {
        return false;
    }
    else if (parentList.size() == 1)
    {
        const MCParticle *pParent(parentList.front());

        if (std::find(filteredMCParticles.begin(), filteredMCParticles.end(), pParent) != filteredMCParticles.end())
        {
            return true;
        }
        else
        {
            return this->CheckParentInList(pParent, filteredMCParticles);
        }
    }
    else
    {
        return false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::IsVisibleTarget(const MCParticle *pParticle, const MCParticleList &targets, const MCParticle *&pTarget)
{
    // Particle is a target
    if (std::find(targets.begin(), targets.end(), pParticle) != targets.end())
    {
        pTarget = pParticle;
        return;
    }
    // Particle not target, check parents
    else
    {
        const MCParticleList parentList(pParticle->GetParentList());
        if (parentList.empty())
        {
            pTarget = nullptr;
            return;
        }

        if (parentList.size() > 1)
        {
            pTarget = nullptr;
            return;
        }

        const MCParticle *pParent(parentList.front());
        this->IsVisibleTarget(pParent, targets, pTarget);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamTargetsAlgorithm::AddCaloHitToMap(const CaloHit *pCaloHit, const MCParticle *pTarget, LArMCParticleHelper::MCContributionMap &mcToCaloHitMap)
{
    if (pTarget == nullptr) return;

    if (mcToCaloHitMap.find(pTarget) == mcToCaloHitMap.end())
    {
        std::cout << "New MC, pdg " << pTarget->GetParticleId() << std::endl;
        CaloHitList primaryCaloHits;
        primaryCaloHits.push_back(pCaloHit);
        mcToCaloHitMap.insert(LArMCParticleHelper::MCContributionMap::value_type(pTarget, primaryCaloHits));
    }
    else
    {
        mcToCaloHitMap.at(pTarget).push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamTargetsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
