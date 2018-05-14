/**
 *  @file   larpandoracontent/LArMonitoring/PiZeroAnalysisAlgorithm.cc
 *
 *  @brief  Implementation of the pi zero analysis algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/PiZeroAnalysisAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

PiZeroAnalysisAlgorithm::PiZeroAnalysisAlgorithm() :
    m_printAllToScreen(false),
    m_printMatchingToScreen(true),
    m_writeToTree(false),
    m_visualizePiZero(false),
    m_eventNumber(0),
    m_hitToGeV(0.00083333f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PiZeroAnalysisAlgorithm::~PiZeroAnalysisAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "PiZeroAnalysisAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PiZeroAnalysisAlgorithm::Run()
{
std::cout << "PiZeroAnalysisAlgorithm::Run" << std::endl;
    ++m_eventNumber;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;

    if (pMCParticleList && pCaloHitList)
    {
        for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            if (pCaloHit->GetHitType() != TPC_VIEW_U && pCaloHit->GetHitType() != TPC_VIEW_V && pCaloHit->GetHitType() != TPC_VIEW_W)
                continue;

            try
            {
                const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

                if (allMCParticleToHitsMap.find(pHitParticle) == allMCParticleToHitsMap.end())
                {
                    CaloHitList caloHits;
                    allMCParticleToHitsMap.insert(LArMCParticleHelper::MCContributionMap::value_type(pHitParticle, caloHits));
                }

                allMCParticleToHitsMap.at(pHitParticle).push_back(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
//                std::cout << "Hit has no MC information" << std::endl;
            }
        }
    }
/*
for (const auto &iter : allMCParticleToHitsMap)
{
std::cout << "Particle " << iter.first->GetParticleId() << ", nHits " << iter.second.size() << std::endl;
}
*/
    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;

    if (pPfoList)
    {
        PfoList allConnectedPfos;
        LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(allConnectedPfos, allMCParticleToHitsMap, pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {allMCParticleToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);

    AnalysisInfoVector analysisInfoVector;

    this->FillAnalysisInfo(pMCParticleList, allMCParticleToHitsMap, pfoToHitsMap, mcToPfoHitSharingMap, analysisInfoVector);

    if (m_writeToTree)
        this->WriteToTree(analysisInfoVector);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::FillAnalysisInfo(const MCParticleList *const pMCParticleList, LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap,
    AnalysisInfoVector &analysisInfoVector) const
{
std::cout << "PiZeroAnalysisAlgorithm::FillAnalysisInfo" << std::endl;
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (pMCParticle->GetParticleId() != 111) continue;

        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
        const int nuance(LArMCParticleHelper::GetNuanceCode(pParentMCParticle));
        if (nuance != 2001) continue;

        const MCParticleList &daughterMCParticles(pMCParticle->GetDaughterList());

        if (daughterMCParticles.size() != 2)
        {
            std::cout << "PiZero not decaying to two daughters" << std::endl;
            continue;
        }

        const MCParticle *pDaughterMCParticle1(nullptr);
        const MCParticle *pDaughterMCParticle2(nullptr);

        for (const MCParticle *const pDaughterMCParticle : daughterMCParticles)
        {
            if (!pDaughterMCParticle1)
            {
                pDaughterMCParticle1 = pDaughterMCParticle;
            }
            else
            {
                pDaughterMCParticle2 = pDaughterMCParticle;
            }
        }

        if (pDaughterMCParticle1->GetParticleId() != 22 || pDaughterMCParticle2->GetParticleId() != 22)
        {
            std::cout << "Daughter of PiZero not a photon" << std::endl;
            continue;
        }

        MatchedParticle matchedParticle1(this->FillMatchedParticleInfo(pDaughterMCParticle1, mcParticleToHitsMap, pfoToHitsMap, mcParticleToPfoHitSharingMap));
        MatchedParticle matchedParticle2(this->FillMatchedParticleInfo(pDaughterMCParticle2, mcParticleToHitsMap, pfoToHitsMap, mcParticleToPfoHitSharingMap));

        if (m_visualizePiZero)
        {
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

            PfoList photon1;
            photon1.push_back(matchedParticle1.GetPfo());
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &photon1, "Photon1Pfo", CYAN, true, true));

            CaloHitList sharedHitsPhoton1U, sharedHitsPhoton1V, sharedHitsPhoton1W;
            CaloHitList allHitsPhoton1U, allHitsPhoton1V, allHitsPhoton1W;

            for (const auto &pair : mcParticleToPfoHitSharingMap.at(pDaughterMCParticle1))
            {
                for (const CaloHit *pCaloHit: pair.second)
                {
                    this->AddCaloHit(allHitsPhoton1U, allHitsPhoton1V, allHitsPhoton1W, pCaloHit);

                    if (pair.first == matchedParticle1.GetPfo())
                        this->AddCaloHit(sharedHitsPhoton1U, sharedHitsPhoton1V, sharedHitsPhoton1W, pCaloHit);
                }
            }

            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &sharedHitsPhoton1U, "Photon1SharedCaloHitsU", BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &sharedHitsPhoton1V, "Photon1SharedCaloHitsV", BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &sharedHitsPhoton1W, "Photon1SharedCaloHitsW", BLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allHitsPhoton1U, "Photon1AllCaloHitsU", DARKBLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allHitsPhoton1V, "Photon1AllCaloHitsV", DARKBLUE));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allHitsPhoton1W, "Photon1AllCaloHitsW", DARKBLUE));

            PfoList photon2;
            photon2.push_back(matchedParticle2.GetPfo());
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &photon2, "Photon2Pfo", MAGENTA, true, true));

            CaloHitList sharedHitsPhoton2U, sharedHitsPhoton2V, sharedHitsPhoton2W;
            CaloHitList allHitsPhoton2U, allHitsPhoton2V, allHitsPhoton2W;

            for (const auto &pair : mcParticleToPfoHitSharingMap.at(pDaughterMCParticle2))
            {
                for (const CaloHit *pCaloHit: pair.second)
                {
                    this->AddCaloHit(allHitsPhoton2U, allHitsPhoton2V, allHitsPhoton2W, pCaloHit);

                    if (pair.first == matchedParticle2.GetPfo())
                        this->AddCaloHit(sharedHitsPhoton2U, sharedHitsPhoton2V, sharedHitsPhoton2W, pCaloHit);
                }
            }

            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &sharedHitsPhoton2U, "Photon2SharedCaloHitsU", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &sharedHitsPhoton2V, "Photon2SharedCaloHitsV", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &sharedHitsPhoton2W, "Photon2SharedCaloHitsW", RED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allHitsPhoton2U, "Photon2AllCaloHitsU", DARKRED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allHitsPhoton2V, "Photon2AllCaloHitsV", DARKRED));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &allHitsPhoton2W, "Photon2AllCaloHitsW", DARKRED));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        AnalysisInfo analysisInfo(matchedParticle1, matchedParticle2);
        analysisInfo.CalculatePiZeroMasses();
        analysisInfoVector.push_back(analysisInfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::AddCaloHit(CaloHitList &caloHitListU, CaloHitList &caloHitListV, CaloHitList &caloHitListW, const CaloHit *pCaloHit) const
{
    if (pCaloHit->GetHitType() == TPC_VIEW_U)
        caloHitListU.push_back(pCaloHit);

    if (pCaloHit->GetHitType() == TPC_VIEW_V)
        caloHitListV.push_back(pCaloHit);

    if (pCaloHit->GetHitType() == TPC_VIEW_W)
        caloHitListW.push_back(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PiZeroAnalysisAlgorithm::MatchedParticle PiZeroAnalysisAlgorithm::FillMatchedParticleInfo(const MCParticle *const pMCParticle, LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const
{
    if (mcParticleToHitsMap.find(pMCParticle) == mcParticleToHitsMap.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CaloHitList mcCaloHitList(mcParticleToHitsMap.at(pMCParticle));
    const Pfo *pBestMatch(nullptr);
    int nSharedHits(0);
    CaloHitList sharedCaloHitList;

    if (mcParticleToPfoHitSharingMap.find(pMCParticle) == mcParticleToPfoHitSharingMap.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    for (const auto &pair : mcParticleToPfoHitSharingMap.at(pMCParticle))
    {
        if (pair.second.size() > nSharedHits)
        {
            sharedCaloHitList = pair.second;
            nSharedHits = pair.second.size();
            pBestMatch = pair.first;
        }
    }

    CaloHitList pfoCaloHitList(pfoToHitsMap.at(pBestMatch));
    return MatchedParticle(pMCParticle, pBestMatch, mcCaloHitList, pfoCaloHitList, sharedCaloHitList, m_hitToGeV);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::WriteToTree(AnalysisInfoVector &analysisInfoVector) const
{
    for (const auto &analysisInfo : analysisInfoVector)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber - 1));

        // Photon 1
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsPhoton1", analysisInfo.GetMatch1().GetNMCHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsUPhoton1", analysisInfo.GetMatch1().GetNMCHitsU()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsVPhoton1", analysisInfo.GetMatch1().GetNMCHitsV()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsWPhoton1", analysisInfo.GetMatch1().GetNMCHitsW()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsPhoton1", analysisInfo.GetMatch1().GetNPfoHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsUPhoton1", analysisInfo.GetMatch1().GetNPfoHitsU()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsVPhoton1", analysisInfo.GetMatch1().GetNPfoHitsV()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsWPhoton1", analysisInfo.GetMatch1().GetNPfoHitsW()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsPhoton1", analysisInfo.GetMatch1().GetSharedHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsUPhoton1", analysisInfo.GetMatch1().GetSharedHitsU()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsVPhoton1", analysisInfo.GetMatch1().GetSharedHitsV()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsWPhoton1", analysisInfo.GetMatch1().GetSharedHitsW()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1EnergyMC", analysisInfo.GetMatch1().GetMCParticle()->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PxMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PyMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PzMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetMagnitude()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1EnergyCheatedPatRec", analysisInfo.GetMatch1().GetCheatedPatRecEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PxCheatedPatRec", analysisInfo.GetMatch1().GetCheatedPatRecPx()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PyCheatedPatRec", analysisInfo.GetMatch1().GetCheatedPatRecPy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PzCheatedPatRec", analysisInfo.GetMatch1().GetCheatedPatRecPz()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PxCheatedPatRecRecoDir", analysisInfo.GetMatch1().GetCheatedPatRecRecoDirPx()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PyCheatedPatRecRecoDir", analysisInfo.GetMatch1().GetCheatedPatRecRecoDirPy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PzCheatedPatRecRecoDir", analysisInfo.GetMatch1().GetCheatedPatRecRecoDirPz()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1EnergyReco", analysisInfo.GetMatch1().GetRecoEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PxReco", analysisInfo.GetMatch1().GetRecoPx()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PyReco", analysisInfo.GetMatch1().GetRecoPy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PzReco", analysisInfo.GetMatch1().GetRecoPz()));

        // Photon 2
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsPhoton2", analysisInfo.GetMatch2().GetNMCHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsUPhoton2", analysisInfo.GetMatch2().GetNMCHitsU()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsVPhoton2", analysisInfo.GetMatch2().GetNMCHitsV()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsWPhoton2", analysisInfo.GetMatch2().GetNMCHitsW()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsPhoton2", analysisInfo.GetMatch2().GetNPfoHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsUPhoton2", analysisInfo.GetMatch2().GetNPfoHitsU()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsVPhoton2", analysisInfo.GetMatch2().GetNPfoHitsV()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsWPhoton2", analysisInfo.GetMatch2().GetNPfoHitsW()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsPhoton2", analysisInfo.GetMatch2().GetSharedHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsUPhoton2", analysisInfo.GetMatch2().GetSharedHitsU()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsVPhoton2", analysisInfo.GetMatch2().GetSharedHitsV()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsWPhoton2", analysisInfo.GetMatch2().GetSharedHitsW()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2EnergyMC", analysisInfo.GetMatch2().GetMCParticle()->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PxMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PyMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PzMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetZ()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2EnergyCheatedPatRec", analysisInfo.GetMatch2().GetCheatedPatRecEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PxCheatedPatRec", analysisInfo.GetMatch2().GetCheatedPatRecPx()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PyCheatedPatRec", analysisInfo.GetMatch2().GetCheatedPatRecPy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PzCheatedPatRec", analysisInfo.GetMatch2().GetCheatedPatRecPz()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PxCheatedPatRecRecoDir", analysisInfo.GetMatch2().GetCheatedPatRecRecoDirPx()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PyCheatedPatRecRecoDir", analysisInfo.GetMatch2().GetCheatedPatRecRecoDirPy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PzCheatedPatRecRecoDir", analysisInfo.GetMatch2().GetCheatedPatRecRecoDirPz()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2EnergyReco", analysisInfo.GetMatch2().GetRecoEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PxReco", analysisInfo.GetMatch2().GetRecoPx()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PyReco", analysisInfo.GetMatch2().GetRecoPy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PzReco", analysisInfo.GetMatch2().GetRecoPz()));

        // Pion
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroMassMC", analysisInfo.GetPiZeroMassMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroEnergyMC", analysisInfo.GetPiZeroEnergyMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPxMC", analysisInfo.GetPiZeroPxMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPyMC", analysisInfo.GetPiZeroPyMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPzMC", analysisInfo.GetPiZeroPzMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPMC", analysisInfo.GetPiZeroPMC()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroMassCheatedPatRec", analysisInfo.GetPiZeroMassCheatedPatRec()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroEnergyCheatedPatRec", analysisInfo.GetPiZeroEnergyCheatedPatRec()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPxCheatedPatRec", analysisInfo.GetPiZeroPxCheatedPatRec()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPyCheatedPatRec", analysisInfo.GetPiZeroPyCheatedPatRec()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPzCheatedPatRec", analysisInfo.GetPiZeroPzCheatedPatRec()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPCheatedPatRec", analysisInfo.GetPiZeroPCheatedPatRec()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroMassCheatedPatRecRecoDir", analysisInfo.GetPiZeroMassCheatedPatRecRecoDir()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroEnergyCheatedPatRecRecoDir", analysisInfo.GetPiZeroEnergyCheatedPatRecRecoDir()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPxCheatedPatRecRecoDir", analysisInfo.GetPiZeroPxCheatedPatRecRecoDir()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPyCheatedPatRecRecoDir", analysisInfo.GetPiZeroPyCheatedPatRecRecoDir()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPzCheatedPatRecRecoDir", analysisInfo.GetPiZeroPzCheatedPatRecRecoDir()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPCheatedPatRecRecoDir", analysisInfo.GetPiZeroPCheatedPatRecRecoDir()));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroMassReco", analysisInfo.GetPiZeroMassReco()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroEnergyReco", analysisInfo.GetPiZeroEnergyReco()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPxReco", analysisInfo.GetPiZeroPxReco()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPyReco", analysisInfo.GetPiZeroPyReco()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPzReco", analysisInfo.GetPiZeroPzReco()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPReco", analysisInfo.GetPiZeroPReco()));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PiZeroAnalysisAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizePiZero", m_visualizePiZero));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitToGeVCalibration", m_hitToGeV));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PiZeroAnalysisAlgorithm::MatchedParticle::MatchedParticle(const MCParticle *pMCParticle, const Pfo *pPfo, CaloHitList &mcCaloHitList, CaloHitList &pfoCaloHitList, CaloHitList &sharedCaloHitList, const float hitToGeV) :
    m_pMCParticle(pMCParticle),
    m_pMatchedPfo(pPfo),
    m_nMCHits(mcCaloHitList.size()),
    m_nMCHitsU(0),
    m_nMCHitsV(0),
    m_nMCHitsW(0),
    m_nPfoHits(pfoCaloHitList.size()),
    m_nPfoHitsU(0),
    m_nPfoHitsV(0),
    m_nPfoHitsW(0),
    m_nSharedHits(sharedCaloHitList.size()),
    m_nSharedHitsU(0),
    m_nSharedHitsV(0),
    m_nSharedHitsW(0),
    m_recoEnergy(std::numeric_limits<float>::max()),
    m_recoPx(std::numeric_limits<float>::max()),
    m_recoPy(std::numeric_limits<float>::max()),
    m_recoPz(std::numeric_limits<float>::max()),
    m_cheatedPatRecEnergy(std::numeric_limits<float>::max()),
    m_cheatedPatRecPx(std::numeric_limits<float>::max()),
    m_cheatedPatRecPy(std::numeric_limits<float>::max()),
    m_cheatedPatRecPz(std::numeric_limits<float>::max()),
    m_cheatedPatRecRecoDirPx(std::numeric_limits<float>::max()),
    m_cheatedPatRecRecoDirPy(std::numeric_limits<float>::max()),
    m_cheatedPatRecRecoDirPz(std::numeric_limits<float>::max()),
    m_hitToGeV(hitToGeV)
{
    this->CountHits(mcCaloHitList, m_nMCHitsU, m_nMCHitsV, m_nMCHitsW);
    this->CountHits(pfoCaloHitList, m_nPfoHitsU, m_nPfoHitsV, m_nPfoHitsW);
    this->CountHits(sharedCaloHitList, m_nSharedHitsU, m_nSharedHitsV, m_nSharedHitsW);
    this->CalculateRecoEnergy(mcCaloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::MatchedParticle::CountHits(CaloHitList &caloHitList, int &nHitsU, int &nHitsV, int &nHitsW)
{
    for (const CaloHit *pCaloHit : caloHitList)
    {
        if (pCaloHit->GetHitType() == TPC_VIEW_U)
            nHitsU++;

        if (pCaloHit->GetHitType() == TPC_VIEW_V)
            nHitsV++;

        if (pCaloHit->GetHitType() == TPC_VIEW_W)
            nHitsW++;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::MatchedParticle::CalculateRecoEnergy(CaloHitList &mcCaloHitList)
{
    // ATTN: Assume collection plane W for energy, could be better
    // ATTN2: Also assuming photons here so no masses
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(m_pMatchedPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(m_pMatchedPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(m_pMatchedPfo, TPC_VIEW_W, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(m_pMatchedPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(m_pMatchedPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetIsolatedCaloHits(m_pMatchedPfo, TPC_VIEW_W, caloHitList);

    // Reco Energy
    m_recoEnergy = 0.f;
    m_recoEnergy = (float)(caloHitList.size()) * m_hitToGeV;

    // Cheated Energy
    m_cheatedPatRecEnergy = 0.f;
    m_cheatedPatRecEnergy = (float)(mcCaloHitList.size()) * m_hitToGeV;

    // Cheated Direction
    m_cheatedPatRecPx = m_pMCParticle->GetMomentum().GetUnitVector().GetX() * m_cheatedPatRecEnergy;
    m_cheatedPatRecPy = m_pMCParticle->GetMomentum().GetUnitVector().GetY() * m_cheatedPatRecEnergy;
    m_cheatedPatRecPz = m_pMCParticle->GetMomentum().GetUnitVector().GetZ() * m_cheatedPatRecEnergy;

    // Reco Direction
    CaloHitList caloHitList3D;
    LArPfoHelper::GetCaloHits(m_pMatchedPfo, TPC_3D, caloHitList3D);
    LArPfoHelper::GetIsolatedCaloHits(m_pMatchedPfo, TPC_3D, caloHitList3D);

    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);

    try
    {
        LArPcaHelper::RunPca(caloHitList3D, centroid, eigenValues, eigenVecs);
        CartesianVector direction(eigenVecs.at(0));
        m_recoPx = direction.GetUnitVector().GetX() * m_recoEnergy;
        m_recoPy = direction.GetUnitVector().GetY() * m_recoEnergy;
        m_recoPz = direction.GetUnitVector().GetZ() * m_recoEnergy;

        m_cheatedPatRecRecoDirPx = direction.GetUnitVector().GetX() * m_cheatedPatRecEnergy;
        m_cheatedPatRecRecoDirPy = direction.GetUnitVector().GetY() * m_cheatedPatRecEnergy;
        m_cheatedPatRecRecoDirPz = direction.GetUnitVector().GetZ() * m_cheatedPatRecEnergy;
    }
    catch (const StatusCodeException &) {}
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PiZeroAnalysisAlgorithm::AnalysisInfo::AnalysisInfo(MatchedParticle &photon1, MatchedParticle &photon2) :
    m_photon1(photon1),
    m_photon2(photon2),
    m_piZeroMassMC(-1.f),
    m_piZeroEnergyMC(std::numeric_limits<float>::max()),
    m_piZeroPxMC(std::numeric_limits<float>::max()),
    m_piZeroPyMC(std::numeric_limits<float>::max()),
    m_piZeroPzMC(std::numeric_limits<float>::max()),
    m_piZeroPMC(std::numeric_limits<float>::max()),
    m_piZeroMassCheatedPatRec(-1.f),
    m_piZeroEnergyCheatedPatRec(std::numeric_limits<float>::max()),
    m_piZeroPxCheatedPatRec(std::numeric_limits<float>::max()),
    m_piZeroPyCheatedPatRec(std::numeric_limits<float>::max()),
    m_piZeroPzCheatedPatRec(std::numeric_limits<float>::max()),
    m_piZeroPCheatedPatRec(std::numeric_limits<float>::max()),
    m_piZeroMassCheatedPatRecRecoDir(-1.f),
    m_piZeroEnergyCheatedPatRecRecoDir(std::numeric_limits<float>::max()),
    m_piZeroPxCheatedPatRecRecoDir(std::numeric_limits<float>::max()),
    m_piZeroPyCheatedPatRecRecoDir(std::numeric_limits<float>::max()),
    m_piZeroPzCheatedPatRecRecoDir(std::numeric_limits<float>::max()),
    m_piZeroPCheatedPatRecRecoDir(std::numeric_limits<float>::max()),
    m_piZeroMassReco(-1.f),
    m_piZeroEnergyReco(std::numeric_limits<float>::max()),
    m_piZeroPxReco(std::numeric_limits<float>::max()),
    m_piZeroPyReco(std::numeric_limits<float>::max()),
    m_piZeroPzReco(std::numeric_limits<float>::max()),
    m_piZeroPReco(std::numeric_limits<float>::max())
{
    if (photon1.GetNMCHits() > photon2.GetNMCHits())
    {
        m_photon1 = photon1;
        m_photon2 = photon2;
    }
    else
    {
        m_photon1 = photon2;
        m_photon2 = photon1;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::AnalysisInfo::CalculatePiZeroMasses()
{
    // MC
    CartesianVector momentum1(m_photon1.GetMCParticle()->GetMomentum());
    const float energy1(m_photon1.GetMCParticle()->GetEnergy());

    CartesianVector momentum2(m_photon2.GetMCParticle()->GetMomentum());
    const float energy2(m_photon2.GetMCParticle()->GetEnergy());

    CartesianVector piZeroMomentum(momentum1 + momentum2);

    m_piZeroEnergyMC = energy1 + energy2;
    m_piZeroPxMC = piZeroMomentum.GetX();
    m_piZeroPyMC = piZeroMomentum.GetY();
    m_piZeroPzMC = piZeroMomentum.GetZ();
    m_piZeroPMC = piZeroMomentum.GetMagnitude();

    if (m_piZeroEnergyMC > m_piZeroPMC)
        m_piZeroMassMC = std::sqrt(m_piZeroEnergyMC*m_piZeroEnergyMC - m_piZeroPMC*m_piZeroPMC);

    // Reco E, Cheated Pat Rec, Cheated Direction
    CartesianVector cheatedPatRecMomentum1(m_photon1.GetCheatedPatRecPx(), m_photon1.GetCheatedPatRecPy(), m_photon1.GetCheatedPatRecPz());
    const float cheatedPatRecEnergy1(m_photon1.GetCheatedPatRecEnergy());

    CartesianVector cheatedPatRecMomentum2(m_photon2.GetCheatedPatRecPx(), m_photon2.GetCheatedPatRecPy(), m_photon2.GetCheatedPatRecPz());
    const float cheatedPatRecEnergy2(m_photon2.GetCheatedPatRecEnergy());

    CartesianVector cheatedPatRecPiZeroMomentum(cheatedPatRecMomentum1 + cheatedPatRecMomentum2);

    m_piZeroEnergyCheatedPatRec = cheatedPatRecEnergy1 + cheatedPatRecEnergy2;
    m_piZeroPxCheatedPatRec = cheatedPatRecPiZeroMomentum.GetX();
    m_piZeroPyCheatedPatRec = cheatedPatRecPiZeroMomentum.GetY();
    m_piZeroPzCheatedPatRec = cheatedPatRecPiZeroMomentum.GetZ();
    m_piZeroPCheatedPatRec = cheatedPatRecPiZeroMomentum.GetMagnitude();

    if (m_piZeroEnergyCheatedPatRec > m_piZeroPCheatedPatRec)
        m_piZeroMassCheatedPatRec = std::sqrt(m_piZeroEnergyCheatedPatRec*m_piZeroEnergyCheatedPatRec - m_piZeroPCheatedPatRec*m_piZeroPCheatedPatRec);

    // Reco E, Cheated Pat Rec, Reco Direction
    CartesianVector cheatedPatRecRecoDirMomentum1(m_photon1.GetCheatedPatRecRecoDirPx(), m_photon1.GetCheatedPatRecRecoDirPy(), m_photon1.GetCheatedPatRecRecoDirPz());
    const float cheatedPatRecRecoDirEnergy1(m_photon1.GetCheatedPatRecEnergy());

    CartesianVector cheatedPatRecRecoDirMomentum2(m_photon2.GetCheatedPatRecRecoDirPx(), m_photon2.GetCheatedPatRecRecoDirPy(), m_photon2.GetCheatedPatRecRecoDirPz());
    const float cheatedPatRecRecoDirEnergy2(m_photon2.GetCheatedPatRecEnergy());

    CartesianVector cheatedPatRecRecoDirPiZeroMomentum(cheatedPatRecRecoDirMomentum1 + cheatedPatRecRecoDirMomentum2);

    m_piZeroEnergyCheatedPatRecRecoDir = cheatedPatRecRecoDirEnergy1 + cheatedPatRecRecoDirEnergy2;
    m_piZeroPxCheatedPatRecRecoDir = cheatedPatRecRecoDirPiZeroMomentum.GetX();
    m_piZeroPyCheatedPatRecRecoDir = cheatedPatRecRecoDirPiZeroMomentum.GetY();
    m_piZeroPzCheatedPatRecRecoDir = cheatedPatRecRecoDirPiZeroMomentum.GetZ();
    m_piZeroPCheatedPatRecRecoDir = cheatedPatRecRecoDirPiZeroMomentum.GetMagnitude();

    if (m_piZeroEnergyCheatedPatRecRecoDir > m_piZeroPCheatedPatRecRecoDir)
        m_piZeroMassCheatedPatRecRecoDir = std::sqrt(m_piZeroEnergyCheatedPatRecRecoDir*m_piZeroEnergyCheatedPatRecRecoDir - m_piZeroPCheatedPatRecRecoDir*m_piZeroPCheatedPatRecRecoDir);

    // Reco E + Reco MC
    CartesianVector recoMomentum1(m_photon1.GetRecoPx(), m_photon1.GetRecoPy(), m_photon1.GetRecoPz());
    const float recoEnergy1(m_photon1.GetRecoEnergy());

    CartesianVector recoMomentum2(m_photon2.GetRecoPx(), m_photon2.GetRecoPy(), m_photon2.GetRecoPz());
    const float recoEnergy2(m_photon2.GetRecoEnergy());

    CartesianVector recoPiZeroMomentum(recoMomentum1 + recoMomentum2);

    m_piZeroEnergyReco = recoEnergy1 + recoEnergy2;
    m_piZeroPxReco = recoPiZeroMomentum.GetX();
    m_piZeroPyReco = recoPiZeroMomentum.GetY();
    m_piZeroPzReco = recoPiZeroMomentum.GetZ();
    m_piZeroPReco = recoPiZeroMomentum.GetMagnitude();

    if (m_piZeroEnergyReco > m_piZeroPReco)
        m_piZeroMassReco = std::sqrt(m_piZeroEnergyReco*m_piZeroEnergyReco - m_piZeroPReco*m_piZeroPReco);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
