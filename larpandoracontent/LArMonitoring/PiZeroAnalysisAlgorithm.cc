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
    m_eventNumber(0)
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

            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            if (allMCParticleToHitsMap.find(pHitParticle) == allMCParticleToHitsMap.end())
            {
                CaloHitList caloHits;
                allMCParticleToHitsMap.insert(LArMCParticleHelper::MCContributionMap::value_type(pHitParticle, caloHits));
            }

            allMCParticleToHitsMap.at(pHitParticle).push_back(pCaloHit);
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
/*
        for (const Pfo *const pPfo : allConnectedPfos)
        {
        }
*/
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
//std::cout << "PiZeroAnalysisAlgorithm::FillAnalysisInfo" << std::endl;
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
                    if (pCaloHit->GetHitType() == TPC_VIEW_U)
                    {
                        allHitsPhoton1U.push_back(pCaloHit);
                    }
                    else if (pCaloHit->GetHitType() == TPC_VIEW_V)
                    {
                        allHitsPhoton1V.push_back(pCaloHit);
                    }
                    else if (pCaloHit->GetHitType() == TPC_VIEW_W)
                    {
                        allHitsPhoton1W.push_back(pCaloHit);
                    }

                    if (pair.first == matchedParticle1.GetPfo())
                    {
                        if (pCaloHit->GetHitType() == TPC_VIEW_U)
                        {
                            sharedHitsPhoton1U.push_back(pCaloHit);
                        }
                        else if (pCaloHit->GetHitType() == TPC_VIEW_V)
                        {
                            sharedHitsPhoton1V.push_back(pCaloHit);
                        }
                        else if (pCaloHit->GetHitType() == TPC_VIEW_W)
                        {
                            sharedHitsPhoton1W.push_back(pCaloHit);
                        }
                    }
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
                    if (pCaloHit->GetHitType() == TPC_VIEW_U)
                    {
                        allHitsPhoton2U.push_back(pCaloHit);
                    }
                    else if (pCaloHit->GetHitType() == TPC_VIEW_V)
                    {
                        allHitsPhoton2V.push_back(pCaloHit);
                    }
                    else if (pCaloHit->GetHitType() == TPC_VIEW_W)
                    {
                        allHitsPhoton2W.push_back(pCaloHit);
                    }

                    if (pair.first == matchedParticle2.GetPfo())
                    {
                        if (pCaloHit->GetHitType() == TPC_VIEW_U)
                        {
                            sharedHitsPhoton2U.push_back(pCaloHit);
                        }
                        else if (pCaloHit->GetHitType() == TPC_VIEW_V)
                        {
                            sharedHitsPhoton2V.push_back(pCaloHit);
                        }
                        else if (pCaloHit->GetHitType() == TPC_VIEW_W)
                        {
                            sharedHitsPhoton2W.push_back(pCaloHit);
                        }
                    }
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

PiZeroAnalysisAlgorithm::MatchedParticle PiZeroAnalysisAlgorithm::FillMatchedParticleInfo(const MCParticle *const pMCParticle, LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
    LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const
{
//std::cout << "PiZeroAnalysisAlgorithm::FillMatchedParticleInfo" << std::endl;
    if (mcParticleToHitsMap.find(pMCParticle) == mcParticleToHitsMap.end())
    {
        std::cout << "Missing MC particle in map" << std::endl;
    }
    const int nMCHits(mcParticleToHitsMap.at(pMCParticle).size());

    const Pfo *pBestMatch(nullptr);
    int nPfoHits(0);

    if (mcParticleToPfoHitSharingMap.find(pMCParticle) == mcParticleToPfoHitSharingMap.end())
    {
        std::cout << "Missing MC particle sharing in map" << std::endl;
    }

    for (const auto &pair : mcParticleToPfoHitSharingMap.at(pMCParticle))
    {
        if (pair.second.size() > nPfoHits)
        {
            nPfoHits = pair.second.size();
            pBestMatch = pair.first;
        }
    }

    const int nSharedHits(pfoToHitsMap.at(pBestMatch).size());
    return MatchedParticle(pMCParticle, pBestMatch, nMCHits, nPfoHits, nSharedHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PiZeroAnalysisAlgorithm::WriteToTree(AnalysisInfoVector &analysisInfoVector) const
{
    for (const auto &analysisInfo : analysisInfoVector)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber - 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsPhoton1", analysisInfo.GetMatch1().GetNMCHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsPhoton1", analysisInfo.GetMatch1().GetNPfoHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsPhoton1", analysisInfo.GetMatch1().GetSharedHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1EnergyMC", analysisInfo.GetMatch1().GetMCParticle()->GetEnergy()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PxMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PyMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PzMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon1PMC", analysisInfo.GetMatch1().GetMCParticle()->GetMomentum().GetMagnitude()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsPhoton2", analysisInfo.GetMatch2().GetNMCHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfoHitsPhoton2", analysisInfo.GetMatch2().GetNPfoHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SharedHitsPhoton2", analysisInfo.GetMatch2().GetSharedHits()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PxMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PyMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PzMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Photon2PMC", analysisInfo.GetMatch2().GetMCParticle()->GetMomentum().GetMagnitude()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroMassMC", analysisInfo.GetPiZeroMassMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroEnergyMC", analysisInfo.GetPiZeroEnergyMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPxMC", analysisInfo.GetPiZeroPxMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPyMC", analysisInfo.GetPiZeroPyMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPzMC", analysisInfo.GetPiZeroPzMC()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PiZeroPMC", analysisInfo.GetPiZeroPMC()));

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

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PiZeroAnalysisAlgorithm::MatchedParticle::MatchedParticle(const MCParticle *pMCParticle, const Pfo *pPfo, const int nMCHits, const int nPfoHits, const int nSharedHits) :
    m_pMCParticle(pMCParticle),
    m_pMatchedPfo(pPfo),
    m_nMCHits(nMCHits),
    m_nPfoHits(nPfoHits),
    m_nSharedHits(nSharedHits)
{
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
    m_piZeroMassReco(-1.f)
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
