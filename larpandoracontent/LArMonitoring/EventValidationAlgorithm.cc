/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/EventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

EventValidationAlgorithm::EventValidationAlgorithm() :
    m_useTrueNeutrinosOnly(false),
    m_testBeamMode(false),
    m_selectInputHits(true),
    m_minHitSharingFraction(0.9f),
    m_maxPhotonPropagation(2.5f),
    m_printAllToScreen(false),
    m_printMatchingToScreen(true),
    m_writeToTree(false),
    m_useSmallPrimaries(true),
    m_matchingMinSharedHits(5),
    m_matchingMinCompleteness(0.1f),
    m_matchingMinPurity(0.5f),
    m_fileIdentifier(0),
    m_eventNumber(0),
    m_gridSize(16),
    m_gridDimensions(50),
    m_verbose(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventValidationAlgorithm::~EventValidationAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::Run()
{
    ++m_eventNumber;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

    for (const Pfo *pPfo : allConnectedPfos)
    {
        CaloHitList caloHitList;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);
        LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_U, caloHitList);
        LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_V, caloHitList);
        LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_W, caloHitList);

        const int id(pPfo->GetParticleId());
        int label(std::numeric_limits<int>::max());

        if (id == 22 || std::abs(id) == 11)
        {
            label = 0;
        }
        else
        {
            label = 1;
        }

        for (const CaloHit *pCaloHit : caloHitList)
            m_hitToIntMap.insert(HitToIntMap::value_type(pCaloHit, label));
    }

    ValidationInfo validationInfo;
    this->FillValidationInfo(pMCParticleList, pCaloHitList, pPfoList, validationInfo);

    if (m_printAllToScreen)
        this->PrintAllMatches(validationInfo, pCaloHitList);

    if (m_printMatchingToScreen)
        this->PrintInterpretedMatches(validationInfo, pCaloHitList);

    if (m_writeToTree)
        this->WriteInterpretedMatches(validationInfo, pCaloHitList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::PrimaryParameters parameters;

        parameters.m_selectInputHits = m_selectInputHits;
        parameters.m_minHitSharingFraction = m_minHitSharingFraction;
        parameters.m_maxPhotonPropagation = m_maxPhotonPropagation;
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, allMCParticleToHitsMap);
        if (!m_useTrueNeutrinosOnly) LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {
        PfoList allConnectedPfos;
        LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

        PfoList finalStatePfos;
        for (const ParticleFlowObject *const pPfo : allConnectedPfos)
        {
            if ((!m_testBeamMode && LArPfoHelper::IsFinalState(pPfo)) || (m_testBeamMode && pPfo->GetParentPfoList().empty()))
                finalStatePfos.push_back(pPfo);
        }

        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap);
        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const CaloHitList *const pCaloHitList, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    int nNeutrinoPrimaries(0);
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary) && validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary)) ++nNeutrinoPrimaries;

    PfoVector primaryPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(validationInfo.GetPfoToHitsMap(), primaryPfoVector);

    int pfoIndex(0), neutrinoPfoIndex(0);
    PfoToIdMap pfoToIdMap, neutrinoPfoToIdMap;
    for (const Pfo *const pPrimaryPfo : primaryPfoVector)
    {
        pfoToIdMap.insert(PfoToIdMap::value_type(pPrimaryPfo, ++pfoIndex));
        const Pfo *const pRecoNeutrino(LArPfoHelper::IsNeutrinoFinalState(pPrimaryPfo) ? LArPfoHelper::GetParentNeutrino(pPrimaryPfo) : nullptr);

        if (pRecoNeutrino && !neutrinoPfoToIdMap.count(pRecoNeutrino))
            neutrinoPfoToIdMap.insert(PfoToIdMap::value_type(pRecoNeutrino, ++neutrinoPfoIndex));
    }

    PfoSet recoNeutrinos;
    MCParticleList associatedMCPrimaries;

    int nCorrectNu(0), nTotalNu(0), nCorrectTB(0), nTotalTB(0), nCorrectCR(0), nTotalCR(0), nFakeNu(0), nFakeCR(0), nSplitNu(0), nSplitCR(0), nLost(0);
    int mcPrimaryIndex(0), nTargetMatches(0), nTargetNuMatches(0), nTargetCRMatches(0), nTargetGoodNuMatches(0), nTargetNuSplits(0), nTargetNuLosses(0);
    IntVector mcPrimaryId, mcPrimaryPdg, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW, nMCTrkHitsU, nMCTrkHitsV, nMCTrkHitsW, nMCShwHitsU, nMCShwHitsV, nMCShwHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedNuPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoIsRecoNu, bestMatchPfoRecoNuId, bestMatchPfoIsTestBeam;
    IntVector bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoTrueTrkNHitsU, bestMatchPfoTrueTrkNHitsV, bestMatchPfoTrueTrkNHitsW, bestMatchPfoTrueShwNHitsU, bestMatchPfoTrueShwNHitsV, bestMatchPfoTrueShwNHitsW;
    IntVector bestMatchPfoCNNTrkNHitsU, bestMatchPfoCNNTrkNHitsV, bestMatchPfoCNNTrkNHitsW, bestMatchPfoPandoraTrkNHitsU, bestMatchPfoPandoraTrkNHitsV, bestMatchPfoPandoraTrkNHitsW;
    IntVector bestMatchPfoCNNShwNHitsU, bestMatchPfoCNNShwNHitsV, bestMatchPfoCNNShwNHitsW, bestMatchPfoPandoraShwNHitsU, bestMatchPfoPandoraShwNHitsV, bestMatchPfoPandoraShwNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;

    std::stringstream targetSS;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        const bool hasMatch(mcToPfoHitSharingMap.count(pMCPrimary) && !mcToPfoHitSharingMap.at(pMCPrimary).empty());
        const bool isTargetPrimary(validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary));

        if (!hasMatch && !isTargetPrimary)
            continue;

        associatedMCPrimaries.push_back(pMCPrimary);
        const int nTargetPrimaries(associatedMCPrimaries.size());
        const bool isLastNeutrinoPrimary(++mcPrimaryIndex == nNeutrinoPrimaries);
        const CaloHitList &mcPrimaryHitList(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary));

        const int mcNuanceCode(LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)));
        const int isBeamNeutrinoFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCPrimary));
        const int isBeamParticle(LArMCParticleHelper::IsBeamParticle(pMCPrimary));
        const int isCosmicRay(LArMCParticleHelper::IsCosmicRay(pMCPrimary));
#ifdef MONITORING
        const CartesianVector &targetVertex(LArMCParticleHelper::GetParentMCParticle(pMCPrimary)->GetVertex());
        const float targetVertexX(targetVertex.GetX()), targetVertexY(targetVertex.GetY()), targetVertexZ(targetVertex.GetZ());
#endif
        targetSS << (!isTargetPrimary ? "(Non target) " : "")
                 << "PrimaryId " << mcPrimaryIndex
                 << ", Nu " << isBeamNeutrinoFinalState
                 << ", TB " << isBeamParticle
                 << ", CR " << isCosmicRay
                 << ", MCPDG " << pMCPrimary->GetParticleId()
                 << ", Energy " << pMCPrimary->GetEnergy()
                 << ", Dist. " << (pMCPrimary->GetEndpoint() - pMCPrimary->GetVertex()).GetMagnitude()
                 << ", nMCHits " << mcPrimaryHitList.size()
                 << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList) << ")" << std::endl;

        mcPrimaryId.push_back(mcPrimaryIndex);
        mcPrimaryPdg.push_back(pMCPrimary->GetParticleId());
        mcPrimaryE.push_back(pMCPrimary->GetEnergy());
        mcPrimaryPX.push_back(pMCPrimary->GetMomentum().GetX());
        mcPrimaryPY.push_back(pMCPrimary->GetMomentum().GetY());
        mcPrimaryPZ.push_back(pMCPrimary->GetMomentum().GetZ());
        mcPrimaryVtxX.push_back(pMCPrimary->GetVertex().GetX());
        mcPrimaryVtxY.push_back(pMCPrimary->GetVertex().GetY());
        mcPrimaryVtxZ.push_back(pMCPrimary->GetVertex().GetZ());
        mcPrimaryEndX.push_back(pMCPrimary->GetEndpoint().GetX());
        mcPrimaryEndY.push_back(pMCPrimary->GetEndpoint().GetY());
        mcPrimaryEndZ.push_back(pMCPrimary->GetEndpoint().GetZ());
        nMCHitsTotal.push_back(mcPrimaryHitList.size());
        nMCHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, mcPrimaryHitList));
        nMCHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, mcPrimaryHitList));
        nMCHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, mcPrimaryHitList));

        int nTrkU(0), nTrkV(0), nTrkW(0), nShwU(0), nShwV(0), nShwW(0);
        this->CountTrkShwHitsByType(TPC_VIEW_U, mcPrimaryHitList, nTrkU, nShwU);
        this->CountTrkShwHitsByType(TPC_VIEW_V, mcPrimaryHitList, nTrkV, nShwV);
        this->CountTrkShwHitsByType(TPC_VIEW_W, mcPrimaryHitList, nTrkW, nShwW);

        nMCTrkHitsU.push_back(nTrkU);
        nMCTrkHitsV.push_back(nTrkV);
        nMCTrkHitsW.push_back(nTrkW);
        nMCShwHitsU.push_back(nShwU);
        nMCShwHitsV.push_back(nShwV);
        nMCShwHitsW.push_back(nShwW);

        int matchIndex(0), nPrimaryMatches(0), nPrimaryNuMatches(0), nPrimaryCRMatches(0), nPrimaryGoodNuMatches(0), nPrimaryNuSplits(0);
#ifdef MONITORING
        float recoVertexX(std::numeric_limits<float>::max()), recoVertexY(std::numeric_limits<float>::max()), recoVertexZ(std::numeric_limits<float>::max());
#endif
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pMCPrimary))
        {
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isRecoNeutrinoFinalState(LArPfoHelper::IsNeutrinoFinalState(pfoToSharedHits.first));
            const bool isRecoTestBeam(LArPfoHelper::IsTestBeam(pfoToSharedHits.first));
            const bool isGoodMatch(this->IsGoodMatch(mcPrimaryHitList, pfoHitList, sharedHitList));

            const int pfoId(pfoToIdMap.at(pfoToSharedHits.first));
            const int recoNuId(isRecoNeutrinoFinalState ? neutrinoPfoToIdMap.at(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first)) : -1);

            int nTrueTrkU(0), nTrueTrkV(0), nTrueTrkW(0), nTrueShwU(0), nTrueShwV(0), nTrueShwW(0);
            this->CountTrkShwHitsByType(TPC_VIEW_U, pfoHitList, nTrueTrkU, nTrueShwU);
            this->CountTrkShwHitsByType(TPC_VIEW_V, pfoHitList, nTrueTrkV, nTrueShwV);
            this->CountTrkShwHitsByType(TPC_VIEW_W, pfoHitList, nTrueTrkW, nTrueShwW);

            int nCNNTrkU(0), nCNNTrkV(0), nCNNTrkW(0), nCNNShwU(0), nCNNShwV(0), nCNNShwW(0);
            this->CountCNNTrkShwHitsByType(TPC_VIEW_U, pfoHitList, *pCaloHitList, nCNNTrkU, nCNNShwU);
            this->CountCNNTrkShwHitsByType(TPC_VIEW_V, pfoHitList, *pCaloHitList, nCNNTrkV, nCNNShwV);
            this->CountCNNTrkShwHitsByType(TPC_VIEW_W, pfoHitList, *pCaloHitList, nCNNTrkW, nCNNShwW);

            int nPandoraTrkU(0), nPandoraTrkV(0), nPandoraTrkW(0), nPandoraShwU(0), nPandoraShwV(0), nPandoraShwW(0);
            this->CountPandoraTrkShwHitsByType(TPC_VIEW_U, pfoHitList, nPandoraTrkU, nPandoraShwU);
            this->CountPandoraTrkShwHitsByType(TPC_VIEW_V, pfoHitList, nPandoraTrkV, nPandoraShwV);
            this->CountPandoraTrkShwHitsByType(TPC_VIEW_W, pfoHitList, nPandoraTrkW, nPandoraShwW);

            if (0 == matchIndex++)
            {
                bestMatchPfoId.push_back(pfoId);
                bestMatchPfoPdg.push_back(pfoToSharedHits.first->GetParticleId());
                bestMatchPfoIsRecoNu.push_back(isRecoNeutrinoFinalState ? 1 : 0);
                bestMatchPfoRecoNuId.push_back(recoNuId);
                bestMatchPfoIsTestBeam.push_back(isRecoTestBeam ? 1 : 0);
                bestMatchPfoNHitsTotal.push_back(pfoHitList.size());
                bestMatchPfoNHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                bestMatchPfoNHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                bestMatchPfoNHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                bestMatchPfoTrueTrkNHitsU.push_back(nTrueTrkU);
                bestMatchPfoTrueShwNHitsU.push_back(nTrueShwU);
                bestMatchPfoTrueTrkNHitsV.push_back(nTrueTrkV);
                bestMatchPfoTrueShwNHitsV.push_back(nTrueShwV);
                bestMatchPfoTrueTrkNHitsW.push_back(nTrueTrkW);
                bestMatchPfoTrueShwNHitsW.push_back(nTrueShwW);
                bestMatchPfoCNNTrkNHitsU.push_back(nCNNTrkU);
                bestMatchPfoCNNShwNHitsU.push_back(nCNNShwU);
                bestMatchPfoCNNTrkNHitsV.push_back(nCNNTrkV);
                bestMatchPfoCNNShwNHitsV.push_back(nCNNShwV);
                bestMatchPfoCNNTrkNHitsW.push_back(nCNNTrkW);
                bestMatchPfoCNNShwNHitsW.push_back(nCNNShwW);
                bestMatchPfoPandoraTrkNHitsU.push_back(nPandoraTrkU);
                bestMatchPfoPandoraShwNHitsU.push_back(nPandoraShwU);
                bestMatchPfoPandoraTrkNHitsV.push_back(nPandoraTrkV);
                bestMatchPfoPandoraShwNHitsV.push_back(nPandoraShwV);
                bestMatchPfoPandoraTrkNHitsW.push_back(nPandoraTrkW);
                bestMatchPfoPandoraShwNHitsW.push_back(nPandoraShwW);
                bestMatchPfoNSharedHitsTotal.push_back(sharedHitList.size());
                bestMatchPfoNSharedHitsU.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                bestMatchPfoNSharedHitsV.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                bestMatchPfoNSharedHitsW.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));
#ifdef MONITORING
                try
                {
                    const Vertex *const pRecoVertex(LArPfoHelper::GetVertex(isRecoNeutrinoFinalState ? LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first) : pfoToSharedHits.first));
                    recoVertexX = pRecoVertex->GetPosition().GetX();
                    recoVertexY = pRecoVertex->GetPosition().GetY();
                    recoVertexZ = pRecoVertex->GetPosition().GetZ();
                }
                catch (const StatusCodeException &) {}
#endif
            }

            if (isGoodMatch) ++nPrimaryMatches;

            if (isRecoNeutrinoFinalState)
            {
                const Pfo *const pRecoNeutrino(LArPfoHelper::GetParentNeutrino(pfoToSharedHits.first));
                const bool isSplitRecoNeutrino(!recoNeutrinos.empty() && !recoNeutrinos.count(pRecoNeutrino));
                if (!isSplitRecoNeutrino && isGoodMatch) ++nPrimaryGoodNuMatches;
                if (isSplitRecoNeutrino && isBeamNeutrinoFinalState && isGoodMatch) ++nPrimaryNuSplits;
                recoNeutrinos.insert(pRecoNeutrino);
            }

            if (!m_testBeamMode)
            {
                if (isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryNuMatches;
                if (!isRecoNeutrinoFinalState && isGoodMatch) ++nPrimaryCRMatches;
            }
            else
            {
                bool isTestBeam(LArPfoHelper::IsTestBeam(pfoToSharedHits.first));
                if (isTestBeam && isGoodMatch) ++nPrimaryNuMatches;
                if (!isTestBeam && isGoodMatch) ++nPrimaryCRMatches;
            }

            targetSS << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                     << "MatchedPfoId " << pfoId
                     << ", Nu " << isRecoNeutrinoFinalState;
            if (isRecoNeutrinoFinalState) targetSS << " [NuId: " << recoNuId << "]";
            targetSS << ", TB " << isRecoTestBeam
                     << ", CR " << (!isRecoNeutrinoFinalState && !isRecoTestBeam)
                     << ", PDG " << pfoToSharedHits.first->GetParticleId()
                     << ", nMatchedHits " << sharedHitList.size()
                     << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                     << ", nPfoHits " << pfoHitList.size()
                     << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                     << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;
        }

        if (mcToPfoHitSharingMap.at(pMCPrimary).empty())
        {
            targetSS << "-No matched Pfo" << std::endl;
            bestMatchPfoId.push_back(-1); bestMatchPfoPdg.push_back(0); bestMatchPfoIsRecoNu.push_back(0); bestMatchPfoRecoNuId.push_back(-1); bestMatchPfoIsTestBeam.push_back(0);
            bestMatchPfoNHitsTotal.push_back(0); bestMatchPfoNHitsU.push_back(0); bestMatchPfoNHitsV.push_back(0); bestMatchPfoNHitsW.push_back(0);
            bestMatchPfoTrueTrkNHitsU.push_back(0);
            bestMatchPfoTrueShwNHitsU.push_back(0);
            bestMatchPfoTrueTrkNHitsV.push_back(0);
            bestMatchPfoTrueShwNHitsV.push_back(0);
            bestMatchPfoTrueTrkNHitsW.push_back(0);
            bestMatchPfoTrueShwNHitsW.push_back(0);
            bestMatchPfoCNNTrkNHitsU.push_back(0);
            bestMatchPfoCNNShwNHitsU.push_back(0);
            bestMatchPfoCNNTrkNHitsV.push_back(0);
            bestMatchPfoCNNShwNHitsV.push_back(0);
            bestMatchPfoCNNTrkNHitsW.push_back(0);
            bestMatchPfoCNNShwNHitsW.push_back(0);
            bestMatchPfoPandoraTrkNHitsU.push_back(0);
            bestMatchPfoPandoraShwNHitsU.push_back(0);
            bestMatchPfoPandoraTrkNHitsV.push_back(0);
            bestMatchPfoPandoraShwNHitsV.push_back(0);
            bestMatchPfoPandoraTrkNHitsW.push_back(0);
            bestMatchPfoPandoraShwNHitsW.push_back(0);
            bestMatchPfoNSharedHitsTotal.push_back(0); bestMatchPfoNSharedHitsU.push_back(0); bestMatchPfoNSharedHitsV.push_back(0); bestMatchPfoNSharedHitsW.push_back(0);
        }

        nPrimaryMatchedPfos.push_back(nPrimaryMatches);
        nPrimaryMatchedNuPfos.push_back(nPrimaryNuMatches);
        nPrimaryMatchedCRPfos.push_back(nPrimaryCRMatches);
        nTargetMatches += nPrimaryMatches;
        nTargetNuMatches += nPrimaryNuMatches;
        nTargetCRMatches += nPrimaryCRMatches;
        nTargetGoodNuMatches += nPrimaryGoodNuMatches;
        nTargetNuSplits += nPrimaryNuSplits;
        if (0 == nPrimaryMatches) ++nTargetNuLosses;

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcNuanceCode", mcNuanceCode));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isNeutrino", isBeamNeutrinoFinalState));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isBeamParticle", isBeamParticle));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCosmicRay", isCosmicRay));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetPrimaries", nTargetPrimaries));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexX", targetVertexX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexY", targetVertexY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "targetVertexZ", targetVertexZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexX", recoVertexX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexY", recoVertexY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "recoVertexZ", recoVertexZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryId", &mcPrimaryId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPdg", &mcPrimaryPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryE", &mcPrimaryE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPX", &mcPrimaryPX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPY", &mcPrimaryPY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryPZ", &mcPrimaryPZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxX", &mcPrimaryVtxX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxY", &mcPrimaryVtxY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryVtxZ", &mcPrimaryVtxZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndX", &mcPrimaryEndX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndY", &mcPrimaryEndY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryEndZ", &mcPrimaryEndZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsTotal", &nMCHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsU", &nMCHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNTrkHitsU", &nMCTrkHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNShwHitsU", &nMCShwHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsV", &nMCHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNTrkHitsV", &nMCTrkHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNShwHitsV", &nMCShwHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNHitsW", &nMCHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNTrkHitsW", &nMCTrkHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPrimaryNShwHitsW", &nMCShwHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedPfos", &nPrimaryMatchedPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedNuPfos", &nPrimaryMatchedNuPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPrimaryMatchedCRPfos", &nPrimaryMatchedCRPfos));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoId", &bestMatchPfoId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPdg", &bestMatchPfoPdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsRecoNu", &bestMatchPfoIsRecoNu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoRecoNuId", &bestMatchPfoRecoNuId));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsTotal", &bestMatchPfoNHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsU", &bestMatchPfoNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTrueTrkNHitsU", &bestMatchPfoTrueTrkNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTrueShwNHitsU", &bestMatchPfoTrueShwNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoCNNTrkNHitsU", &bestMatchPfoCNNTrkNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoCNNShwNHitsU", &bestMatchPfoCNNShwNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPandoraTrkNHitsU", &bestMatchPfoPandoraTrkNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPandoraShwNHitsU", &bestMatchPfoPandoraShwNHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsV", &bestMatchPfoNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTrueTrkNHitsV", &bestMatchPfoTrueTrkNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTrueShwNHitsV", &bestMatchPfoTrueShwNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoCNNTrkNHitsV", &bestMatchPfoCNNTrkNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoCNNShwNHitsV", &bestMatchPfoCNNShwNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPandoraTrkNHitsV", &bestMatchPfoPandoraTrkNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPandoraShwNHitsV", &bestMatchPfoPandoraShwNHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNHitsW", &bestMatchPfoNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTrueTrkNHitsW", &bestMatchPfoTrueTrkNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoTrueShwNHitsW", &bestMatchPfoTrueShwNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoCNNTrkNHitsW", &bestMatchPfoCNNTrkNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoCNNShwNHitsW", &bestMatchPfoCNNShwNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPandoraTrkNHitsW", &bestMatchPfoPandoraTrkNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoPandoraShwNHitsW", &bestMatchPfoPandoraShwNHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsTotal", &bestMatchPfoNSharedHitsTotal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsU", &bestMatchPfoNSharedHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsV", &bestMatchPfoNSharedHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoNSharedHitsW", &bestMatchPfoNSharedHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetMatches", nTargetMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuMatches", nTargetNuMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetCRMatches", nTargetCRMatches));

            if (m_testBeamMode)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchPfoIsTestBeam", &bestMatchPfoIsTestBeam));
            }

            if (!m_testBeamMode)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetGoodNuMatches", nTargetGoodNuMatches));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuSplits", nTargetNuSplits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTargetNuLosses", nTargetNuLosses));
            }
        }

        if (isLastNeutrinoPrimary || isBeamParticle || isCosmicRay)
        {
            const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(associatedMCPrimaries));
#ifdef MONITORING
            const int interactionTypeInt(static_cast<int>(interactionType));
#endif
            // ATTN Some redundancy introduced to contributing variables
            const int isCorrectNu(isBeamNeutrinoFinalState && (nTargetGoodNuMatches == nTargetNuMatches) && (nTargetGoodNuMatches == nTargetPrimaries) && (nTargetCRMatches == 0) && (nTargetNuSplits == 0) && (nTargetNuLosses == 0));
            const int isCorrectTB(isBeamParticle && (nTargetNuMatches == 1) && (nTargetCRMatches == 0));
            const int isCorrectCR(isCosmicRay && (nTargetNuMatches == 0) && (nTargetCRMatches == 1));
            const int isFakeNu(isCosmicRay && (nTargetNuMatches > 0));
            const int isFakeCR(!isCosmicRay && (nTargetCRMatches > 0));
            const int isSplitNu(!isCosmicRay && ((nTargetNuMatches > nTargetPrimaries) || (nTargetNuSplits > 0)));
            const int isSplitCR(isCosmicRay && (nTargetCRMatches > 1));
            const int isLost(nTargetMatches == 0);

            std::stringstream outcomeSS;
            outcomeSS << LArInteractionTypeHelper::ToString(interactionType) << " (Nuance " << mcNuanceCode << ", Nu " << isBeamNeutrinoFinalState << ", TB " << isBeamParticle << ", CR " << isCosmicRay << ")" << std::endl;

            if (isLastNeutrinoPrimary) ++nTotalNu;
            if (isBeamParticle) ++nTotalTB;
            if (isCosmicRay) ++nTotalCR;
            if (isCorrectNu) ++nCorrectNu;
            if (isCorrectTB) ++nCorrectTB;
            if (isCorrectCR) ++nCorrectCR;
            if (isFakeNu) ++nFakeNu;
            if (isFakeCR) ++nFakeCR;
            if (isSplitNu) ++nSplitNu;
            if (isSplitCR) ++nSplitCR;
            if (isLost) ++nLost;

            if (isCorrectNu) outcomeSS << "IsCorrectNu ";
            if (isCorrectTB) outcomeSS << "IsCorrectTB ";
            if (isCorrectCR) outcomeSS << "IsCorrectCR ";
            if (isFakeNu) outcomeSS << "IsFakeNu ";
            if (isFakeCR) outcomeSS << "IsFakeCR ";
            if (isSplitNu) outcomeSS << "isSplitNu ";
            if (isSplitCR) outcomeSS << "IsSplitCR ";
            if (isLost) outcomeSS << "IsLost ";
            if (nTargetNuMatches > 0) outcomeSS << "(NNuMatches: " << nTargetNuMatches << ") ";
            if (nTargetNuLosses > 0) outcomeSS << "(NNuLosses: " << nTargetNuLosses << ") ";
            if (nTargetNuSplits > 0) outcomeSS << "(NNuSplits: " << nTargetNuSplits << ") ";
            if (nTargetCRMatches > 0) outcomeSS << "(NCRMatches: " << nTargetCRMatches << ") ";
            if (printToScreen) std::cout << outcomeSS.str() << std::endl << targetSS.str() << std::endl;

            if (fillTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "interactionType", interactionTypeInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectNu", isCorrectNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectTB", isCorrectTB));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectCR", isCorrectCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeNu", isFakeNu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isFakeCR", isFakeCR));
                if (!m_testBeamMode)
                {
                    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitNu", isSplitNu));
                }
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isSplitCR", isSplitCR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isLost", isLost));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            }

            targetSS.str(std::string()); targetSS.clear();
            recoNeutrinos.clear(); associatedMCPrimaries.clear();
            nTargetMatches = 0; nTargetNuMatches = 0; nTargetCRMatches = 0; nTargetGoodNuMatches = 0; nTargetNuSplits = 0; nTargetNuLosses = 0;
            mcPrimaryId.clear(); mcPrimaryPdg.clear(); nMCHitsTotal.clear(); nMCHitsU.clear(); nMCHitsV.clear(); nMCHitsW.clear(); nMCTrkHitsU.clear(); nMCTrkHitsV.clear(); nMCTrkHitsW.clear(); nMCShwHitsU.clear(); nMCShwHitsV.clear(); nMCShwHitsW.clear();
            mcPrimaryE.clear(); mcPrimaryPX.clear(); mcPrimaryPY.clear(); mcPrimaryPZ.clear();
            mcPrimaryVtxX.clear(); mcPrimaryVtxY.clear(); mcPrimaryVtxZ.clear(); mcPrimaryEndX.clear(); mcPrimaryEndY.clear(); mcPrimaryEndZ.clear();
            nPrimaryMatchedPfos.clear(); nPrimaryMatchedNuPfos.clear(); nPrimaryMatchedCRPfos.clear();
            bestMatchPfoId.clear(); bestMatchPfoPdg.clear(); bestMatchPfoIsRecoNu.clear(); bestMatchPfoRecoNuId.clear(); bestMatchPfoIsTestBeam.clear();
            bestMatchPfoNHitsTotal.clear(); bestMatchPfoNHitsU.clear(); bestMatchPfoNHitsV.clear(); bestMatchPfoNHitsW.clear();
            bestMatchPfoTrueTrkNHitsU.clear();
            bestMatchPfoTrueShwNHitsU.clear();
            bestMatchPfoTrueTrkNHitsV.clear();
            bestMatchPfoTrueShwNHitsV.clear();
            bestMatchPfoTrueTrkNHitsW.clear();
            bestMatchPfoTrueShwNHitsW.clear();
            bestMatchPfoCNNTrkNHitsU.clear();
            bestMatchPfoCNNShwNHitsU.clear();
            bestMatchPfoPandoraTrkNHitsU.clear();
            bestMatchPfoPandoraShwNHitsU.clear();
            bestMatchPfoCNNTrkNHitsV.clear();
            bestMatchPfoCNNShwNHitsV.clear();
            bestMatchPfoPandoraTrkNHitsV.clear();
            bestMatchPfoPandoraShwNHitsV.clear();
            bestMatchPfoCNNTrkNHitsW.clear();
            bestMatchPfoCNNShwNHitsW.clear();
            bestMatchPfoPandoraTrkNHitsW.clear();
            bestMatchPfoPandoraShwNHitsW.clear();
            bestMatchPfoNSharedHitsTotal.clear(); bestMatchPfoNSharedHitsU.clear(); bestMatchPfoNSharedHitsV.clear(); bestMatchPfoNSharedHitsW.clear();
        }
    }

    if (useInterpretedMatching)
    {
        std::stringstream summarySS;
        summarySS << "---SUMMARY--------------------------------------------------------------------------------------" << std::endl;
        if (nTotalNu > 0) summarySS << "#CorrectNu: " << nCorrectNu << "/" << nTotalNu << ", Fraction: " << (nTotalNu > 0 ? static_cast<float>(nCorrectNu) / static_cast<float>(nTotalNu) : 0.f) << std::endl;
        if (nTotalTB > 0) summarySS << "#CorrectTB: " << nCorrectTB << "/" << nTotalTB << ", Fraction: " << (nTotalTB > 0 ? static_cast<float>(nCorrectTB) / static_cast<float>(nTotalTB) : 0.f) << std::endl;
        if (nTotalCR > 0) summarySS << "#CorrectCR: " << nCorrectCR << "/" << nTotalCR << ", Fraction: " << (nTotalCR > 0 ? static_cast<float>(nCorrectCR) / static_cast<float>(nTotalCR) : 0.f) << std::endl;
        if (nFakeNu > 0) summarySS << "#FakeNu: " << nFakeNu << " ";
        if (nFakeCR > 0) summarySS << "#FakeCR: " << nFakeCR << " ";
        if (nSplitNu > 0) summarySS << "#SplitNu: " << nSplitNu << " ";
        if (nSplitCR > 0) summarySS << "#SplitCR: " << nSplitCR << " ";
        if (nLost > 0) summarySS << "#Lost: " << nLost << " ";
        if (nFakeNu || nFakeCR || nSplitNu || nSplitCR || nLost) summarySS << std::endl;
        if (printToScreen) std::cout << summarySS.str();
    }

    if (printToScreen) std::cout << "------------------------------------------------------------------------------------------------" << std::endl << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::CountTrkShwHitsByType(const HitType hitType, const CaloHitList &caloHitList, int &nTrack, int &nShower) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (hitType == pCaloHit->GetHitType())
        {
            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                if (pMCParticle->GetParticleId() == 22 || std::abs(pMCParticle->GetParticleId()) == 11)
                {
                    nShower++;
                }
                else
                {
                    nTrack++;
                }
            }
            catch (...) {}
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::CountCNNTrkShwHitsByType(const HitType hitType, const CaloHitList &targetCaloHitList, const CaloHitList &allCaloHitList, int &nTrack, int &nShower) const
{
    for (const CaloHit *pTargetCaloHit : targetCaloHitList)
    {
        if (hitType != pTargetCaloHit->GetHitType())
            continue;

        TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);

        for (const CaloHit *pNeighbourCaloHit : allCaloHitList)
        {
            if (hitType != pNeighbourCaloHit->GetHitType())
                continue;

            CartesianVector relativePosition(pNeighbourCaloHit->GetPositionVector() - pTargetCaloHit->GetPositionVector());
            twoDHistogram.Fill(relativePosition.GetX(), relativePosition.GetZ(), pNeighbourCaloHit->GetInputEnergy());
        }

        KerasModel::DataBlock2D dataBlock2D;
        this->HistogramToDataBlock(twoDHistogram, dataBlock2D);
        Data1D outputData1D;
        m_kerasModel.CalculateOutput(&dataBlock2D, outputData1D, this);

        if (m_verbose)
        {
            for (unsigned int counter = 0; counter < outputData1D.GetSizeI(); counter++)
                std::cout << "Class " << counter << ", outcome " << outputData1D.Get(counter) << std::endl;
        }

        if (outputData1D.Get(0) > outputData1D.Get(1))
        {
            nShower++;
        }
        else
        {
            nTrack++;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::HistogramToDataBlock(const TwoDHistogram &twoDHistogram, KerasModel::DataBlock2D &dataBlock2D) const
{
    Data3D data3D;
    Data2D data2D;
    for (int yBin = 0; yBin < twoDHistogram.GetNBinsY(); yBin++)
    {
        Data1D data1D;
        for (int xBin = 0; xBin < twoDHistogram.GetNBinsX(); xBin++)
        {
            data1D.Append(twoDHistogram.GetBinContent(xBin, yBin) * 256.f / 10000.f );
        }
        data2D.Append(data1D);
    }
    data3D.Append(data2D);
    dataBlock2D.SetData(data3D);
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::CountPandoraTrkShwHitsByType(const HitType hitType, const CaloHitList &caloHitList, int &nTrack, int &nShower) const
{
    for (const CaloHit *pCaloHit : caloHitList)
    {
        if (hitType != pCaloHit->GetHitType())
            continue;

        if (m_hitToIntMap.find(pCaloHit) != m_hitToIntMap.end())
        {
            if (m_hitToIntMap.at(pCaloHit) == 0)
            {
                nShower++;
            }
            else if (m_hitToIntMap.at(pCaloHit) == 1)
            {
                nTrack++;
            }
        }
        else
        {
            std::cout << "Hit missing from map" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::InterpretMatching(const ValidationInfo &validationInfo, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetAllMCParticleToHitsMap()}, mcPrimaryVector);

    PfoSet usedPfos;
    while (this->GetStrongestPfoMatch(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap)) {}
    this->GetRemainingPfoMatches(validationInfo, mcPrimaryVector, usedPfos, interpretedMCToPfoHitSharingMap);

    // Ensure all primaries have an entry, and sorting is as desired
    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        LArMCParticleHelper::PfoToSharedHitsVector &pfoHitPairs(interpretedMCToPfoHitSharingMap[pMCPrimary]);
        std::sort(pfoHitPairs.begin(), pfoHitPairs.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool {
            return ((a.second.size() != b.second.size()) ? a.second.size() > b.second.size() : LArPfoHelper::SortByNHits(a.first, b.first)); });
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::GetStrongestPfoMatch(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    const MCParticle *pBestMCParticle(nullptr);
    LArMCParticleHelper::PfoCaloHitListPair bestPfoHitPair(nullptr, CaloHitList());

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (interpretedMCToPfoHitSharingMap.count(pMCPrimary))
            continue;

        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            if (!this->IsGoodMatch(validationInfo.GetAllMCParticleToHitsMap().at(pMCPrimary), validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first), pfoToSharedHits.second))
                continue;

            if (pfoToSharedHits.second.size() > bestPfoHitPair.second.size())
            {
                pBestMCParticle = pMCPrimary;
                bestPfoHitPair = pfoToSharedHits;
            }
        }
    }

    if (!pBestMCParticle || !bestPfoHitPair.first)
        return false;

    interpretedMCToPfoHitSharingMap[pBestMCParticle].push_back(bestPfoHitPair);
    usedPfos.insert(bestPfoHitPair.first);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventValidationAlgorithm::GetRemainingPfoMatches(const ValidationInfo &validationInfo, const MCParticleVector &mcPrimaryVector,
    const PfoSet &usedPfos, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const
{
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;

    for (const MCParticle *const pMCPrimary : mcPrimaryVector)
    {
        if (!m_useSmallPrimaries && !validationInfo.GetTargetMCParticleToHitsMap().count(pMCPrimary))
            continue;

        if (!validationInfo.GetMCToPfoHitSharingMap().count(pMCPrimary))
            continue;

        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : validationInfo.GetMCToPfoHitSharingMap().at(pMCPrimary))
        {
            if (usedPfos.count(pfoToSharedHits.first))
                continue;

            const LArMCParticleHelper::MCParticleCaloHitListPair mcParticleToHits(pMCPrimary, pfoToSharedHits.second);
            LArMCParticleHelper::PfoToMCParticleHitSharingMap::iterator iter(pfoToMCParticleHitSharingMap.find(pfoToSharedHits.first));

            if (pfoToMCParticleHitSharingMap.end() == iter)
            {
                pfoToMCParticleHitSharingMap[pfoToSharedHits.first].push_back(mcParticleToHits);
            }
            else
            {
                if (1 != iter->second.size())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                LArMCParticleHelper::MCParticleCaloHitListPair &originalMCParticleToHits(iter->second.at(0));

                if (mcParticleToHits.second.size() > originalMCParticleToHits.second.size())
                    originalMCParticleToHits = mcParticleToHits;
            }
        }
    }

    for (const auto &mapEntry : pfoToMCParticleHitSharingMap)
    {
        const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleToHits(mapEntry.second.at(0));
        interpretedMCToPfoHitSharingMap[mcParticleToHits.first].push_back(LArMCParticleHelper::PfoCaloHitListPair(mapEntry.first, mcParticleToHits.second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventValidationAlgorithm::IsGoodMatch(const CaloHitList &trueHits, const CaloHitList &recoHits, const CaloHitList &sharedHits) const
{
    const float purity((recoHits.size() > 0) ? static_cast<float>(sharedHits.size()) / static_cast<float>(recoHits.size()) : 0.f);
    const float completeness((trueHits.size() > 0) ? static_cast<float>(sharedHits.size()) / static_cast<float>(trueHits.size()) : 0.f);

    return ((sharedHits.size() >= m_matchingMinSharedHits) && (purity >= m_matchingMinPurity) && (completeness >= m_matchingMinCompleteness));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridSize", m_gridSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridDimensions", m_gridDimensions));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbose", m_verbose));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CNNModelName", m_cnnModelName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CNNModelXml", m_cnnModelXml));

    m_kerasModel.Initialize(m_cnnModelXml, m_cnnModelName);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrueNeutrinosOnly", m_useTrueNeutrinosOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TestBeamMode", m_testBeamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintAllToScreen", m_printAllToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintMatchingToScreen", m_printMatchingToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseSmallPrimaries", m_useSmallPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinSharedHits", m_matchingMinSharedHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinCompleteness", m_matchingMinCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MatchingMinPurity", m_matchingMinPurity));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "FileIdentifier", m_fileIdentifier));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
