/**
 *  @file   larpandoracontent/LArTrackShowerId/CutPfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/CutPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CutPfoCharacterisationAlgorithm::CutPfoCharacterisationAlgorithm() :
    m_postBranchAddition(false),
    m_slidingFitWindow(5),
    m_slidingShowerFitWindow(10),
    m_maxShowerLengthCut(80.f),
    m_dTdLWidthRatioCut(0.045f),
    m_vertexDistanceRatioCut(0.6f),
    m_showerWidthRatioCut(0.2f),
    m_writeToTree(false),
    m_visualizeTruth(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CutPfoCharacterisationAlgorithm::~CutPfoCharacterisationAlgorithm()
{
    if (m_writeToTree)
         PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    float straightLineLength(-1.f);
    float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude();

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            dTdLMin = std::min(dTdLMin, static_cast<float>(mapEntry.second.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(mapEntry.second.GetGradient()));
        }
    }
    catch (const StatusCodeException &)
    {
    }

    float dTdLWidthRatio((dTdLMax - dTdLMin) / straightLineLength);
    bool straightLineLenghtOk(true);

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        straightLineLenghtOk = false;

    bool straightLineLengthCutPass(false);

    if (straightLineLength > m_maxShowerLengthCut)
        straightLineLengthCutPass = true;

    bool dTdLWidthRatioCutPass(false);
    if (straightLineLenghtOk)
    {
        if ((dTdLMax - dTdLMin) / straightLineLength > m_dTdLWidthRatioCut)
            dTdLWidthRatioCutPass = true;
    }

    const float vertexDistance(CutClusterCharacterisationAlgorithm::GetVertexDistance(this, pCluster));
    bool vertexDistanceRatioCutPass(false);

    if (straightLineLenghtOk)
    {
        if ((vertexDistance > std::numeric_limits<float>::epsilon()) && ((vertexDistance / straightLineLength) > m_vertexDistanceRatioCut))
            vertexDistanceRatioCutPass = true;
    }

    const float showerFitWidth(CutClusterCharacterisationAlgorithm::GetShowerFitWidth(this, pCluster, m_slidingShowerFitWindow));
    const float showerWidthRatio(showerFitWidth / straightLineLength);
    bool showerWidthRatioCutPass(false);

    if (straightLineLenghtOk)
    {
        if ((showerFitWidth < std::numeric_limits<float>::epsilon()) || (showerWidthRatio > m_showerWidthRatioCut))
            showerWidthRatioCutPass = true;
    }

    int isRecoTrack(1);

    if (!straightLineLenghtOk || dTdLWidthRatioCutPass || vertexDistanceRatioCutPass || showerWidthRatioCutPass)
    {
        isRecoTrack = 0;
    }

    std::cout << "straightLineLength   : " << straightLineLength << std::endl;
    std::cout << "dTdLWidthRatio       : " << (straightLineLenghtOk ? dTdLWidthRatio : -1.f) << std::endl;
    std::cout << "vertexDistanceRatio  : " << (straightLineLenghtOk ? (vertexDistance / straightLineLength) : -1.f) << std::endl;
    std::cout << "showerWidthRatio     : " << (straightLineLenghtOk ? showerWidthRatio : -1.f) << std::endl;
    std::cout << "Cuts-Shower          : LT " << m_maxShowerLengthCut << ", GT " << m_dTdLWidthRatioCut << ", GT " << m_vertexDistanceRatioCut << ", GT " << m_showerWidthRatioCut << ". IsRecoTrack " << isRecoTrack << std::endl;

    int nTrackHits(0), nShowerHits(0), nBadHits(0);
    const int isTrack(this->GetTruthIsTrack(pCluster, nTrackHits, nShowerHits, nBadHits));

    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTrack", isTrack));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isRecoTrack", isRecoTrack));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrackHits", nTrackHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nShowerHits", nShowerHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nBadHits", nBadHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength", straightLineLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLMax", dTdLMax));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLMin", dTdLMin));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLWidthRatio", dTdLWidthRatio));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance", vertexDistance));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth", showerFitWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerWidthRatio", showerWidthRatio));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }

    if (!straightLineLenghtOk || dTdLWidthRatioCutPass || vertexDistanceRatioCutPass || showerWidthRatioCutPass)
    {
        return false;
    }
    else if (straightLineLengthCutPass)
    {
        return true;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutPfoCharacterisationAlgorithm::GetTruthIsTrack(const Cluster *const pCluster, int &nTrackHits, int &nShowerHits, int &nBadHits) const
{
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    nTrackHits = 0;
    nShowerHits = 0;
    nBadHits = 0;

    CaloHitList caloHitListTrack;
    CaloHitList caloHitListShower;
    CaloHitList caloHitListBad;

    HitType hitType(TPC_VIEW_U);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        hitType = pCaloHit->GetHitType();
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const int mcParticleId(pHitParticle->GetParticleId());

            if (std::fabs(mcParticleId) == 11 || mcParticleId == 22)
            {
                nShowerHits++;
                caloHitListShower.push_back(pCaloHit);
            }
            else
            {
                nTrackHits++;
                caloHitListTrack.push_back(pCaloHit);
            }
        }
        catch (StatusCodeException &)
        {
            nBadHits++;
            caloHitListBad.push_back(pCaloHit);
        }
    }

    if (m_visualizeTruth)
    {
        std::cout << "nShowerHits : " << nShowerHits << std::endl;
        std::cout << "nTrackHits  : " << nTrackHits << std::endl;
        std::cout << "nBadHits    : " << nBadHits << std::endl;

        ClusterList clusterList;
        clusterList.push_back(pCluster);

        const CaloHitList *pCaloHitList(nullptr);
        if (hitType == TPC_VIEW_U)
        {
            PandoraContentApi::GetList(*this, "CaloHitListU", pCaloHitList);
        }
        else if (hitType == TPC_VIEW_V)
        {
            PandoraContentApi::GetList(*this, "CaloHitListV", pCaloHitList);
        }
        else if (hitType == TPC_VIEW_W)
        {
            PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitList);
        }

        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitList, "AllHits", BLACK));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "ACluster", RED, true));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListShower, "ShowerLikeHits", BLUE));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListTrack, "TrackLikeHits", RED));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListBad, "BadHits", GRAY));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    if (nShowerHits > nTrackHits)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const /*pPfo*/) const
{
	throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CutPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PostBranchAddition", m_postBranchAddition));

    // Allow change in default values via a single xml tag, can subsequently override all individual values below, if required
    if (m_postBranchAddition)
    {
        m_maxShowerLengthCut = 80.f;
        m_dTdLWidthRatioCut = 0.03f;
        m_vertexDistanceRatioCut = 1.f;
        m_showerWidthRatioCut = 0.3f;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLengthCut", m_maxShowerLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DTDLWidthRatioCut", m_dTdLWidthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDistanceRatioCut", m_vertexDistanceRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerWidthRatioCut", m_showerWidthRatioCut));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "WriteToTree", m_writeToTree));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeTruth", m_visualizeTruth));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
