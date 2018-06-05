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

    int isTrack(this->GetTruthIsTrack(pCluster));

    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTrack", isTrack));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength", straightLineLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLMax", dTdLMax));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLMin", dTdLMin));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "dTdLWidthRatio", dTdLWidthRatio));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance", vertexDistance));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth", showerFitWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerWidthRatio", showerWidthRatio));
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

bool CutPfoCharacterisationAlgorithm::GetTruthIsTrack(const Cluster *const pCluster) const
{
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    int nHitsTrack(0), nHitsShower(0);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
        const int mcParticleId(pHitParticle->GetParticleId());

        if (std::fabs(mcParticleId) == 11 || mcParticleId == 22)
        {
            nHitsShower++;
        }
        else
        {
            nHitsTrack++;
        }
    }

    if (m_visualizeTruth)
    {
        std::cout << "nHitsShower : " << nHitsShower << std::endl;
        std::cout << "nHitsTrack  : " << nHitsTrack << std::endl;

        ClusterList clusterList;
        clusterList.push_back(pCluster);

        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "ACluster", RED, true));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    if (nHitsShower > nHitsTrack)
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
