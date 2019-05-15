/**
 *  @file   larpandoracontent/LArVertex/BdtVertexSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the Bdt vertex selection algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArVertex/BdtVertexSelectionAlgorithm.h"

#include <random>

using namespace pandora;

namespace lar_content
{

BdtVertexSelectionAlgorithm::BdtVertexSelectionAlgorithm() :
    MLVertexSelectionBaseAlgorithm(),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Vertex * BdtVertexSelectionAlgorithm::CompareVertices(const VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const LArMvaHelper::MvaFeatureVector &eventFeatureList, const bool useRPhi, const bool isRegion) const
{
    AdaBoostDecisionTree adaBoostDecisionTree(isRegion ? m_bdtRegion : m_bdtVertex);

    const Vertex *pBestVertex(vertexVector.front());
    LArMvaHelper::MvaFeatureVector chosenFeatureList;

    VertexFeatureInfo chosenVertexFeatureInfo(vertexFeatureInfoMap.at(pBestVertex));
    this->AddVertexFeaturesToVector(chosenVertexFeatureInfo, chosenFeatureList, useRPhi);

    for (const Vertex *const pVertex : vertexVector)
    {
        if (pVertex == pBestVertex)
            continue;

        LArMvaHelper::MvaFeatureVector featureList;
        VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
        this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList, useRPhi);

        if (LArMvaHelper::Classify(adaBoostDecisionTree, eventFeatureList, featureList, chosenFeatureList))
        {
            pBestVertex = pVertex;
            chosenFeatureList = featureList;
        }
    }

    return pBestVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BdtVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // ATTN : Need to know if in training mode before reading machine learning specific parameters
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MLVertexSelectionBaseAlgorithm::ReadSettings(xmlHandle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BdtFileName", m_bdtFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RegionBdtName", m_regionBdtName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexBdtName", m_vertexBdtName));

    if ((!m_trainingSetMode || m_allowClassifyDuringTraining))
    {
        if (m_bdtFileName.empty() || m_regionBdtName.empty() || m_vertexBdtName.empty())
        {
            std::cout << "BdtVertexSelectionAlgorithm: BdtFileName, RegionBdtName and VertexBdtName must be set if training set mode is" <<
                         "off or we allow classification during training" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullBdtFileName(LArFileHelper::FindFileInPath(m_bdtFileName, m_filePathEnvironmentVariable));
        m_bdtRegion.Initialize(fullBdtFileName, m_regionBdtName);
        m_bdtVertex.Initialize(fullBdtFileName, m_vertexBdtName);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
