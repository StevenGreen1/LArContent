/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/DeepLearningSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the layer splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArDeepLearningHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/DeepLearningSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeepLearningSplittingAlgorithm::DeepLearningSplittingAlgorithm() :
    m_gridSize(16),
    m_gridDimensions(50)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    // ATTN: Assumes correct list is passed
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    CaloHitList caloHitsInCluster;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitsInCluster);
    caloHitsInCluster.insert(caloHitsInCluster.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

    for (const CaloHit *pCaloHit : caloHitsInCluster)
    {
        TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);

        for (const CaloHit *pNeighbourCaloHit : *pCaloHitList)
        {
            CartesianVector relativePosition(pNeighbourCaloHit->GetPositionVector() - pCaloHit->GetPositionVector());
            twoDHistogram.Fill(relativePosition.GetX(), relativePosition.GetZ(), pNeighbourCaloHit->GetInputEnergy());
        }

        KerasModel::DataBlock2D dataBlock2D;
        LArDeepLearningHelper::HistogramToDataBlock(twoDHistogram, dataBlock2D);
        Data1D outputData1D;
        m_kerasModel.CalculateOutput(&dataBlock2D, outputData1D, this);

        if (outputData1D.Get(0) > outputData1D.Get(1))
        {
            firstHitList.push_back(pCaloHit);
        }
        else
        {
            secondHitList.push_back(pCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeepLearningSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridSize", m_gridSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridDimensions", m_gridDimensions));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CNNModelName", m_cnnModelName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CNNModelXml", m_cnnModelXml));

    m_kerasModel.Initialize(m_cnnModelXml, m_cnnModelName);

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
