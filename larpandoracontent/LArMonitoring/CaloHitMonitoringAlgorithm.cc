/**
 *  @file   larpandoracontent/LArMonitoring/CaloHitMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the calo hit monitoring algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArDeepLearningHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/CaloHitMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CaloHitMonitoringAlgorithm::CaloHitMonitoringAlgorithm() :
    m_treeName("DeepLearningMonitoringTree"),
    m_fileName("DeepLearningMonitoring.root"),
    m_gridSize(16),
    m_gridDimensions(50),
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitMonitoringAlgorithm::~CaloHitMonitoringAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitMonitoringAlgorithm::Run()
{
    m_eventNumber++;

    const PfoList *pPfoList = nullptr;
    (void) PandoraContentApi::GetList(*this, m_pfoListName, pPfoList);

    if (!pPfoList)
        return STATUS_CODE_SUCCESS;

    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

    for (const Pfo *pPfo : allConnectedPfos)
    {
        CaloHitList caloHitList;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);
        LArPfoHelper::GetIsolatedCaloHits(pPfo, TPC_VIEW_W, caloHitList);

        int isShowerPandora(0);
        int isTestBeamPfo(LArPfoHelper::IsTestBeam(pPfo) ? 1 : 0);

        if (pPfo->GetParticleId() == 22 || std::abs(pPfo->GetParticleId()) == 11)
            isShowerPandora = 1;

        for (const CaloHit *pCaloHit : caloHitList)
        {
            int isShowerTruth(0);
            int isTestBeamTruth(0);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

                if (pMCParticle->GetParticleId() == 22 || std::abs(pMCParticle->GetParticleId()) == 11)
                    isShowerTruth = 1;

                if (LArMCParticleHelper::IsBeamParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle)))
                    isTestBeamTruth = 1;
            }
            catch(...)
            {
                continue;
            }

            TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);

            for (const CaloHit *pNeighbourCaloHit : caloHitList)
            {
                CartesianVector relativePosition(pNeighbourCaloHit->GetPositionVector() - pCaloHit->GetPositionVector());
                twoDHistogram.Fill(relativePosition.GetX(), relativePosition.GetZ(), pNeighbourCaloHit->GetInputEnergy());
            }

            KerasModel::DataBlock2D dataBlock2D;
            LArDeepLearningHelper::HistogramToDataBlock(twoDHistogram, dataBlock2D);
            Data1D outputData1D;
            m_kerasModel.CalculateOutput(&dataBlock2D, outputData1D, this);

            int isShowerCNN(0);

            if (outputData1D.Get(0) > outputData1D.Get(1))
                isShowerCNN = 1;

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_eventNumber));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTestBeamPfo", isTestBeamPfo));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTestBeamTruth", isTestBeamTruth));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isShowerTruth", isShowerTruth));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isShowerPandora", isShowerPandora));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isShowerCNN", isShowerCNN));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CNNModelName", m_cnnModelName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CNNModelXml", m_cnnModelXml));

    m_kerasModel.Initialize(m_cnnModelXml, m_cnnModelName);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridSize", m_gridSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridDimensions", m_gridDimensions));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
