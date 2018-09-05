/**
 *  @file   larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.cc
 *
 *  @brief  Implementation of the pfo validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CaloHitEventDumpAlgorithm::CaloHitEventDumpAlgorithm() :
    m_verbose(false),
    m_treeName("DeepLearningMonitoringTree"),
    m_fileName("DeepLearningMonitoring.root"),
    m_gridSize(16),
    m_gridDimensions(50),
    m_useTrainingMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitEventDumpAlgorithm::~CaloHitEventDumpAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitEventDumpAlgorithm::Run()
{
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    CaloHitList wCaloHitList;

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        if (TPC_VIEW_W == pCaloHit->GetHitType())
            wCaloHitList.push_back(pCaloHit);
    }

    if (m_useTrainingMode)
    {
        for (const CaloHit *pTargetCaloHit : wCaloHitList)
        {
            LArMvaHelper::MvaFeatureVector featureVector;
            int targetParticleId(0);

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pTargetCaloHit));
                targetParticleId = pMCParticle->GetParticleId();
            }
            catch(...)
            {
                continue;
            }

            featureVector.push_back(static_cast<double>(targetParticleId));

            TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);

            for (const CaloHit *pNeighbourCaloHit : wCaloHitList)
            {
                CartesianVector relativePosition(pNeighbourCaloHit->GetPositionVector() - pTargetCaloHit->GetPositionVector());
                twoDHistogram.Fill(relativePosition.GetX(), relativePosition.GetZ(), pNeighbourCaloHit->GetInputEnergy());
            }

            for (int xBin = 0; xBin < twoDHistogram.GetNBinsX(); xBin++)
            {
                for (int yBin = 0; yBin < twoDHistogram.GetNBinsY(); yBin++)
                {
                    featureVector.push_back(static_cast<double>(twoDHistogram.GetBinContent(xBin, yBin)));
                }
            }

            // ATTN: The hard coded in true is a redundant variable for the deep learning cases
            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, true, featureVector);
        }
        return STATUS_CODE_SUCCESS;
    }

    KerasModel kerasModel;
    kerasModel.Initialize("Testing.xml", "TrackShowerId");

// SingleImage.txt Test
//    TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);
    TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridSize/2.f,  m_gridSize/2.f, m_gridSize, -1.f * m_gridSize/2.f,  m_gridSize/2.f);
    const int xOffset(8);
    const int yOffset(8);
/*
    twoDHistogram.Fill(8-xOffset, 8-yOffset, 66 * 10000.f / 256.f);
    twoDHistogram.Fill(8-xOffset, 9-yOffset, 187 * 10000.f / 256.f);
    twoDHistogram.Fill(9-xOffset, 8-yOffset, 11 * 10000.f / 256.f);
    twoDHistogram.Fill(9-xOffset, 9-yOffset, 9 * 10000.f / 256.f);
    twoDHistogram.Fill(9-xOffset, 10-yOffset, 63 * 10000.f / 256.f);
    twoDHistogram.Fill(9-xOffset, 11-yOffset, 2 * 10000.f / 256.f);
    twoDHistogram.Fill(8-xOffset, 10-yOffset, 5 * 10000.f / 256.f);
    twoDHistogram.Fill(10-xOffset, 10-yOffset, 37 * 10000.f / 256.f);
    twoDHistogram.Fill(10-xOffset, 11-yOffset, 17 * 10000.f / 256.f);
    twoDHistogram.Fill(11-xOffset, 10-yOffset, 10 * 10000.f / 256.f);
    twoDHistogram.Fill(11-xOffset, 11-yOffset, 65 * 10000.f / 256.f);
    twoDHistogram.Fill(11-xOffset, 12-yOffset, 10 * 10000.f / 256.f);
    twoDHistogram.Fill(12-xOffset, 12-yOffset, 60 * 10000.f / 256.f);
    twoDHistogram.Fill(13-xOffset, 12-yOffset, 29 * 10000.f / 256.f);
    twoDHistogram.Fill(13-xOffset, 13-yOffset, 44 * 10000.f / 256.f);
    twoDHistogram.Fill(14-xOffset, 13-yOffset, 50 * 10000.f / 256.f);
    twoDHistogram.Fill(14-xOffset, 14-yOffset, 10 * 10000.f / 256.f);
    twoDHistogram.Fill(15-xOffset, 14-yOffset, 78 * 10000.f / 256.f);
*/

    twoDHistogram.Fill(8-xOffset, 8-yOffset, 66 * 10000.f / 256.f);
    twoDHistogram.Fill(9-xOffset, 8-yOffset, 187 * 10000.f / 256.f);
    twoDHistogram.Fill(8-xOffset, 9-yOffset, 11 * 10000.f / 256.f);
    twoDHistogram.Fill(9-xOffset, 9-yOffset, 9 * 10000.f / 256.f);
    twoDHistogram.Fill(10-xOffset, 9-yOffset, 63 * 10000.f / 256.f);
    twoDHistogram.Fill(11-xOffset, 9-yOffset, 2 * 10000.f / 256.f);
    twoDHistogram.Fill(10-xOffset, 8-yOffset, 5 * 10000.f / 256.f);
    twoDHistogram.Fill(10-xOffset, 10-yOffset, 37 * 10000.f / 256.f);
    twoDHistogram.Fill(11-xOffset, 10-yOffset, 17 * 10000.f / 256.f);
    twoDHistogram.Fill(10-xOffset, 11-yOffset, 10 * 10000.f / 256.f);
    twoDHistogram.Fill(11-xOffset, 11-yOffset, 65 * 10000.f / 256.f);
    twoDHistogram.Fill(12-xOffset, 11-yOffset, 10 * 10000.f / 256.f);
    twoDHistogram.Fill(12-xOffset, 12-yOffset, 60 * 10000.f / 256.f);
    twoDHistogram.Fill(12-xOffset, 13-yOffset, 29 * 10000.f / 256.f);
    twoDHistogram.Fill(13-xOffset, 13-yOffset, 44 * 10000.f / 256.f);
    twoDHistogram.Fill(13-xOffset, 14-yOffset, 50 * 10000.f / 256.f);
    twoDHistogram.Fill(14-xOffset, 14-yOffset, 10 * 10000.f / 256.f);
    twoDHistogram.Fill(14-xOffset, 15-yOffset, 78 * 10000.f / 256.f);

    PANDORA_MONITORING_API(DrawPandoraHistogram(this->GetPandora(), twoDHistogram, "COLZ"));

    KerasModel::DataBlock2D dataBlock2D;
    this->HistogramToDataBlock(twoDHistogram, dataBlock2D);
    std::cout << "Input data : " << std::endl;
    dataBlock2D.ShowValues();
    Data1D outputData1D;
    kerasModel.CalculateOutput(&dataBlock2D, outputData1D, this);
    for (int i = 0; i < outputData1D.GetSizeI(); i++)
        std::cout << "Class " << i << ", outcome " << outputData1D.Get(i) << std::endl;

    return STATUS_CODE_SUCCESS;
/*
    CaloHitList caloHitListTruthTrack, caloHitListTruthShower;
    CaloHitList caloHitListDeepTrack, caloHitListDeepShower;

    for (const CaloHit *pTargetCaloHit : wCaloHitList)
    {
        int targetParticleId(0);

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pTargetCaloHit));
            targetParticleId = pMCParticle->GetParticleId();
        }
        catch(...)
        {
            continue;
        }

        TwoDHistogram twoDHistogram(m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f, m_gridSize, -1.f * m_gridDimensions/2.f,  m_gridDimensions/2.f);

        for (const CaloHit *pNeighbourCaloHit : wCaloHitList)
        {
            CartesianVector relativePosition(pNeighbourCaloHit->GetPositionVector() - pTargetCaloHit->GetPositionVector());
            twoDHistogram.Fill(relativePosition.GetX(), relativePosition.GetZ(), pNeighbourCaloHit->GetInputEnergy() * 256.f / 10000.f);
        }

        KerasModel::DataBlock2D dataBlock2D;
        this->HistogramToDataBlock(twoDHistogram, dataBlock2D);
        Data1D outputData1D;
        kerasModel.CalculateOutput(&dataBlock2D, outputData1D);

        if (m_verbose)
        {
//            std::cout << "Calo hit type is " << targetParticleId << std::endl;
            int counter(0);
            for (const float &i : outputData1D)
            {
                std::cout << "Class " << counter << ", outcome " << i << std::endl;
                counter++;
            }
        }

        int isShowerDeepLearning(0);
        int isShowerTruth(0);
        const float x(pTargetCaloHit->GetPositionVector().GetX());
        const float y(pTargetCaloHit->GetPositionVector().GetY());
        const float z(pTargetCaloHit->GetPositionVector().GetZ());

        if (std::abs(targetParticleId) == 11 || targetParticleId == 22)
        {
            caloHitListTruthShower.push_back(pTargetCaloHit);
            isShowerTruth = 0;
        }
        else
        {
            caloHitListTruthTrack.push_back(pTargetCaloHit);
            isShowerTruth = 1;
        }

        if (outputData1D.at(0) > outputData1D.at(1))
        {
            caloHitListDeepShower.push_back(pTargetCaloHit);
            isShowerDeepLearning = 1;
        }
        else
        {
            caloHitListDeepTrack.push_back(pTargetCaloHit);
            isShowerDeepLearning = 0;
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isShowerTruth", isShowerTruth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isShowerDeepLearning", isShowerDeepLearning));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "x", x));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "y", y));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "z", z));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wCaloHitList, "WCaloHits_All", BLACK));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListTruthShower, "WCaloHits_TrueShower", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListTruthTrack, "WCaloHits_TrueTracks", RED));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListDeepShower, "WCaloHits_DeepLearningShower", BLUE));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitListDeepTrack, "WCaloHits_DeepLearningTrack", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CaloHitEventDumpAlgorithm::HistogramToDataBlock(const TwoDHistogram &twoDHistogram, KerasModel::DataBlock2D &dataBlock2D)
{
    Data3D data3D;
    Data2D data2D;
    for (int yBin = 0; yBin < twoDHistogram.GetNBinsY(); yBin++)
//    for (int xBin = 0; xBin < twoDHistogram.GetNBinsX(); xBin++)
    {
        Data1D data1D;
//        for (int yBin = 0; yBin < twoDHistogram.GetNBinsY(); yBin++)
        for (int xBin = 0; xBin < twoDHistogram.GetNBinsX(); xBin++)
        {
            data1D.Append(twoDHistogram.GetBinContent(xBin, yBin) * 256.f / 10000.f ); // I don't know why I did this but it worked, possibly helps the fitting to work with ints
        }
        data2D.Append(data1D);
    }
    data3D.Append(data2D);
    dataBlock2D.SetData(data3D);
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitEventDumpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridSize", m_gridSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GridDimensions", m_gridDimensions));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrainingMode", m_useTrainingMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Verbose", m_verbose));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    if (m_useTrainingMode)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
