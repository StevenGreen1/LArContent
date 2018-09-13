/**
 *  @file   larpandoracontent/LArMonitoring/CaloHitMonitoringAlgorithm.h
 *
 *  @brief  Header file for the calo hit monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_HIT_MONITORING_ALGORITHM_H
#define LAR_CALO_HIT_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArDeepLearning/LArKerasModel.h"

namespace lar_content
{

/**
 *  @brief  CaloHitMonitoringAlgorithm class
 */
class CaloHitMonitoringAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CaloHitMonitoringAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~CaloHitMonitoringAlgorithm();

    /**
     *  @brief  Function to convert Pandora histogram to Keras data block
     *
     *  @param  twoDHistogram histogram to convert
     *  @param  dataBlock2D data block to populate
     */
    void HistogramToDataBlock(const pandora::TwoDHistogram &twoDHistogram, KerasModel::DataBlock2D &dataBlock2D);

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    KerasModel    m_kerasModel;               ///< Keras model
    std::string   m_pfoListName;              ///< Name of input pfo list
    std::string   m_treeName;                 ///< Monitoring tree name
    std::string   m_fileName;                 ///< Monitoring file name
    int           m_gridSize;                 ///< Number of bins in grid
    float         m_gridDimensions;           ///< Physical dimension of grid
    std::string   m_cnnModelXml;              ///< Xml model file
    std::string   m_cnnModelName;             ///< CNN model name
    int           m_eventNumber;              ///< Event number
};

} // namespace lar_content

#endif // LAR_CALO_HIT_MONITORING_ALGORITHM_H
