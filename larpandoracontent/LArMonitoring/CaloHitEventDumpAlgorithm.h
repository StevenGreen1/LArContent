/**
 *  @file   larpandoracontent/LArMonitoring/CaloHitEventDumpAlgorithm.h
 *
 *  @brief  Header file for the calo hit event dump algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_HIT_EVENT_DUMP_ALGORITHM_H
#define LAR_CALO_HIT_EVENT_DUMP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArDeepLearning/LArKerasModel.h"

namespace lar_content
{

/**
 *  @brief  CaloHitEventDumpAlgorithm class
 */
class CaloHitEventDumpAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CaloHitEventDumpAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~CaloHitEventDumpAlgorithm();

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

    std::string   m_caloHitListName;          ///< Name of input calo hit list
    bool          m_verbose;                  ///< Print model information
    std::string   m_treeName;                 ///< Monitoring tree name
    std::string   m_fileName;                 ///< Monitoring file name
    int           m_gridSize;                 ///< Number of bins in grid
    float         m_gridDimensions;           ///< Physical dimension of grid
    bool          m_useTrainingMode;          ///< Should use training mode. If true, training examples will be written to the output file
    std::string   m_trainingOutputFile;       ///< Output file name for training examples
};

} // namespace lar_content

#endif // LAR_CALO_HIT_EVENT_DUMP_ALGORITHM_H
