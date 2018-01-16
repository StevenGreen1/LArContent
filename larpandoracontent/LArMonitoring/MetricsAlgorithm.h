/**
 *  @file   larpandoracontent/LArMonitoring/MetricsAlgorithm.h
 * 
 *  @brief  Header file for the metrics algorithm class.
 * 
 *  $Log: $
 */
#ifndef METRICS_ALGORITHM_H
#define METRICS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

/**
 *  @brief  MetricsAlgorithm class
 */
class MetricsAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  Constructor
     */
    MetricsAlgorithm();

    /**
     *  Destructor
     */
    ~MetricsAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables 
    std::string           m_fileName;              ///< Analysis root file name
    std::string           m_treeName;              ///< Analysis root tree name
    std::string           m_mcParticleListName;    ///< MCParticle List Name
    std::string           m_pfoListName;           ///< PFO List Name
};

#endif // #ifndef METRICS_ALGORITHM_H
