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
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *MetricsAlgorithm::Factory::CreateAlgorithm() const
{
    return new MetricsAlgorithm();
}

#endif // #ifndef METRICS_ALGORITHM_H
