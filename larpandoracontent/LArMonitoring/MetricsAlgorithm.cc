/**
 *  @file   larpandoracontent/LArMonitoring/MetricsAlgorithm.cc
 * 
 *  @brief  Implementation of the metrics algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/MetricsAlgorithm.h"

using namespace pandora;

StatusCode MetricsAlgorithm::Run()
{
    // Algorithm code here

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MetricsAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}
