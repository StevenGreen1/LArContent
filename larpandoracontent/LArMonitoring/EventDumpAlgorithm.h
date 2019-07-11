/**
 *  @file   larpandoracontent/LArMonitoring/EventDumpAlgorithm.h
 *
 *  @brief  Header file for the event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_DUMP_ALGORITHM_H
#define LAR_EVENT_DUMP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <map>

namespace lar_content
{

/**
 *  @brief  EventDumpAlgorithm class
 */
class EventDumpAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventDumpAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventDumpAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list
    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file
};

} // namespace lar_content

#endif // LAR_EVENT_DUMP_ALGORITHM_H
