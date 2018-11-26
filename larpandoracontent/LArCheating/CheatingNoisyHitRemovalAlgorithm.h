/**
 *  @file   larpandoracontent/LArCheating/CheatingNoisyHitRemovalAlgorithm.h
 * 
 *  @brief  Header file for the cheating noisy hit removal algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_NOISY_HIT_REMOVAL_ALGORITHM_H
#define LAR_CHEATING_NOISY_HIT_REMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingNoisyHitRemovalAlgorithm::Algorithm class
 */
class CheatingNoisyHitRemovalAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingNoisyHitRemovalAlgorithm() = default;

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName;     ///< Input calo hit list name
    std::string m_outputCaloHitListName;    ///< Output calo hit list name
    std::string m_mcParticleListName;       ///< MC Particle list name
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NOISY_HIT_REMOVAL_ALGORITHM_H
