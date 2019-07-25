/**
 *  @file   larpandoracontent/LArDeepLearning/CaloHitPropertiesAlgorithm.h
 *
 *  @brief  Header file for the calo hit properties algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_HIT_PROPERTIES_ALGORITHM_H
#define LAR_CALO_HIT_PROPERTIES_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CaloHitPropertiesAlgorithm class
 */
class CaloHitPropertiesAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CaloHitPropertiesAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_caloHitListNames; ///< Names of calo hit lists to show
    bool                    m_trainingMode;     ///< Training mode
};

} // namespace lar_content

#endif // LAR_CALO_HIT_PROPERTIES_ALGORITHM_H
