/**
 *  @file   larpandoracontent/LArDeepLearning/CaloHitPropertiesAlgorithm.cc
 *
 *  @brief  Implementation of the calo hit properties algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArDeepLearning/CaloHitPropertiesAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CaloHitPropertiesAlgorithm::CaloHitPropertiesAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitPropertiesAlgorithm::Run()
{
//    const CaloHitList *pCaloHitList = nullptr;
//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CaloHitPropertiesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
