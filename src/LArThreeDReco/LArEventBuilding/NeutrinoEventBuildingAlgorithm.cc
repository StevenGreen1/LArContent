/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/NeutrinoEventBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino event building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArEventBuilding/NeutrinoEventBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode NeutrinoEventBuildingAlgorithm::Run()
{
    const PfoList *pPfoList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (NULL == pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoEventBuildingAlgorithm: unable to build neutrinos" << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Assume there is just one neutrino! And assume that it has some daughters! (TODO: Fix these assumptions)
    const ParticleFlowObject *const pNeutrinoPfo = ((1 == pPfoList->size()) ? *(pPfoList->begin()) : NULL);

    if ((NULL == pNeutrinoPfo) || (pNeutrinoPfo->GetVertexList().empty()))
        return STATUS_CODE_FAILURE;

    try
    {
        PfoList daughterPfoList;
        this->GetDaughterPfoList(daughterPfoList);
        this->AddDaughters(pNeutrinoPfo, daughterPfoList);
        this->SetNeutrinoId(pNeutrinoPfo);
    }    
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoEventBuildingAlgorithm::GetDaughterPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_daughterPfoListNames.begin(), iterEnd = m_daughterPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList(NULL);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, *iter, pPfoList))
        {
            pfoList.insert(pPfoList->begin(), pPfoList->end());
        }
        else
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "NeutrinoEventBuildingAlgorithm: pfo list " << *iter << " unavailable." << std::endl;
        }
    }

    if (pfoList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoEventBuildingAlgorithm::AddDaughters(const ParticleFlowObject *const pNeutrinoPfo, const PfoList &daughterPfoList) const
{
    for (PfoList::const_iterator iter = daughterPfoList.begin(), iterEnd = daughterPfoList.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pDaughterPfo(*iter);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pDaughterPfo))
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoEventBuildingAlgorithm::SetNeutrinoId(const ParticleFlowObject *const pNeutrinoPfo) const
{
    if (pNeutrinoPfo->GetDaughterPfoList().empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    unsigned int nPrimaryTwoDHits(0);
    const ParticleFlowObject *pPrimaryDaughter(NULL);

    for (PfoList::const_iterator dIter = pNeutrinoPfo->GetDaughterPfoList().begin(), dIterEnd = pNeutrinoPfo->GetDaughterPfoList().end();
        dIter != dIterEnd; ++dIter)
    {
        const ParticleFlowObject *const pDaughterPfo(*dIter);
        const unsigned int nTwoDHits(this->GetNTwoDHitsInPfo(pDaughterPfo));

        if (!pPrimaryDaughter || (nTwoDHits > nPrimaryTwoDHits))
        {
            nPrimaryTwoDHits = nTwoDHits;
            pPrimaryDaughter = pDaughterPfo;
        }
    }

    if (NULL == pPrimaryDaughter)
      throw StatusCodeException(STATUS_CODE_FAILURE);

    PandoraContentApi::ParticleFlowObject::Metadata metadata;

    if (E_MINUS == std::abs(pPrimaryDaughter->GetParticleId()))
    {
        metadata.m_particleId = NU_E;
    }
    else if (MU_MINUS == std::abs(pPrimaryDaughter->GetParticleId()))
    {
        metadata.m_particleId = NU_MU;
    }

    if (metadata.m_particleId.IsInitialized())
    {
        metadata.m_charge = PdgTable::GetParticleCharge(metadata.m_particleId.Get());
        metadata.m_mass = PdgTable::GetParticleMass(metadata.m_particleId.Get());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AlterMetadata(*this, pNeutrinoPfo, metadata));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoEventBuildingAlgorithm::GetNTwoDHitsInPfo(const ParticleFlowObject *const pPfo) const
{
    unsigned int nTwoDHits(0);

    const ClusterList &pfoClusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = pfoClusterList.begin(), iterEnd = pfoClusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster(*iter);

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            continue;

        nTwoDHits += pCluster->GetNCaloHits();
    }

    const PfoList &daughterList(pPfo->GetDaughterPfoList());

    for (PfoList::const_iterator iter = daughterList.begin(), iterEnd = daughterList.end(); iter != iterEnd; ++iter)
    {
        nTwoDHits += this->GetNTwoDHitsInPfo(*iter);
    }

    return nTwoDHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoEventBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterPfoListNames", m_daughterPfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
