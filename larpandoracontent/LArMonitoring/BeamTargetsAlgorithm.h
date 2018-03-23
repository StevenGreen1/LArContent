/**
 *  @file   larpandoracontent/LArMonitoring/BeamTargetsAlgorithm.h
 * 
 *  @brief  Header file for the beam target algorithm class.
 * 
 *  $Log: $
 */
#ifndef BEAM_TARGETS_ALGORITHM_H
#define BEAM_TARGETS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  BeamTargetsAlgorithm class
 */
class BeamTargetsAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Get beam MC particles from given list
     */
    void GetBeamMCParticles(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleList &beamMCParticleList);

void GetCaloHitToTargetMap(const pandora::CaloHitList &caloHitList, const pandora::MCParticle *pMCPrimary, LArMCParticleHelper::MCContributionMap &mcToCaloHitMap);
//void GetVisibleDaughter(const pandora::MCParticle *pParent, const pandora::MCParticle *pTarget);
void IsVisibleTarget(const pandora::MCParticle *pParticle, const pandora::MCParticleList &targets, const pandora::MCParticle *&pTarget);

void AddCaloHitToMap(const pandora::CaloHit *pCaloHit, const pandora::MCParticle *pTarget, LArMCParticleHelper::MCContributionMap &mcToCaloHitMap);

void GetAllDownstreamMCParticles(const pandora::MCParticleList &inputMCParticleList, pandora::MCParticleList &outputMCParticleList);
void GetAllDownstreamMCParticles(const pandora::MCParticle *const pMCParticle, pandora::MCParticleList &outputMCParticleList);
void FilterMCParticles(const pandora::MCParticleList &allMCParticles, pandora::MCParticleList &filteredMCParticles, const pandora::CaloHitList &caloHitList);
bool CheckParentInList(const pandora::MCParticle *pMCParticle, const pandora::MCParticleList &filteredMCParticles);


    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BeamTargetsAlgorithm::Factory::CreateAlgorithm() const
{
    return new BeamTargetsAlgorithm();
}

} // namespace lar_content

#endif // #ifndef BEAM_TARGETS_ALGORITHM_H
