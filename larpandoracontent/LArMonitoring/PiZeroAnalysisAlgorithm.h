/**
 *  @file   larpandoracontent/LArMonitoring/PiZeroAnalysisAlgorithm.h
 *
 *  @brief  Header file for the pi zero analysis algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PI_ZERO_ANALYSIS_ALGORITHM_H
#define LAR_PI_ZERO_ANALYSIS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <map>

namespace lar_content
{

/**
 *  @brief  PiZeroAnalysisAlgorithm class
 */
class PiZeroAnalysisAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PiZeroAnalysisAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~PiZeroAnalysisAlgorithm();

private:
    /**
     *  @brief  Matched particle class
     */
    class MatchedParticle
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param pMCParticle target mc particle
         */
        MatchedParticle(const pandora::MCParticle *pMCParticle);

        /**
         *  @brief Get mc particle
         */
        const pandora::MCParticle *GetMCParticle();

    private:
        const pandora::MCParticle *m_pMCParticle; ///< Target MCParticle
        const pandora::Pfo        *m_pMatchedPfo; ///< Best matched reconstructed particle
        int                        nMCHits;       ///< nHits MC
        int                        nPFOHits;      ///< nHits Reco
        int                        nSharedHits;   ///< nHits Shared
    };

    /**
     *  @brief  AnalysisInfo class
     */
    class AnalysisInfo
    {
    public:
        AnalysisInfo(MatchedParticle &photon1, MatchedParticle &photon2);

    private:
        MatchedParticle m_photon1; ///< Matched photon 1 info
        MatchedParticle m_photon2; ///< Matched photon 2 info
    };

    pandora::StatusCode Run();

    /**
     *  @brief  Fill the analysis info containers
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  analysisInfo to receive the analysis info
     */
    void FillAnalysisInfo(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, AnalysisInfo &analysisInfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list

    bool                    m_printAllToScreen;             ///< Whether to print all/raw matching details to screen
    bool                    m_printMatchingToScreen;        ///< Whether to print matching output to screen
    bool                    m_writeToTree;                  ///< Whether to write all/raw matching details to tree

    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file

    int                     m_eventNumber;                  ///< The event number
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *PiZeroAnalysisAlgorithm::MatchedParticle::GetMCParticle()
{
    return m_pMCParticle;
}

} // namespace lar_content

#endif // LAR_PI_ZERO_ANALYSIS_ALGORITHM_H
