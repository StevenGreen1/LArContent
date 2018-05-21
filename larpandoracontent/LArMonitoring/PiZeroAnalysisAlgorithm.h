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
         *  @param pPfo best matched pfo
         *  @param nMCHits number of mc hits
         *  @param nPfoHits number of pfo hits
         *  @param nSharedHits number of shared hits
         */
        MatchedParticle(const pandora::MCParticle *pMCParticle, const pandora::Pfo *pPfo, const int nMCHits, const int nPfoHits, const int nSharedHits);

        /**
         *  @brief  Calculate reco energy based on calo hits
         */
        void CalculateRecoEnergy();

        /**
         *  @brief Get mc particle
         */
        const pandora::MCParticle *GetMCParticle() const;

        /**
         *  @brief Get pfo
         */
        const pandora::Pfo *GetPfo() const;

        /**
         *  @brief Get m_nMCHits
         */
        int GetNMCHits() const;

        /**
         *  @brief Get m_nPfoHits
         */
        int GetNPfoHits() const;

        /**
         *  @brief Get m_nSharedHits
         */
        int GetSharedHits() const;

        /**
         *  @brief Get reconstructed energy
         */
        float GetRecoEnergy() const;

        /**
         *  @brief Get reconstructed momentum along x
         */
        float GetRecoPx() const;

        /**
         *  @brief Get reconstructed momentum along y
         */
        float GetRecoPy() const;

        /**
         *  @brief Get reconstructed momentum along z
         */
        float GetRecoPz() const;

    private:
        const pandora::MCParticle *m_pMCParticle; ///< Target MCParticle
        const pandora::Pfo        *m_pMatchedPfo; ///< Best matched reconstructed particle
        int                        m_nMCHits;     ///< nHits MC
        int                        m_nPfoHits;    ///< nHits Reco
        int                        m_nSharedHits; ///< nHits Shared
        float                      m_recoEnergy;  ///< Reco energy
        float                      m_recoPx;      ///< Reco Px
        float                      m_recoPy;      ///< Reco Py
        float                      m_recoPz;      ///< Reco Pz
    };

    /**
     *  @brief  AnalysisInfo class
     */
    class AnalysisInfo
    {
    public:
        /**
         *  @brief Constructor
         *
         *  @param photon1 matched particle 1
         *  @param photon2 matched particle 2
         */
        AnalysisInfo(MatchedParticle &photon1, MatchedParticle &photon2);

        /**
         *  @brief Calculate mc and reco pi zero masses
         */
        void CalculatePiZeroMasses();

        /**
         *  @brief Get matched particle 1
         */
        MatchedParticle GetMatch1() const;

        /**
         *  @brief Get matched particle 2
         */
        MatchedParticle GetMatch2() const;

        /**
         *  @brief Get mc pi zero mass
         */
        float GetPiZeroMassMC() const;

        /**
         *  @brief Get mc pi zero energy
         */
        float GetPiZeroEnergyMC() const;

        /**
         *  @brief Get mc pi zero px
         */
        float GetPiZeroPxMC() const;

        /**
         *  @brief Get mc pi zero py
         */
        float GetPiZeroPyMC() const;

        /**
         *  @brief Get mc pi zero pz
         */
        float GetPiZeroPzMC() const;

        /**
         *  @brief Get mc pi zero p
         */
        float GetPiZeroPMC() const;

        /**
         *  @brief Get reco pi zero mass
         */
        float GetPiZeroMassReco() const;

        /**
         *  @brief Get reco pi zero energy
         */
        float GetPiZeroEnergyReco() const;

        /**
         *  @brief Get reco pi zero px
         */
        float GetPiZeroPxReco() const;

        /**
         *  @brief Get reco pi zero py
         */
        float GetPiZeroPyReco() const;

        /**
         *  @brief Get reco pi zero pz
         */
        float GetPiZeroPzReco() const;

        /**
         *  @brief Get reco pi zero p
         */
        float GetPiZeroPReco() const;

    private:
        MatchedParticle m_photon1;           ///< Matched photon 1 info
        MatchedParticle m_photon2;           ///< Matched photon 2 info
        float           m_piZeroMassMC;      ///< MC pion mass
        float           m_piZeroEnergyMC;    ///< MC pion energy
        float           m_piZeroPxMC;        ///< MC pion momentum x
        float           m_piZeroPyMC;        ///< MC pion momentum y
        float           m_piZeroPzMC;        ///< MC pion momentum z
        float           m_piZeroPMC;         ///< MC pion total momentum
        float           m_piZeroMassReco;    ///< Reco pion mass
        float           m_piZeroEnergyReco;  ///< Reco pion energy
        float           m_piZeroPxReco;      ///< Reco pion momentum x
        float           m_piZeroPyReco;      ///< Reco pion momentum y
        float           m_piZeroPzReco;      ///< Reco pion momentum z
        float           m_piZeroPReco;       ///< Reco pion total momentum
    };

    typedef std::vector<AnalysisInfo> AnalysisInfoVector;

    pandora::StatusCode Run();

    /**
     *  @brief  Fill the analysis info containers
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  mcParticleToHitsMap map of mc particle to hits made by mc particle
     *  @param  pfoToHitsMap map of pfo to hits in pfo
     *  @param  mcParticleToPfoHitSharingMap  map of mc particle to pfo - calo hit list pair
     *  @param  analysisInfoVector list of pi zero interactions in event
     */
    void FillAnalysisInfo(const pandora::MCParticleList *const pMCParticleList, LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
        LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap,
        AnalysisInfoVector &analysisInfoVector) const;

    /**
     *  @brief  Fill matched particle info
     *
     *  @param  pMCParticle
     *  @param  mcParticleToHitsMap map of mc particle to hits made by mc particle
     *  @param  pfoToHitsMap map of pfo to hits in pfo
     *  @param  mcParticleToPfoHitSharingMap  map of mc particle to pfo - calo hit list pair
     */
    MatchedParticle FillMatchedParticleInfo(const pandora::MCParticle *const pMCParticle, LArMCParticleHelper::MCContributionMap &mcParticleToHitsMap,
        LArMCParticleHelper::PfoContributionMap &pfoToHitsMap, LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap) const;

    /**
     *  @brief Write variables to a tree
     */
    void WriteToTree(AnalysisInfoVector &analysisInfoVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list
    bool                    m_printAllToScreen;             ///< Whether to print all/raw matching details to screen
    bool                    m_printMatchingToScreen;        ///< Whether to print matching output to screen
    bool                    m_writeToTree;                  ///< Whether to write all/raw matching details to tree
    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file
    bool                    m_visualizePiZero;              ///< Draw event
    int                     m_eventNumber;                  ///< The event number
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::MCParticle *PiZeroAnalysisAlgorithm::MatchedParticle::GetMCParticle() const
{
    return m_pMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Pfo *PiZeroAnalysisAlgorithm::MatchedParticle::GetPfo() const
{
    return m_pMatchedPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNMCHits() const
{
    return m_nMCHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNPfoHits() const
{
    return m_nPfoHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetSharedHits() const
{
    return m_nSharedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetRecoEnergy() const
{
    return m_recoEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetRecoPx() const
{
   return m_recoPx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetRecoPy() const
{
   return m_recoPy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetRecoPz() const
{
   return m_recoPz;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline PiZeroAnalysisAlgorithm::MatchedParticle PiZeroAnalysisAlgorithm::AnalysisInfo::GetMatch1() const
{
    return m_photon1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline PiZeroAnalysisAlgorithm::MatchedParticle PiZeroAnalysisAlgorithm::AnalysisInfo::GetMatch2() const
{
    return m_photon2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroMassMC() const
{
    return m_piZeroMassMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroEnergyMC() const
{
    return m_piZeroEnergyMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPxMC() const
{
    return m_piZeroPxMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPyMC() const
{
    return m_piZeroPyMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPzMC() const
{
    return m_piZeroPzMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPMC() const
{
    return m_piZeroPMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroMassReco() const
{
    return m_piZeroMassReco;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroEnergyReco() const
{
    return m_piZeroEnergyReco;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPxReco() const
{
    return m_piZeroPxReco;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPyReco() const
{
    return m_piZeroPyReco;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPzReco() const
{
    return m_piZeroPzReco;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPReco() const
{
    return m_piZeroPReco;
}

} // namespace lar_content

#endif // LAR_PI_ZERO_ANALYSIS_ALGORITHM_H
