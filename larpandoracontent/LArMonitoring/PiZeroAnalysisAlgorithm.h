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
         *  @param mcCaloHitList mc hits
         *  @param pfoCaloHitList pfo hits
         *  @param sharedCaloHitList shared hits
         *  @param hitToGeV calibration
         */
        MatchedParticle(const pandora::MCParticle *pMCParticle, const pandora::Pfo *pPfo, pandora::CaloHitList &mcCaloHitList, pandora::CaloHitList &pfoCaloHitList, pandora::CaloHitList &sharedCaloHitList, const float hitToGeV);

        /**
         *  @brief  Count hits in calo hit list by type
         *
         *  @param  caloHitList
         *  @param  nHitsU
         *  @param  nHitsV
         *  @param  nHitsW
         */
        void CountHits(pandora::CaloHitList &caloHitList, int &nHitsU, int &nHitsV, int &nHitsW);

        /**
         *  @brief  Calculate reco energy based on calo hits
         *
         *  @param mcCaloHitList mc hits
         */
        void CalculateRecoEnergy(pandora::CaloHitList &mcCaloHitList);

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
         *  @brief Get m_nMCHitsU
         */
        int GetNMCHitsU() const;

        /**
         *  @brief Get m_nMCHitsV
         */
        int GetNMCHitsV() const;

        /**
         *  @brief Get m_nMCHitsW
         */
        int GetNMCHitsW() const;

        /**
         *  @brief Get m_nPfoHits
         */
        int GetNPfoHits() const;

        /**
         *  @brief Get m_nPfoHitsU
         */
        int GetNPfoHitsU() const;

        /**
         *  @brief Get m_nPfoHitsV
         */
        int GetNPfoHitsV() const;

        /**
         *  @brief Get m_nPfoHitsW
         */
        int GetNPfoHitsW() const;

        /**
         *  @brief Get m_nSharedHits
         */
        int GetSharedHits() const;

        /**
         *  @brief Get m_nSharedHitsU
         */
        int GetSharedHitsU() const;

        /**
         *  @brief Get m_nSharedHitsV
         */
        int GetSharedHitsV() const;

        /**
         *  @brief Get m_nSharedHitsW
         */
        int GetSharedHitsW() const;

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

        /**
         *  @brief Get reconstructed energy with cheated pat rec
         */
        float GetCheatedPatRecEnergy() const;

        /**
         *  @brief Get reconstructed momentum along x with cheated pat rec
         */
        float GetCheatedPatRecPx() const;

        /**
         *  @brief Get reconstructed momentum along y with cheated pat rec
         */
        float GetCheatedPatRecPy() const;

        /**
         *  @brief Get reconstructed momentum along z with cheated pat rec
         */
        float GetCheatedPatRecPz() const;

        /**
         *  @brief Get reconstructed momentum along x with cheated pat rec, reco direction
         */
        float GetCheatedPatRecRecoDirPx() const;

        /**
         *  @brief Get reconstructed momentum along y with cheated pat rec, reco direction
         */
        float GetCheatedPatRecRecoDirPy() const;

        /**
         *  @brief Get reconstructed momentum along z with cheated pat rec, reco direction
         */
        float GetCheatedPatRecRecoDirPz() const;

    private:
        const pandora::MCParticle *m_pMCParticle;            ///< Target MCParticle
        const pandora::Pfo        *m_pMatchedPfo;            ///< Best matched reconstructed particle
        int                        m_nMCHits;                ///< nHits MC
        int                        m_nMCHitsU;               ///< nHits MC U
        int                        m_nMCHitsV;               ///< nHits MC V
        int                        m_nMCHitsW;               ///< nHits MC W
        int                        m_nPfoHits;               ///< nHits Reco
        int                        m_nPfoHitsU;              ///< nHits Reco U
        int                        m_nPfoHitsV;              ///< nHits Reco V
        int                        m_nPfoHitsW;              ///< nHits Reco W
        int                        m_nSharedHits;            ///< nHits Shared
        int                        m_nSharedHitsU;           ///< nHits Shared U
        int                        m_nSharedHitsV;           ///< nHits Shared V
        int                        m_nSharedHitsW;           ///< nHits Shared W
        float                      m_recoEnergy;             ///< Reco energy
        float                      m_recoPx;                 ///< Reco Px
        float                      m_recoPy;                 ///< Reco Py
        float                      m_recoPz;                 ///< Reco Pz
        float                      m_cheatedPatRecEnergy;    ///< Reconstructed hits cheated pattern recognition energy
        float                      m_cheatedPatRecPx;        ///< Reconstructed hits cheated pattern recognition Px
        float                      m_cheatedPatRecPy;        ///< Reconstructed hits cheated pattern recognition Py
        float                      m_cheatedPatRecPz;        ///< Reconstructed hits cheated pattern recognition Pz
        float                      m_cheatedPatRecRecoDirPx; ///< Reconstructed hits cheated pattern recognition reconstructed directon Px
        float                      m_cheatedPatRecRecoDirPy; ///< Reconstructed hits cheated pattern recognition reconstructed directon Py
        float                      m_cheatedPatRecRecoDirPz; ///< Reconstructed hits cheated pattern recognition reconstructed directon Pz
        float                      m_hitToGeV;               ///< Calibration
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
         *  @brief Get cheated pat rec pi zero mass
         */
        float GetPiZeroMassCheatedPatRec() const;

        /**
         *  @brief Get cheated pat rec pi zero energy
         */
        float GetPiZeroEnergyCheatedPatRec() const;

        /**
         *  @brief Get cheated pat rec pi zero px
         */
        float GetPiZeroPxCheatedPatRec() const;

        /**
         *  @brief Get cheated pat rec pi zero py
         */
        float GetPiZeroPyCheatedPatRec() const;

        /**
         *  @brief Get cheated pat rec pi zero pz
         */
        float GetPiZeroPzCheatedPatRec() const;

        /**
         *  @brief Get cheated pat rec pi zero p
         */
        float GetPiZeroPCheatedPatRec() const;

        /**
         *  @brief Get cheated pat rec, reco direction pi zero mass
         */
        float GetPiZeroMassCheatedPatRecRecoDir() const;

        /**
         *  @brief Get cheated pat rec, reco direction pi zero energy
         */
        float GetPiZeroEnergyCheatedPatRecRecoDir() const;

        /**
         *  @brief Get cheated pat rec, reco direction pi zero px
         */
        float GetPiZeroPxCheatedPatRecRecoDir() const;

        /**
         *  @brief Get cheated pat rec, reco direction pi zero py
         */
        float GetPiZeroPyCheatedPatRecRecoDir() const;

        /**
         *  @brief Get cheated pat rec, reco direction pi zero pz
         */
        float GetPiZeroPzCheatedPatRecRecoDir() const;

        /**
         *  @brief Get cheated pat rec, reco direction pi zero p
         */
        float GetPiZeroPCheatedPatRecRecoDir() const;

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
        MatchedParticle m_photon1;                           ///< Matched photon 1 info
        MatchedParticle m_photon2;                           ///< Matched photon 2 info

        float           m_piZeroMassMC;                      ///< MC pion mass
        float           m_piZeroEnergyMC;                    ///< MC pion energy
        float           m_piZeroPxMC;                        ///< MC pion momentum x
        float           m_piZeroPyMC;                        ///< MC pion momentum y
        float           m_piZeroPzMC;                        ///< MC pion momentum z
        float           m_piZeroPMC;                         ///< MC pion total momentum

        float           m_piZeroMassCheatedPatRec;           ///< Cheated pat rec pion mass
        float           m_piZeroEnergyCheatedPatRec;         ///< Cheated pat rec pion energy
        float           m_piZeroPxCheatedPatRec;             ///< Cheated pat rec pion momentum x
        float           m_piZeroPyCheatedPatRec;             ///< Cheated pat rec pion momentum y
        float           m_piZeroPzCheatedPatRec;             ///< Cheated pat rec pion momentum z
        float           m_piZeroPCheatedPatRec;              ///< Cheated pat rec pion total momentum

        float           m_piZeroMassCheatedPatRecRecoDir;    ///< Cheated pat rec, reco direction pion mass
        float           m_piZeroEnergyCheatedPatRecRecoDir;  ///< Cheated pat rec, reco direction pion energy
        float           m_piZeroPxCheatedPatRecRecoDir;      ///< Cheated pat rec, reco direction pion momentum x
        float           m_piZeroPyCheatedPatRecRecoDir;      ///< Cheated pat rec, reco direction pion momentum y
        float           m_piZeroPzCheatedPatRecRecoDir;      ///< Cheated pat rec, reco direction pion momentum z
        float           m_piZeroPCheatedPatRecRecoDir;       ///< Cheated pat rec, reco direction pion total momentum

        float           m_piZeroMassReco;                    ///< Reco pion mass
        float           m_piZeroEnergyReco;                  ///< Reco pion energy
        float           m_piZeroPxReco;                      ///< Reco pion momentum x
        float           m_piZeroPyReco;                      ///< Reco pion momentum y
        float           m_piZeroPzReco;                      ///< Reco pion momentum z
        float           m_piZeroPReco;                       ///< Reco pion total momentum
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
     *  @brief  Add calo hit to relevant u, v or w list
     *
     *  @param  caloHitListU
     *  @param  caloHitListV
     *  @param  caloHitListW
     *  @param  pCaloHit to add
     */
    void AddCaloHit(pandora::CaloHitList &caloHitListU, pandora::CaloHitList &caloHitListV, pandora::CaloHitList &caloHitListW, const pandora::CaloHit *pCaloHit) const;

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
    float                   m_hitToGeV;                     ///< Calibration
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

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNMCHitsU() const
{
    return m_nMCHitsU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNMCHitsV() const
{
    return m_nMCHitsV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNMCHitsW() const
{
    return m_nMCHitsW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNPfoHits() const
{
    return m_nPfoHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNPfoHitsU() const
{
    return m_nPfoHitsU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNPfoHitsV() const
{
    return m_nPfoHitsV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetNPfoHitsW() const
{
    return m_nPfoHitsW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetSharedHits() const
{
    return m_nSharedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetSharedHitsU() const
{
    return m_nSharedHitsU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetSharedHitsV() const
{
    return m_nSharedHitsV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PiZeroAnalysisAlgorithm::MatchedParticle::GetSharedHitsW() const
{
    return m_nSharedHitsW;
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

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecEnergy() const
{
    return m_cheatedPatRecEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecPx() const
{
   return m_cheatedPatRecPx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecPy() const
{
   return m_cheatedPatRecPy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecPz() const
{
   return m_cheatedPatRecPz;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecRecoDirPx() const
{
   return m_cheatedPatRecRecoDirPx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecRecoDirPy() const
{
   return m_cheatedPatRecRecoDirPy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::MatchedParticle::GetCheatedPatRecRecoDirPz() const
{
   return m_cheatedPatRecRecoDirPz;
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

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroMassCheatedPatRec() const
{
    return m_piZeroMassCheatedPatRec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroEnergyCheatedPatRec() const
{
    return m_piZeroEnergyCheatedPatRec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPxCheatedPatRec() const
{
    return m_piZeroPxCheatedPatRec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPyCheatedPatRec() const
{
    return m_piZeroPyCheatedPatRec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPzCheatedPatRec() const
{
    return m_piZeroPzCheatedPatRec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPCheatedPatRec() const
{
    return m_piZeroPCheatedPatRec;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroMassCheatedPatRecRecoDir() const
{
    return m_piZeroMassCheatedPatRecRecoDir;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroEnergyCheatedPatRecRecoDir() const
{
    return m_piZeroEnergyCheatedPatRecRecoDir;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPxCheatedPatRecRecoDir() const
{
    return m_piZeroPxCheatedPatRecRecoDir;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPyCheatedPatRecRecoDir() const
{
    return m_piZeroPyCheatedPatRecRecoDir;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPzCheatedPatRecRecoDir() const
{
    return m_piZeroPzCheatedPatRecRecoDir;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PiZeroAnalysisAlgorithm::AnalysisInfo::GetPiZeroPCheatedPatRecRecoDir() const
{
    return m_piZeroPCheatedPatRecRecoDir;
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