/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MASTER_ALGORITHM_H
#define LAR_MASTER_ALGORITHM_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/ExternallyConfiguredAlgorithm.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include <unordered_map>

namespace lar_content
{

class StitchingBaseTool;
class CosmicRayTaggingBaseTool;
class SliceIdBaseTool;

typedef std::vector<pandora::CaloHitList> SliceVector;
typedef std::vector<pandora::PfoList> SliceHypotheses;
typedef std::unordered_map<const pandora::ParticleFlowObject*, const pandora::LArTPC*> PfoToLArTPCMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MasterAlgorithm class
 */
class MasterAlgorithm : public pandora::ExternallyConfiguredAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MasterAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~MasterAlgorithm();

    /**
     *  @brief  External steering parameters class
     */
    class ExternalSteeringParameters : public pandora::ExternalParameters
    {
    public:
        pandora::InputBool      m_shouldRunAllHitsCosmicReco;       ///< Whether to run all hits cosmic-ray reconstruction
        pandora::InputBool      m_shouldRunStitching;               ///< Whether to stitch cosmic-ray muons crossing between volumes
        pandora::InputBool      m_shouldRunCosmicHitRemoval;        ///< Whether to remove hits from tagged cosmic-rays
        pandora::InputBool      m_shouldRunSlicing;                 ///< Whether to slice events into separate regions for processing
        pandora::InputBool      m_shouldRunNeutrinoRecoOption;      ///< Whether to run neutrino reconstruction for each slice
        pandora::InputBool      m_shouldRunCosmicRecoOption;        ///< Whether to run cosmic-ray reconstruction for each slice
        pandora::InputBool      m_shouldPerformSliceId;             ///< Whether to identify slices and select most appropriate pfos
        pandora::InputBool      m_printOverallRecoStatus;           ///< Whether to print current operation status messages
    };

    typedef std::unordered_map<const pandora::ParticleFlowObject*, const pandora::LArTPC*> PfoToLArTPCMap;

    /**
     *  @brief  Shift a Pfo hierarchy by a specified x0 value
     *
     *  @param  pPfo the address of the parent pfo
     *  @param  stitchingInfo  the source for additional, local, stitching information
     *  @param  x0 the x0 correction relative to the input pfo
     */
    void ShiftPfoHierarchy(const pandora::ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const;

    /**
     *  @brief  Stitch together a pair of pfos
     *
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     *  @param  stitchingInfo the source for additional, local, stitching information
     */
    void StitchPfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete,
        PfoToLArTPCMap &pfoToLArTPCMap) const;

    // Analysis Variables Added To Test Slicing ID
    class SliceProperties
    {
    public:
        SliceProperties();
        void Reset();

        int m_n2DTargetBeamHits;
        int m_n2DRemnantBeamHits;
        int m_n2DCosmicHits;
        int m_n2DNoMCParticleHits;
        int m_nPfos;
        int m_nUniqueMCPrimaryParticles;
        int m_isBeam;
        pandora::CartesianVector m_centroid;
        pandora::CartesianVector m_majorAxis;
        float m_majorAxisEigenvalue;
        pandora::CartesianVector m_minorAxisOne;
        float m_minorAxisOneEigenvalue;
        pandora::CartesianVector m_minorAxisTwo;
        float m_minorAxisTwoEigenvalue;
        float m_bestSeparation;
        float m_angularSeparationToBeam;
        float m_closestDistance;
    };

    std::string                 m_fileName;                         ///< Root file name
    std::string                 m_treeName;                         ///< Tree name
    int                         m_eventNumber;                      ///< Event number
    float                       m_trueMomentum;                     ///< True momentum of beam sample
    float                       m_momentumAcceptance;               ///< Fractional momentum acceptance

    class Plane
    {
    public:
        /**
         *  @brief Constructor for the plane class.  Equation of plane used here m_a*x + m_b*y + m_c*z + m_d = 0.
         *
         *  @param normal a Cartesian vector that points in a direction that is normal to the plane
         *  @param point a Cartesian vector that corresponds to any point on the plane
         */
        Plane(pandora::CartesianVector &normal, pandora::CartesianVector &point);

       /**
         *  @brief Return the intersection between the plane and a line 
         *
         *  @param a0 point on the line
         *  @param a vector pointing along the line
         */
        pandora::CartesianVector GetLineIntersection(pandora::CartesianVector &a0, pandora::CartesianVector &a) const;

    private:
        float                         m_a;                       ///< Parameter defining a plane
        float                         m_b;                       ///< Parameter defining a plane
        float                         m_c;                       ///< Parameter defining a plane
        float                         m_d;                       ///< Parameter defining a plane
        pandora::CartesianVector      m_unitNormal;              ///< Unit normal to plane
        pandora::CartesianVector      m_point;                   ///< A point on the plane
    };

    typedef std::vector<const Plane*> PlaneVector;
    PlaneVector               m_tpcPlanes;                            ///< Vector of all planes making up global TPC volume
    float                     m_tpcMinX;                              ///< Global TPC volume minimum x extent 
    float                     m_tpcMaxX;                              ///< Global TPC volume maximum x extent
    float                     m_tpcMinY;                              ///< Global TPC volume minimum y extent
    float                     m_tpcMaxY;                              ///< Global TPC volume maximum y extent
    float                     m_tpcMinZ;                              ///< Global TPC volume minimum z extent
    float                     m_tpcMaxZ;                              ///< Global TPC volume maximum z extent
    pandora::CartesianVector  m_beamTPCIntersection;                  ///< Intersection of beam and global TPC volume
    pandora::CartesianVector  m_beamDirection;                        ///< Beam direction

    /**
     *  @brief  Function to find the number of hits belonging to beam/cr for a pfo
     */
    void GetSliceProperties(const pandora::PfoList *const pPfoList, SliceProperties *pSliceProperties);

    /**
     *  @brief Select a given fraction of a slice's calo hits that are closest to the beam spot
     *
     *  @param &inputCaloHitList all calo hits in slice
     *  @param &outputCaloHitList selected calo hits
     *  @param &closestHitToFaceDistance distance of closest hit to beam spot
     */
    pandora::StatusCode GetSelectedCaloHits(pandora::CaloHitList &inputCaloHitList, pandora::CaloHitList &outputCaloHitList, float &closestHitToFaceDistance);

    /**
     *  @brief Find the intercepts of a line with the protoDUNE detector
     *
     *  @param a0 a point on the line in question
     *  @param majorAxis the direction of the line in question
     *  @param interceptOne the first intersection between line and protoDUNE detector 
     *  @param interceptTwo the second intersection between line and protoDUNE detector 
     */
    pandora::StatusCode GetTPCIntercepts(pandora::CartesianVector &a0, pandora::CartesianVector &majorAxis, pandora::CartesianVector &interceptOne,  pandora::CartesianVector &interceptTwo) const;

    /**
     *  @brief Check if a given 3D spacepoint is inside the global TPC volume
     *
     *  @param spacePoint
     */
    bool IsContained(pandora::CartesianVector &spacePoint) const;

private:
    /**
     *  @brief  LArTPCHitList class
     */
    class LArTPCHitList
    {
    public:
        pandora::CaloHitList    m_allHitList;                       ///< The list of all hits originating from a given LArTPC
        pandora::CaloHitList    m_truncatedHitList;                 ///< The list of hits confined within LArTPC boundaries for given beam t0 value
    };

    typedef std::map<unsigned int, LArTPCHitList> VolumeIdToHitListMap;

    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    /**
     *  @brief  Get the mapping from lar tpc volume id to lists of all hits, and truncated hits
     *
     *  @param  volumeIdToHitListMap to receive the populated volume id to hit list map
     */
    pandora::StatusCode GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const;

    /**
     *  @brief  Run the cosmic-ray reconstruction worker instances
     *
     *  @param  volumeIdToHitListMap the volume id to hit list map
     */
    pandora::StatusCode RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const;

    /**
     *  @brief  Recreate cosmic-ray pfos (created by worker instances) in the master instance
     *
     *  @param  pfoToLArTPCMap to receive the populated pfo to lar tpc map
     */
    pandora::StatusCode RecreateCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const;

    /**
     *  @brief  Stitch together cosmic-ray pfos crossing between adjacent lar tpcs
     *
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     */
    pandora::StatusCode StitchCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const;

    /**
     *  @brief  Tag clear, unambiguous cosmic-ray pfos 
     *
     *  @param  clearCosmicRayPfos to receive the list of clear cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous cosmic-ray pfos for further analysis
     */
    pandora::StatusCode TagCosmicRayPfos(pandora::PfoList &clearCosmicRayPfos, pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Run cosmic-ray hit removal, freeing hits in ambiguous pfos for further processing
     *
     *  @param  ambiguousPfos the list of ambiguous cosmic-ray pfos
     */
    pandora::StatusCode RunCosmicRayHitRemoval(const pandora::PfoList &ambiguousPfos) const;

    /**
     *  @brief  Run the event slicing procedures, dividing available hits up into distinct 3D regions
     *
     *  @param  volumeIdToHitListMap the volume id to hit list map
     *  @param  sliceVector to receive the populated slice vector
     */
    pandora::StatusCode RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const;

    /**
     *  @brief  Process each slice under different reconstruction hypotheses
     *
     *  @param  sliceVector the slice vector
     *  @param  nuSliceHypotheses to receive the vector of slice neutrino hypotheses
     *  @param  crSliceHypotheses to receive the vector of slice cosmic-ray hypotheses
     */
    pandora::StatusCode RunSliceReconstruction(SliceVector &sliceVector, SliceHypotheses &nuSliceHypotheses, SliceHypotheses &crSliceHypotheses);

    /**
     *  @brief  Examine slice hypotheses to identify the most appropriate to provide in final event output
     *
     *  @param  nuSliceHypotheses the vector of slice neutrino hypotheses
     *  @param  crSliceHypotheses the vector of slice cosmic-ray hypotheses
     */
    pandora::StatusCode SelectBestSliceHypotheses(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses);

    /**
     *  @brief  Reset all worker instances
     */
    pandora::StatusCode Reset();

    /**
     *  @brief  Copy a specified calo hit to the provided pandora instance
     *
     *  @param  pPandora the address of the target pandora instance
     *  @param  pCaloHit the address of the calo hit
     */
    pandora::StatusCode Copy(const pandora::Pandora *const pPandora, const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Recreate a specified list of pfos in the current pandora instance
     *
     *  @param  inputPfoList the input pfo list
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::PfoList &inputPfoList, pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Recreate a specified pfo in the current pandora instance
     *
     *  @param  pInputPfo the input pfo
     *  @param  pNewParentPfo the new parent of the new output pfo (nullptr if none)
     *  @param  newPfoList to receive the list of new pfos
     */
    pandora::StatusCode Recreate(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject *const pNewParentPfo,
        pandora::PfoList &newPfoList) const;

    /**
     *  @brief  Create a new calo hit in the current pandora instance, based upon the provided input calo hit
     *
     *  @param  pInputCaloHit the address of the input calo hit
     *  @param  pParentCaloHit the address of the parent calo hit
     *
     *  @return the address of the new calo hit
     */
    const pandora::CaloHit *CreateCaloHit(const pandora::CaloHit *const pInputCaloHit, const pandora::CaloHit *const pParentCaloHit) const;

    /**
     *  @brief  Create a new cluster in the current pandora instance, based upon the provided input cluster
     *
     *  @param  pInputCluster the address of the input cluster
     *  @param  newCaloHitList the list of calo hits for the new cluster
     *  @param  newIsolatedCaloHitList the list of isolated calo hits for the new cluster
     *
     *  @return the address of the new cluster
     */
    const pandora::Cluster *CreateCluster(const pandora::Cluster *const pInputCluster, const pandora::CaloHitList &newCaloHitList,
        const pandora::CaloHitList &newIsolatedCaloHitList) const;

    /**
     *  @brief  Create a new vertex in the current pandora instance, based upon the provided input vertex
     *
     *  @param  pInputVertex the address of the input vertex
     *
     *  @return the address of the new vertex
     */
    const pandora::Vertex *CreateVertex(const pandora::Vertex *const pInputVertex) const;

    /**
     *  @brief  Create a new pfo in the current pandora instance, based upon the provided input pfo
     *
     *  @param  pInputPfo the address of the input pfo
     *  @param  newClusterList the list of clusters for the new pfo
     *  @param  newVertexList the list of vertices for the new pfo
     *
     *  @return the address of the new pfo
     */
    const pandora::ParticleFlowObject *CreatePfo(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ClusterList &newClusterList,
        const pandora::VertexList &newVertexList) const;

    /**
     *  @brief  Create a pandora worker instance to handle a single LArTPC
     *
     *  @param  larTPC the lar tpc
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPC &larTPC, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile) const;

    /**
     *  @brief  Create a pandora worker instance to handle a number of LArTPCs
     *
     *  @param  larTPCMap the lar tpc map
     *  @param  gapList the gap list
     *  @param  settingsFile the pandora settings file
     *
     *  @return the address of the pandora instance
     */
    const pandora::Pandora *CreateWorkerInstance(const pandora::LArTPCMap &larTPCMap, const pandora::DetectorGapList &gapList,
        const std::string &settingsFile) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Read settings from external steering parameters block, if present, otherwise from xml as standard
     *
     *  @param  pExternalParameters the address of the external parameters (if present)
     *  @param  inputBool the input boolean value, from the external parameters (if present)
     *  @param  xmlHandle the xml handle
     *  @param  xmlTag the xml tag
     *  @param  outputBool to receive the output boolean value
     */
    pandora::StatusCode ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const pandora::InputBool inputBool,
        const pandora::TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool);

    bool                        m_shouldRunAllHitsCosmicReco;       ///< Whether to run all hits cosmic-ray reconstruction
    bool                        m_shouldRunStitching;               ///< Whether to stitch cosmic-ray muons crossing between volumes
    bool                        m_shouldRunCosmicHitRemoval;        ///< Whether to remove hits from tagged cosmic-rays
    bool                        m_shouldRunSlicing;                 ///< Whether to slice events into separate regions for processing
    bool                        m_shouldRunNeutrinoRecoOption;      ///< Whether to run neutrino reconstruction for each slice
    bool                        m_shouldRunCosmicRecoOption;        ///< Whether to run cosmic-ray reconstruction for each slice
    bool                        m_shouldPerformSliceId;             ///< Whether to identify slices and select most appropriate pfos
    bool                        m_printOverallRecoStatus;           ///< Whether to print current operation status messages
    bool                        m_visualizeOverallRecoStatus;       ///< Whether to display results of current operations

    PandoraInstanceList         m_crWorkerInstances;                ///< The list of cosmic-ray reconstruction worker instances
    const pandora::Pandora     *m_pSlicingWorkerInstance;           ///< The slicing worker instance
    const pandora::Pandora     *m_pSliceNuWorkerInstance;           ///< The per-slice neutrino reconstruction worker instance
    const pandora::Pandora     *m_pSliceCRWorkerInstance;           ///< The per-slice cosmic-ray reconstruction worker instance

    bool                        m_fullWidthCRWorkerWireGaps;        ///< Whether wire-type line gaps in cosmic-ray worker instances should cover all drift time

    typedef std::vector<StitchingBaseTool*> StitchingToolVector;
    typedef std::vector<CosmicRayTaggingBaseTool*> CosmicRayTaggingToolVector;
    typedef std::vector<SliceIdBaseTool*> SliceIdToolVector;

    StitchingToolVector         m_stitchingToolVector;              ///< The stitching tool vector
    CosmicRayTaggingToolVector  m_cosmicRayTaggingToolVector;       ///< The cosmic-ray tagging tool vector
    SliceIdToolVector           m_sliceIdToolVector;                ///< The slice id tool vector

    std::string                 m_filePathEnvironmentVariable;      ///< The environment variable providing a list of paths to xml files
    std::string                 m_crSettingsFile;                   ///< The cosmic-ray reconstruction settings file
    std::string                 m_nuSettingsFile;                   ///< The neutrino reconstruction settings file
    std::string                 m_slicingSettingsFile;              ///< The slicing settings file

    std::string                 m_inputHitListName;                 ///< The input hit list name
    std::string                 m_recreatedPfoListName;             ///< The output recreated pfo list name
    std::string                 m_recreatedClusterListName;         ///< The output recreated cluster list name
    std::string                 m_recreatedVertexListName;          ///< The output recreated vertex list name
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingBaseTool class
 */
class StitchingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pMultiPfoList the list of pfos in multiple lar tpcs
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     */
    virtual void Run(const MasterAlgorithm *const pAlgorithm, const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CosmicRayTaggingBaseTool class
 */
class CosmicRayTaggingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Find the list of ambiguous pfos (could represent cosmic-ray muons or neutrinos)
     *
     *  @param  parentCosmicRayPfos the list of parent cosmic-ray pfos
     *  @param  ambiguousPfos to receive the list of ambiguous pfos
     */
    virtual void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SliceIdBaseTool class
 */
class SliceIdBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Select which reconstruction hypotheses to use; neutrino outcomes or cosmic-ray muon outcomes for each slice
     *
     *  @param  nuSliceHypotheses the parent pfos representing the neutrino outcome for each slice
     *  @param  crSliceHypotheses the parent pfos representing the cosmic-ray muon outcome for each slice
     *  @param  sliceNuPfos to receive the list of selected pfos
     */
    virtual void SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses,
        pandora::PfoList &selectedPfos) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_MASTER_ALGORITHM_H
