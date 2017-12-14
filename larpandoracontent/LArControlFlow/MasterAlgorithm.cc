/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MasterAlgorithm::MasterAlgorithm() :
    m_fileName("SliceAnalysis.root"),
    m_treeName("SliceAnalysisTree"),
    m_eventNumber(0),
    m_trueMomentum(-1.f),
    m_momentumAcceptance(20.f),
    m_beamTPCIntersection(CartesianVector(0.f, 0.f, 0.f)),
    m_beamDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(true),
    m_visualizeOverallRecoStatus(false),
    m_pSlicingWorkerInstance(nullptr),
    m_pSliceNuWorkerInstance(nullptr),
    m_pSliceCRWorkerInstance(nullptr),
    m_fullWidthCRWorkerWireGaps(true),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MasterAlgorithm::~MasterAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::ShiftPfoHierarchy(const ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const
{
    if (!pParentPfo->GetParentPfoList().empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PfoToLArTPCMap::const_iterator larTPCIter(pfoToLArTPCMap.find(pParentPfo));

    if (pfoToLArTPCMap.end() == larTPCIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    PfoList pfoList;
    LArPfoHelper::GetAllDownstreamPfos(pParentPfo, pfoList);
    const float signedX0(larTPCIter->second->IsDriftInPositiveX() ? -x0 : x0);

    if (m_visualizeOverallRecoStatus)
    {
        std::cout << "ShiftPfoHierarchy: signedX0 " << signedX0 << std::endl;
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "BeforeShiftCRPfos", GREEN));
    }

    for (const ParticleFlowObject *const pDaughterPfo : pfoList)
    {
        CaloHitList caloHitList;
        for (const Cluster *const pCluster : pDaughterPfo->GetClusterList())
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
            caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        }

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            PandoraContentApi::CaloHit::Metadata metadata;
            metadata.m_x0 = signedX0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pCaloHit, metadata));
        }

        for (const Vertex *const pVertex : pDaughterPfo->GetVertexList())
        {
            PandoraContentApi::Vertex::Metadata metadata;
            metadata.m_x0 = signedX0;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::AlterMetadata(*this, pVertex, metadata));
        }
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "AfterShiftCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::StitchPfos(const ParticleFlowObject *const pPfoToEnlarge, const ParticleFlowObject *const pPfoToDelete,
    PfoToLArTPCMap &pfoToLArTPCMap) const
{
    if (pPfoToEnlarge == pPfoToDelete)
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN Remove both pfos from stitching information here - could cause problems if stitching across multiple TPCs
    pfoToLArTPCMap.erase(pPfoToEnlarge);
    pfoToLArTPCMap.erase(pPfoToDelete);

    const PfoList daughterPfos(pPfoToDelete->GetDaughterPfoList());
    const ClusterVector daughterClusters(pPfoToDelete->GetClusterList().begin(), pPfoToDelete->GetClusterList().end());
    const VertexVector daughterVertices(pPfoToDelete->GetVertexList().begin(), pPfoToDelete->GetVertexList().end());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete, m_recreatedPfoListName));

    for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoToEnlarge, pDaughterPfo));

    for (const  Vertex *const pDaughterVertex : daughterVertices)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pDaughterVertex, m_recreatedVertexListName));

    for (const Cluster *const pDaughterCluster : daughterClusters)
    {
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *pParentCluster(PfoMopUpBaseAlgorithm::GetParentCluster(pPfoToEnlarge->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
                m_recreatedClusterListName, m_recreatedClusterListName));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToEnlarge, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Initialize()
{
    try
    {
        const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
        const DetectorGapList &gapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

        for (const LArTPCMap::value_type &mapEntry : larTPCMap)
            m_crWorkerInstances.push_back(this->CreateWorkerInstance(*(mapEntry.second), gapList, m_crSettingsFile));

        if (m_shouldRunSlicing)
            m_pSlicingWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_slicingSettingsFile);

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile);

        if (m_shouldRunCosmicRecoOption)
            m_pSliceCRWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "MasterAlgorithm: Exception during initialization of worker instances " << statusCodeException.ToString() << std::endl;
        return statusCodeException.GetStatusCode();
    }

    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());
    
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {   
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }
    
    m_tpcMinX = parentMinX;
    m_tpcMaxX = parentMaxX;
    m_tpcMinY = parentMinY;
    m_tpcMaxY = parentMaxY;
    m_tpcMinZ = parentMinZ;
    m_tpcMaxZ = parentMaxZ;
    
    CartesianVector normalTop(0,0,1), pointTop(0,0,m_tpcMaxZ);
    const Plane *pPlaneTop = new Plane(normalTop, pointTop);
    m_tpcPlanes.push_back(pPlaneTop);
    
    CartesianVector normalBottom(0,0,-1), pointBottom(0,0,m_tpcMinZ);
    const Plane *pPlaneBottom = new Plane(normalBottom, pointBottom);
    m_tpcPlanes.push_back(pPlaneBottom);
    
    CartesianVector normalRight(1,0,0), pointRight(m_tpcMaxX,0,0);
    const Plane *pPlaneRight = new Plane(normalRight, pointRight);
    m_tpcPlanes.push_back(pPlaneRight);
    
    CartesianVector normalLeft(-1,0,0), pointLeft(m_tpcMinX,0,0);
    const Plane *pPlaneLeft = new Plane(normalLeft, pointLeft);
    m_tpcPlanes.push_back(pPlaneLeft);
    
    CartesianVector normalBack(0,1,0), pointBack(0,m_tpcMaxY,0);
    const Plane *pPlaneBack = new Plane(normalBack, pointBack);
    m_tpcPlanes.push_back(pPlaneBack);
    
    CartesianVector normalFront(0,-1,0), pointFront(0,m_tpcMinY,0);
    const Plane *pPlaneFront = new Plane(normalFront, pointFront);
    m_tpcPlanes.push_back(pPlaneFront);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Run()
{
    VolumeIdToHitListMap volumeIdToHitListMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVolumeIdToHitListMap(volumeIdToHitListMap));

    if (m_shouldRunAllHitsCosmicReco)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayReconstruction(volumeIdToHitListMap));

        PfoToLArTPCMap pfoToLArTPCMap;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RecreateCosmicRayPfos(pfoToLArTPCMap));

        if (m_shouldRunStitching)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->StitchCosmicRayPfos(pfoToLArTPCMap));
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        PfoList clearCosmicRayPfos, ambiguousPfos;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->TagCosmicRayPfos(clearCosmicRayPfos, ambiguousPfos));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunCosmicRayHitRemoval(ambiguousPfos));
    }

    SliceVector sliceVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSlicing(volumeIdToHitListMap, sliceVector));

    if (m_shouldRunNeutrinoRecoOption || m_shouldRunCosmicRecoOption)
    {
        SliceHypotheses nuSliceHypotheses, crSliceHypotheses;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunSliceReconstruction(sliceVector, nuSliceHypotheses, crSliceHypotheses));
        const PfoList pfoBeamList(nuSliceHypotheses.at(0));
        const PfoList pfoCRList(crSliceHypotheses.at(0));
        std::cout << "m_shouldRunNeutrinoRecoOption " << m_shouldRunNeutrinoRecoOption << std::endl;
        std::cout << "pfoBeamList.size() " << pfoBeamList.size() << std::endl;
        std::cout << "pfoBeamList.front()->GetClusterList().size() " << pfoBeamList.front()->GetClusterList().size() << std::endl;
        std::cout << "pfoCRList.size() " << pfoCRList.size() << std::endl;
        std::cout << "pfoCRList.front()->GetClusterList().size() " << pfoCRList.front()->GetClusterList().size() << std::endl;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SelectBestSliceHypotheses(nuSliceHypotheses, crSliceHypotheses));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Reset());
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::GetVolumeIdToHitListMap(VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const unsigned int nLArTPCs(larTPCMap.size());

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputHitListName, pCaloHitList));

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit*>(pCaloHit));

        if (!pLArCaloHit && (1 != nLArTPCs))
            return STATUS_CODE_INVALID_PARAMETER;

        const unsigned int volumeId(pLArCaloHit ? pLArCaloHit->GetLArTPCVolumeId() : 0);
        const LArTPC *const pLArTPC(larTPCMap.at(volumeId));

        LArTPCHitList &larTPCHitList(volumeIdToHitListMap[volumeId]);
        larTPCHitList.m_allHitList.push_back(pCaloHit);

        if (((pCaloHit->GetPositionVector().GetX() >= (pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX())) &&
            (pCaloHit->GetPositionVector().GetX() <= (pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX()))))
        {
            larTPCHitList.m_truncatedHitList.push_back(pCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayReconstruction(const VolumeIdToHitListMap &volumeIdToHitListMap) const
{
    unsigned int workerCounter(0);

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());
        VolumeIdToHitListMap::const_iterator iter(volumeIdToHitListMap.find(larTPC.GetLArTPCVolumeId()));

        if (volumeIdToHitListMap.end() == iter)
            continue;

        for (const CaloHit *const pCaloHit : iter->second.m_allHitList)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pCRWorker, pCaloHit));
        }

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance " << ++workerCounter << " of " << m_crWorkerInstances.size() << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pCRWorker));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RecreateCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const PfoList *pCRPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pCRWorker, pCRPfos));

        PfoList newPfoList;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(*pCRPfos, newPfoList));

        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());

        for (const Pfo *const pNewPfo : newPfoList)
            pfoToLArTPCMap[pNewPfo] = &larTPC;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::StitchCosmicRayPfos(PfoToLArTPCMap &pfoToLArTPCMap) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "RecreatedCRPfos", GREEN));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    for (StitchingBaseTool *const pStitchingTool : m_stitchingToolVector)
        pStitchingTool->Run(this, pRecreatedCRPfos, pfoToLArTPCMap);

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "AfterStitchingCRPfos", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::TagCosmicRayPfos(PfoList &clearCosmicRayPfos, PfoList &ambiguousPfos) const
{
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));

    PfoList parentCosmicRayPfos;
    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        if (pPfo->GetParentPfoList().empty())
            parentCosmicRayPfos.push_back(pPfo);
    }

    for (CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool : m_cosmicRayTaggingToolVector)
        pCosmicRayTaggingTool->FindAmbiguousPfos(parentCosmicRayPfos, ambiguousPfos);

    for (const Pfo *const pPfo : parentCosmicRayPfos)
    {
        if (ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo))
            clearCosmicRayPfos.push_back(pPfo);
    }

    if (m_visualizeOverallRecoStatus)
    {
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &clearCosmicRayPfos, "ClearCRPfos", RED));
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &ambiguousPfos, "AmbiguousCRPfos", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunCosmicRayHitRemoval(const PfoList &ambiguousPfos) const
{
    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(ambiguousPfos, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        const ClusterList clusterList(pPfoToDelete->GetClusterList());
        const VertexList vertexList(pPfoToDelete->GetVertexList());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &vertexList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSlicing(const VolumeIdToHitListMap &volumeIdToHitListMap, SliceVector &sliceVector) const
{
    for (const VolumeIdToHitListMap::value_type &mapEntry : volumeIdToHitListMap)
    {
        for (const CaloHit *const pCaloHit : mapEntry.second.m_truncatedHitList)
        {
            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                continue;

            const HitType hitType(pCaloHit->GetHitType());
            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                continue;

            if (m_shouldRunSlicing)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSlicingWorkerInstance, pCaloHit));
            }
            else
            {
                if (sliceVector.empty()) sliceVector.push_back(CaloHitList());
                sliceVector.back().push_back(pCaloHit);
            }
        }
    }

    if (m_shouldRunSlicing)
    {
        if (m_printOverallRecoStatus)
            std::cout << "Running slicing worker instance" << std::endl;

        const PfoList *pSlicePfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSlicingWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSlicingWorkerInstance, pSlicePfos));

        if (m_visualizeOverallRecoStatus)
        {
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pSlicePfos, "OnePfoPerSlice", BLUE));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        for (const Pfo *const pSlicePfo : *pSlicePfos)
        {
            sliceVector.push_back(CaloHitList());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_U, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_V, sliceVector.back());
            LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_W, sliceVector.back());
        }
    }

    if (m_printOverallRecoStatus)
        std::cout << "Identified " << sliceVector.size() << " slice(s)" << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::RunSliceReconstruction(SliceVector &sliceVector, SliceHypotheses &nuSliceHypotheses, SliceHypotheses &crSliceHypotheses) 
{
    unsigned int sliceCounter(0);

    for (const CaloHitList &sliceHits : sliceVector)
    {
        for (const CaloHit *const pSliceCaloHit : sliceHits)
        {
            // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
            const CaloHit *const pCaloHitInMaster(m_shouldRunSlicing ? static_cast<const CaloHit*>(pSliceCaloHit->GetParentAddress()) : pSliceCaloHit);

            if (m_shouldRunNeutrinoRecoOption)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceNuWorkerInstance, pCaloHitInMaster));
            }

            if (m_shouldRunCosmicRecoOption)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceCRWorkerInstance, pCaloHitInMaster));
            }
        }

        ++sliceCounter;

        if (m_shouldRunNeutrinoRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running nu worker instance for slice " << sliceCounter << " of " << sliceVector.size() << std::endl;

            const PfoList *pSliceNuPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceNuWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceNuWorkerInstance, pSliceNuPfos));
            nuSliceHypotheses.push_back(*pSliceNuPfos);
        }

        if (m_shouldRunCosmicRecoOption)
        {
            if (m_printOverallRecoStatus)
                std::cout << "Running cr worker instance for slice " << sliceCounter << " of " << sliceVector.size() << std::endl;

            const PfoList *pSliceCRPfos(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceCRWorkerInstance));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceCRWorkerInstance, pSliceCRPfos));
            crSliceHypotheses.push_back(*pSliceCRPfos);
        }
    }

    if (m_shouldRunNeutrinoRecoOption && m_shouldRunCosmicRecoOption && (nuSliceHypotheses.size() != crSliceHypotheses.size()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::SelectBestSliceHypotheses(const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses) 
{
    if (m_printOverallRecoStatus)
        std::cout << "Select best slice hypotheses" << std::endl;    

    PfoList selectedSlicePfos;

//~~~~
    // Variables, number of PFOs, energy, momentum, direction 
    m_eventNumber++;

    int sliceIndex(0), nSlices(beamSliceHypotheses.size());

    for (sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList pfoBeamList(beamSliceHypotheses.at(sliceIndex));
        const PfoList pfoCRList(crSliceHypotheses.at(sliceIndex));

        SliceProperties *pSlicePropertiesBeam = new SliceProperties();
        SliceProperties *pSlicePropertiesCosmic = new SliceProperties();

//        std::cout << "Beam Slice" << std::endl;
        this->GetSliceProperties(&pfoBeamList, pSlicePropertiesBeam);

//        std::cout << "Cosmic Slice" << std::endl;
        this->GetSliceProperties(&pfoCRList, pSlicePropertiesCosmic);
/*
        PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoBeamList, "BeamSlice", RED, true, true);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoCRList, "CRSlice", BLUE, true, true);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
*/

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "sliceIndex", sliceIndex));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isBeam", pSlicePropertiesBeam->m_isBeam));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DTargetBeamHitsBeamSlice", pSlicePropertiesBeam->m_n2DTargetBeamHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DRemnantBeamHitsBeamSlice", pSlicePropertiesBeam->m_n2DRemnantBeamHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DCosmicBeamHitsBeamSlice", pSlicePropertiesBeam->m_n2DCosmicHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DNoMCParticleHitsBeamSlice", pSlicePropertiesBeam->m_n2DNoMCParticleHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfosBeamSlice", pSlicePropertiesBeam->m_nPfos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nUniqueMCPrimaryParticlesBeamSlice", pSlicePropertiesBeam->m_nUniqueMCPrimaryParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroidXBeamSlice", pSlicePropertiesBeam->m_centroid.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroidYBeamSlice", pSlicePropertiesBeam->m_centroid.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroidZBeamSlice", pSlicePropertiesBeam->m_centroid.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisXBeamSlice", pSlicePropertiesBeam->m_majorAxis.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisYBeamSlice", pSlicePropertiesBeam->m_majorAxis.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisZBeamSlice", pSlicePropertiesBeam->m_majorAxis.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisEigenvalueBeamSlice", pSlicePropertiesBeam->m_majorAxisEigenvalue));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneXBeamSlice", pSlicePropertiesBeam->m_minorAxisOne.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneYBeamSlice", pSlicePropertiesBeam->m_minorAxisOne.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneZBeamSlice", pSlicePropertiesBeam->m_minorAxisOne.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneEigenvalueBeamSlice", pSlicePropertiesBeam->m_minorAxisOneEigenvalue));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoXBeamSlice", pSlicePropertiesBeam->m_minorAxisTwo.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoYBeamSlice", pSlicePropertiesBeam->m_minorAxisTwo.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoZBeamSlice", pSlicePropertiesBeam->m_minorAxisTwo.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoEigenvalueBeamSlice", pSlicePropertiesBeam->m_minorAxisTwoEigenvalue));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "BestSeparationBeamSlice", pSlicePropertiesBeam->m_bestSeparation));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "AngularSeparationToBeam", pSlicePropertiesBeam->m_angularSeparationToBeam));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClosestDistance", pSlicePropertiesBeam->m_closestDistance));

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DTargetBeamHitsCosmicSlice", pSlicePropertiesCosmic->m_n2DTargetBeamHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DRemnantBeamHitsCosmicSlice", pSlicePropertiesCosmic->m_n2DRemnantBeamHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DCosmicBeamHitsCosmicSlice", pSlicePropertiesCosmic->m_n2DCosmicHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "n2DNoMCParticleHitsCosmicSlice", pSlicePropertiesCosmic->m_n2DNoMCParticleHits));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfosCosmicSlice", pSlicePropertiesCosmic->m_nPfos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nUniqueMCPrimaryParticlesCosmicSlice", pSlicePropertiesCosmic->m_nUniqueMCPrimaryParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroidXCosmicSlice", pSlicePropertiesCosmic->m_centroid.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroidYCosmicSlice", pSlicePropertiesCosmic->m_centroid.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "centroidZCosmicSlice", pSlicePropertiesCosmic->m_centroid.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisXCosmicSlice", pSlicePropertiesCosmic->m_majorAxis.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisYCosmicSlice", pSlicePropertiesCosmic->m_majorAxis.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisZCosmicSlice", pSlicePropertiesCosmic->m_majorAxis.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "majorAxisEigenvalueCosmicSlice", pSlicePropertiesCosmic->m_majorAxisEigenvalue));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneXCosmicSlice", pSlicePropertiesCosmic->m_minorAxisOne.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneYCosmicSlice", pSlicePropertiesCosmic->m_minorAxisOne.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneZCosmicSlice", pSlicePropertiesCosmic->m_minorAxisOne.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisOneEigenvalueCosmicSlice", pSlicePropertiesCosmic->m_minorAxisOneEigenvalue));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoXCosmicSlice", pSlicePropertiesCosmic->m_minorAxisTwo.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoYCosmicSlice", pSlicePropertiesCosmic->m_minorAxisTwo.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoZCosmicSlice", pSlicePropertiesCosmic->m_minorAxisTwo.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minorAxisTwoEigenvalueCosmicSlice", pSlicePropertiesCosmic->m_minorAxisTwoEigenvalue));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "BestSeparationBeamSlice", pSlicePropertiesCosmic->m_bestSeparation));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "AngularSeparationToBeam", pSlicePropertiesCosmic->m_angularSeparationToBeam));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClosestDistance", pSlicePropertiesCosmic->m_closestDistance));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }

//~~~~
    if (m_shouldPerformSliceId)
    {
        for (SliceIdBaseTool *const pSliceIdTool : m_sliceIdToolVector)
            pSliceIdTool->SelectOutputPfos(beamSliceHypotheses, crSliceHypotheses, selectedSlicePfos);
    }
    else if (m_shouldRunNeutrinoRecoOption != m_shouldRunCosmicRecoOption)
    {
        const SliceHypotheses &sliceHypotheses(m_shouldRunNeutrinoRecoOption ? beamSliceHypotheses : crSliceHypotheses);

        for (const PfoList &slice : sliceHypotheses)
            selectedSlicePfos.insert(selectedSlicePfos.end(), slice.begin(), slice.end());
    }

    PfoList newSlicePfoList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(selectedSlicePfos, newSlicePfoList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::GetSliceProperties(const PfoList *const pPfoList, SliceProperties *pSliceProperties)
{
    int n2DTargetBeamHits(0), n2DRemnantBeamHits(0), n2DCosmicHits(0), n2DNoMCParticleHits(0);

    CaloHitList caloHitList;
    PfoList allPfoList;

    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allPfoList);
    pSliceProperties->m_nPfos = allPfoList.size();

    for (const ParticleFlowObject *const pPfo : allPfoList)
    {
        // Get Clusters
        ClusterList clusterList;
        LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);

        // Get Calo Hit List
        for (const Cluster *const pCluster : clusterList)
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        }
    }

    // Make MCParticle To Weight Map and populate with calo hit data
    MCParticleWeightMap mcParticleWeightMap;
    MCParticleVector uniqueMCPrimaryParticles;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CaloHit *const pMasterCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        const MCParticleWeightMap &hitMCParticleWeightMap(pMasterCaloHit->GetMCParticleWeightMap());

        MCParticleVector mcParticleVector;
        for (const MCParticleWeightMap::value_type &mapEntry : hitMCParticleWeightMap) mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        if (mcParticleVector.size() == 0)
        {
            n2DNoMCParticleHits++;
            continue;
        }

        const MCParticle *pDominantMCParticle(nullptr);
        float biggestWeight(0.f);

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(hitMCParticleWeightMap.at(pMCParticle));
            if (weight > biggestWeight) 
            {
                pDominantMCParticle = pMCParticle;
            }

            if (!pDominantMCParticle)
            {
                std::cout << "No dominant particle found." << std::endl;
                continue;
            }
        }

        const MCParticle *pPrimaryMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pDominantMCParticle));
        const LArMCParticle *const pLArPrimaryMCParticle(dynamic_cast<const LArMCParticle*>(pPrimaryMCParticle));

        const float momentum(pLArPrimaryMCParticle->GetMomentum().GetMagnitude());
        const float deltaMomentum(std::fabs(momentum - m_trueMomentum));
        const float fracMomentum(deltaMomentum / m_trueMomentum);

        bool isTargetBeam(pLArPrimaryMCParticle->GetNuanceCode() == 2000 && fracMomentum < (m_momentumAcceptance/100.f));
        bool isRemnantBeam(pLArPrimaryMCParticle->GetNuanceCode() == 2000 && fracMomentum > (m_momentumAcceptance/100.f));
        bool isCR(pLArPrimaryMCParticle->GetNuanceCode() == 3000);

        if (uniqueMCPrimaryParticles.end() != std::find(uniqueMCPrimaryParticles.begin(), uniqueMCPrimaryParticles.end(), pPrimaryMCParticle)) 
            uniqueMCPrimaryParticles.push_back(pPrimaryMCParticle);

        if (isTargetBeam)
        {
            n2DTargetBeamHits++;
        }
        else if (isCR)
        {
            n2DCosmicHits++;
        }
        else if (isRemnantBeam)
        {
            n2DRemnantBeamHits++;
        }
    }

    // 3D
    CaloHitList caloHitList3D;
    for (const ParticleFlowObject *const pPfo : allPfoList)
    {
        // Get Clusters
        ClusterList clusterList;
        LArPfoHelper::GetThreeDClusterList(pPfo, clusterList);

        // Get Calo Hit List
        for (const Cluster *const pCluster : clusterList)
        {
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);
        }
    }

    CaloHitList selectedCaloHitList;
    float closestDistance(std::numeric_limits<float>::max());

    this->GetSelectedCaloHits(caloHitList3D, selectedCaloHitList, closestDistance);

    CartesianVector centroidSel(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecsSel;
    LArPcaHelper::EigenValues eigenValuesSel(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(selectedCaloHitList, centroidSel, eigenValuesSel, eigenVecsSel);
    CartesianVector majorAxisSel(eigenVecsSel.at(0));
    CartesianVector interceptOne(0.f,0.f,0.f), interceptTwo(0.f,0.f,0.f);

    try
    {
        this->GetTPCIntercepts(centroidSel, majorAxisSel, interceptOne, interceptTwo);
    }
    catch (...)
    {
        std::cout << "Unable to find intersection of major axis of slice with TPC" << std::endl;
    }

    // Find the intercept closest to the beam TPC intersection
    float separationOne(std::sqrt(interceptOne.GetDistanceSquared(m_beamTPCIntersection))), separationTwo(std::sqrt(interceptTwo.GetDistanceSquared(m_beamTPCIntersection)));
    float angularSeparationToBeam(majorAxisSel.GetOpeningAngle(m_beamDirection));
/*
    std::cout << "angularSeparationToBeam = " << angularSeparationToBeam << std::endl;
    std::cout << "angularSeparationToBeam deg = " << angularSeparationToBeam*180.f/M_PI << std::endl;
    std::cout << "closestDistance = " << closestDistance << std::endl;
*/
    CartesianVector fitDirection(eigenVecsSel.at(0));

    const int nUniqueMCPrimaryParticles = uniqueMCPrimaryParticles.size();

    int isBeam(n2DTargetBeamHits > (n2DRemnantBeamHits + n2DCosmicHits + n2DNoMCParticleHits) ? 1 : 0);

    pSliceProperties->m_n2DTargetBeamHits = n2DTargetBeamHits;
    pSliceProperties->m_n2DRemnantBeamHits = n2DRemnantBeamHits;
    pSliceProperties->m_n2DCosmicHits = n2DCosmicHits;
    pSliceProperties->m_n2DNoMCParticleHits = n2DNoMCParticleHits;
    pSliceProperties->m_nUniqueMCPrimaryParticles = nUniqueMCPrimaryParticles;
    pSliceProperties->m_isBeam = isBeam;
    pSliceProperties->m_majorAxis = fitDirection;
    pSliceProperties->m_majorAxisEigenvalue = eigenValuesSel.GetX();
    pSliceProperties->m_minorAxisOne = eigenVecsSel.at(1);
    pSliceProperties->m_minorAxisOneEigenvalue = eigenValuesSel.GetY();
    pSliceProperties->m_minorAxisTwo = eigenVecsSel.at(2);
    pSliceProperties->m_minorAxisTwoEigenvalue = eigenValuesSel.GetZ();
    pSliceProperties->m_centroid = centroidSel;
    pSliceProperties->m_bestSeparation = std::min(separationOne, separationTwo);
    pSliceProperties->m_angularSeparationToBeam = angularSeparationToBeam*180.f/M_PI;
    pSliceProperties->m_closestDistance = closestDistance;

//    std::cout << "The number of beam hits is " << n2DBeamHits << "." << std::endl;
//    std::cout << "The number of cosmic ray hits is " << n2DCosmicHits << "." << std::endl;
//    std::cout << "The number of noise hits is " << n2DJunkHits << "." << std::endl;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode MasterAlgorithm::GetTPCIntercepts(CartesianVector &a0, CartesianVector &lineDirection, CartesianVector &interceptOne,  CartesianVector &interceptTwo) const
{
    CartesianVector a(lineDirection.GetUnitVector());
    std::vector<CartesianVector> intercepts;

    for (const Plane *const pPlane : m_tpcPlanes)
    {
        CartesianVector intercept(pPlane->GetLineIntersection(a0, a));
        if (this->IsContained(intercept))
        {
            intercepts.push_back(intercept);
        }
    }

    if (intercepts.size() == 2)
    {
        interceptOne = intercepts.at(0);
        interceptTwo = intercepts.at(1);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MasterAlgorithm::IsContained(CartesianVector &spacePoint) const
{
    // TODO: Make cleaner floating point comparisons 
    float safetyFactor(0.0001);
    bool isContainedX(m_tpcMinX - safetyFactor < spacePoint.GetX() && spacePoint.GetX() < m_tpcMaxX + safetyFactor);
    bool isContainedY(m_tpcMinY - safetyFactor < spacePoint.GetY() && spacePoint.GetY() < m_tpcMaxY + safetyFactor);
    bool isContainedZ(m_tpcMinZ - safetyFactor < spacePoint.GetZ() && spacePoint.GetZ() < m_tpcMaxZ + safetyFactor);

    return isContainedX && isContainedY && isContainedZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode MasterAlgorithm::GetSelectedCaloHits(const CaloHitList &inputCaloHitList, CaloHitList &outputCaloHitList, float &closestHitToFaceDistance) const
{
    if (inputCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    typedef std::pair<const CaloHit*, float> HitDistancePair;
    typedef std::vector<HitDistancePair> HitDistanceVector;
    HitDistanceVector hitDistanceVector;
    
    for (const CaloHit *const pCaloHit : inputCaloHitList)
        hitDistanceVector.emplace_back(pCaloHit, (pCaloHit->GetPositionVector() - m_beamTPCIntersection).GetMagnitudeSquared());

    std::sort(hitDistanceVector.begin(), hitDistanceVector.end(), [](const HitDistancePair &lhs, const HitDistancePair &rhs) -> bool {return (lhs.second < rhs.second);});
    closestHitToFaceDistance = std::sqrt(hitDistanceVector.front().second);

    const unsigned int nInputHits(inputCaloHitList.size());
    const unsigned int nSelectedCaloHits(nInputHits < 100 ? nInputHits : static_cast<unsigned int>(std::round(static_cast<float>(nInputHits) * 10.f / 100.f + 0.5f)));

    for (const HitDistancePair &hitDistancePair : hitDistanceVector)
    {
        outputCaloHitList.push_back(hitDistancePair.first);

        if (outputCaloHitList.size() >= nSelectedCaloHits)
            break;
    }
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Reset()
{
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pCRWorker));

    if (m_pSlicingWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSlicingWorkerInstance));

    if (m_pSliceNuWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceNuWorkerInstance));

    if (m_pSliceCRWorkerInstance)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceCRWorkerInstance));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

MasterAlgorithm::Plane::Plane(CartesianVector &normal, CartesianVector &point) :
    m_unitNormal(0.f, 0.f, 0.f),
    m_point(0.f, 0.f, 0.f)
{
    m_point = point;
    m_unitNormal = normal.GetUnitVector();
    m_a = m_unitNormal.GetX();
    m_b = m_unitNormal.GetY();
    m_c = m_unitNormal.GetZ();
    m_d = -1.f * (normal.GetDotProduct(point));
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector MasterAlgorithm::Plane::GetLineIntersection(CartesianVector &a0, CartesianVector &a) const
{
    if (std::fabs(a.GetDotProduct(m_unitNormal)) < std::numeric_limits<float>::min())
        return CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    float t(-1.f * (a0.GetDotProduct(m_unitNormal) + m_d) / (a.GetDotProduct(m_unitNormal)));
    return CartesianVector(a.GetX() * t + a0.GetX(), a.GetY() * t + a0.GetY(), a.GetZ() * t + a0.GetZ());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const CaloHit *const pCaloHit) const
{
    // TODO Useful protection to ensure that only calo hits owned by the master instance are copied
    PandoraApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pCaloHit->GetCellSize1();
    parameters.m_cellThickness = pCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit->GetTime();
    parameters.m_inputEnergy = pCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit->IsDigital();
    parameters.m_hitType = pCaloHit->GetHitType();
    parameters.m_hitRegion = pCaloHit->GetHitRegion();
    parameters.m_layer = pCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void*>(pCaloHit);
    return PandoraApi::CaloHit::Create(*pPandora, parameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const PfoList &inputPfoList, PfoList &newPfoList) const
{
    if (inputPfoList.empty())
        return STATUS_CODE_SUCCESS;

    // TODO if no pfo in input list is primary - raise exception

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Pfo *const pPfo : inputPfoList)
    {
        if (pPfo->GetParentPfoList().empty())
            this->Recreate(pPfo, nullptr, newPfoList);
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_recreatedClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_recreatedClusterListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_recreatedVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_recreatedVertexListName));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, m_recreatedPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, m_recreatedPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *MasterAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const CaloHit *const pParentCaloHit) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void*>(pParentCaloHit);

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *MasterAlgorithm::CreateCluster(const Cluster *const pInputCluster, const CaloHitList &newCaloHitList,
    const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *MasterAlgorithm::CreateVertex(const Vertex *const pInputVertex) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition();
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *MasterAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ClusterList &newClusterList,
    const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPC &larTPC, const DetectorGapList &gapList, const std::string &settingsFile) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The LArTPC
    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = larTPC.GetLArTPCVolumeId();
    larTPCParameters.m_centerX = larTPC.GetCenterX();
    larTPCParameters.m_centerY = larTPC.GetCenterY();
    larTPCParameters.m_centerZ = larTPC.GetCenterZ();
    larTPCParameters.m_widthX = larTPC.GetWidthX();
    larTPCParameters.m_widthY = larTPC.GetWidthY();
    larTPCParameters.m_widthZ = larTPC.GetWidthZ();
    larTPCParameters.m_wirePitchU = larTPC.GetWirePitchU();
    larTPCParameters.m_wirePitchV = larTPC.GetWirePitchV();
    larTPCParameters.m_wirePitchW = larTPC.GetWirePitchW();
    larTPCParameters.m_wireAngleU = larTPC.GetWireAngleU();
    larTPCParameters.m_wireAngleV = larTPC.GetWireAngleV();
    larTPCParameters.m_sigmaUVW = larTPC.GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = larTPC.IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    const float tpcMinX(larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX()), tpcMaxX(larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX());

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap && (((pLineGap->GetLineEndX() >= tpcMinX) && (pLineGap->GetLineEndX() <= tpcMaxX)) ||
            ((pLineGap->GetLineStartX() >= tpcMinX) && (pLineGap->GetLineStartX() <= tpcMaxX))) )
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            const LineGapType lineGapType(pLineGap->GetLineGapType());
            lineGapParameters.m_lineGapType = lineGapType;
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();

            if (m_fullWidthCRWorkerWireGaps && ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W)))
            {
                lineGapParameters.m_lineStartX = -std::numeric_limits<float>::max();
                lineGapParameters.m_lineEndX = std::numeric_limits<float>::max();
            }

            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPCMap &larTPCMap, const DetectorGapList &gapList, const std::string &settingsFile) const
{
    if (larTPCMap.empty())
    {
        std::cout << "MasterAlgorithm::CreateWorkerInstance - no LArTPC details provided" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    // The Pandora instance
    const Pandora *const pPandora(new Pandora());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The Parent LArTPC
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCVolumeId = 0;
    larTPCParameters.m_centerX = 0.5f * (parentMaxX + parentMinX);
    larTPCParameters.m_centerY = 0.5f * (parentMaxY + parentMinY);
    larTPCParameters.m_centerZ = 0.5f * (parentMaxZ + parentMinZ);
    larTPCParameters.m_widthX = parentMaxX - parentMinX;
    larTPCParameters.m_widthY = parentMaxY - parentMinY;
    larTPCParameters.m_widthZ = parentMaxZ - parentMinZ;
    larTPCParameters.m_wirePitchU = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchV = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchW = pFirstLArTPC->GetWirePitchW();
    larTPCParameters.m_wireAngleU = pFirstLArTPC->GetWireAngleU();
    larTPCParameters.m_wireAngleV = pFirstLArTPC->GetWireAngleV();
    larTPCParameters.m_sigmaUVW = pFirstLArTPC->GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = pFirstLArTPC->IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap)
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            lineGapParameters.m_lineGapType = pLineGap->GetLineGapType();
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();
            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    ExternalSteeringParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        pExternalParameters = dynamic_cast<ExternalSteeringParameters*>(this->GetExternalParameters());
        if (!pExternalParameters) return STATUS_CODE_FAILURE;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunAllHitsCosmicReco, xmlHandle, "ShouldRunAllHitsCosmicReco", m_shouldRunAllHitsCosmicReco));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunStitching, xmlHandle, "ShouldRunStitching", m_shouldRunStitching));

    if (m_shouldRunStitching && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunStitching requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunStitching)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "StitchingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            StitchingBaseTool *const pStitchingTool(dynamic_cast<StitchingBaseTool*>(pAlgorithmTool));
            if (!pStitchingTool) return STATUS_CODE_INVALID_PARAMETER;
            m_stitchingToolVector.push_back(pStitchingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicHitRemoval, xmlHandle, "ShouldRunCosmicHitRemoval", m_shouldRunCosmicHitRemoval));

    if (m_shouldRunCosmicHitRemoval && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunCosmicHitRemoval requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "CosmicRayTaggingTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            CosmicRayTaggingBaseTool *const pCosmicRayTaggingTool(dynamic_cast<CosmicRayTaggingBaseTool*>(pAlgorithmTool));
            if (!pCosmicRayTaggingTool) return STATUS_CODE_INVALID_PARAMETER;
            m_cosmicRayTaggingToolVector.push_back(pCosmicRayTaggingTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunSlicing, xmlHandle, "ShouldRunSlicing", m_shouldRunSlicing));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunNeutrinoRecoOption, xmlHandle, "ShouldRunNeutrinoRecoOption", m_shouldRunNeutrinoRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicRecoOption, xmlHandle, "ShouldRunCosmicRecoOption", m_shouldRunCosmicRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldPerformSliceId, xmlHandle, "ShouldPerformSliceId", m_shouldPerformSliceId));

    if (m_shouldPerformSliceId && (!m_shouldRunSlicing || !m_shouldRunNeutrinoRecoOption || !m_shouldRunCosmicRecoOption))
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldPerformSliceId requires ShouldRunSlicing and both neutrino and cosmic reconstruction options" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldPerformSliceId)
    {
        AlgorithmToolVector algorithmToolVector;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "SliceIdTools", algorithmToolVector));

        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        {
            SliceIdBaseTool *const pSliceIdIdTool(dynamic_cast<SliceIdBaseTool*>(pAlgorithmTool));
            if (!pSliceIdIdTool) return STATUS_CODE_INVALID_PARAMETER;
            m_sliceIdToolVector.push_back(pSliceIdIdTool);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_printOverallRecoStatus, xmlHandle, "PrintOverallRecoStatus", m_printOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualizeOverallRecoStatus", m_visualizeOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullWidthCRWorkerWireGaps", m_fullWidthCRWorkerWireGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CRSettingsFile", m_crSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuSettingsFile", m_nuSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SlicingSettingsFile", m_slicingSettingsFile));
    m_crSettingsFile = LArFileHelper::FindFileInPath(m_crSettingsFile, m_filePathEnvironmentVariable);
    m_nuSettingsFile = LArFileHelper::FindFileInPath(m_nuSettingsFile, m_filePathEnvironmentVariable);
    m_slicingSettingsFile = LArFileHelper::FindFileInPath(m_slicingSettingsFile, m_filePathEnvironmentVariable);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputHitListName", m_inputHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedPfoListName", m_recreatedPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedClusterListName", m_recreatedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RecreatedVertexListName", m_recreatedVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SliceAnalysisFileName", m_fileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrueMomentum", m_trueMomentum));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MomentumAcceptance", m_momentumAcceptance));

    FloatVector beamTPCIntersection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamTPCIntersection", beamTPCIntersection));

    if (3 == beamTPCIntersection.size())
    {
        m_beamTPCIntersection.SetValues(beamTPCIntersection.at(0), beamTPCIntersection.at(1), beamTPCIntersection.at(2));
    }
    else
    {
        // Default for protoDUNE.
        m_beamTPCIntersection.SetValues(-33.051, 461.06, 0);
    }

    FloatVector beamDirection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamDirection", beamDirection));

    if (3 == beamDirection.size())
    {
        m_beamDirection.SetValues(beamDirection.at(0), beamDirection.at(1), beamDirection.at(2));
    }
    else
    {
        // Default for protoDUNE.
        float thetaXZ0(-11.844f * M_PI / 180.f);
        m_beamDirection.SetValues(std::sin(thetaXZ0), 0, std::cos(thetaXZ0));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const InputBool inputBool,
    const TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool)
{
    if (pExternalParameters && inputBool.IsInitialized())
    {
        outputBool = inputBool.Get();
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, xmlTag, outputBool));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

MasterAlgorithm::SliceProperties::SliceProperties() : 
    m_n2DTargetBeamHits(std::numeric_limits<int>::max()),
    m_n2DRemnantBeamHits(std::numeric_limits<int>::max()),
    m_n2DCosmicHits(std::numeric_limits<int>::max()),
    m_n2DNoMCParticleHits(std::numeric_limits<int>::max()),
    m_nPfos(std::numeric_limits<int>::max()),
    m_nUniqueMCPrimaryParticles(std::numeric_limits<int>::max()),
    m_isBeam(std::numeric_limits<int>::max()),
    m_centroid(CartesianVector(0.f, 0.f, 0.f)),
    m_majorAxis(CartesianVector(0.f, 0.f, 0.f)),
    m_majorAxisEigenvalue(std::numeric_limits<float>::max()),
    m_minorAxisOne(CartesianVector(0.f, 0.f, 0.f)),
    m_minorAxisOneEigenvalue(std::numeric_limits<float>::max()),
    m_minorAxisTwo(CartesianVector(0.f, 0.f, 0.f)),
    m_minorAxisTwoEigenvalue(std::numeric_limits<float>::max()),
    m_bestSeparation(std::numeric_limits<float>::max()),
    m_angularSeparationToBeam(std::numeric_limits<float>::max()),
    m_closestDistance(std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MasterAlgorithm::SliceProperties::Reset()
{
    m_n2DTargetBeamHits = std::numeric_limits<int>::max();
    m_n2DRemnantBeamHits = std::numeric_limits<int>::max();
    m_n2DCosmicHits = std::numeric_limits<int>::max();
    m_n2DNoMCParticleHits = std::numeric_limits<int>::max();
    m_nPfos = std::numeric_limits<int>::max();
    m_nUniqueMCPrimaryParticles = std::numeric_limits<int>::max();
    m_isBeam = std::numeric_limits<int>::max();

    m_centroid.SetValues(0.f, 0.f, 0.f);
    m_majorAxis.SetValues(0.f, 0.f, 0.f);
    m_majorAxisEigenvalue = std::numeric_limits<float>::max();
    m_minorAxisOne.SetValues(0.f, 0.f, 0.f);
    m_minorAxisOneEigenvalue = std::numeric_limits<float>::max();
    m_minorAxisTwo.SetValues(0.f, 0.f, 0.f);
    m_minorAxisTwoEigenvalue = std::numeric_limits<float>::max();

    m_bestSeparation = std::numeric_limits<float>::max();
    m_angularSeparationToBeam = std::numeric_limits<float>::max();
    m_closestDistance = std::numeric_limits<float>::max();
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
