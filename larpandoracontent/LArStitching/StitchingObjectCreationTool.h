/**
 *  @file   larpandoracontent/LArStitching/StitchingObjectCreationTool.h
 *
 *  @brief  Header file for the stitching object creation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_OBJECT_CREATION_TOOL_H
#define LAR_STITCHING_OBJECT_CREATION_TOOL_H 1

#include "Geometry/LArTPC.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  StitchingObjectCreationTool class
 */
class StitchingObjectCreationTool : public StitchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    StitchingObjectCreationTool();

    typedef StitchingAlgorithm::StitchingInfo StitchingInfo;

    void Run(const StitchingAlgorithm *const pAlgorithm, StitchingInfo &stitchingInfo);

private:
    /**
     *  @brief  Recreate all (3D aspects of) pfos associated with a provided pandora instance
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pPandora the address of the input pandora instance acting as a source for pfos
     *  @param  larTPC the lar tpc description association with the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const pandora::Pandora *const pPandora, const pandora::LArTPC &larTPC,
        StitchingInfo &stitchingInfo) const;

    /**
     *  @brief  Recreate (3D aspects of) a provided pfo
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pInputPfo the address of the input pfo
     *  @param  pNewParentPfo the address of the parent for the new pfo, nullptr if not relevant
     *  @param  larTPC the lar tpc description association with the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ParticleFlowObject *const pNewParentPfo, const pandora::LArTPC &larTPC, StitchingInfo &stitchingInfo) const;

    /**
     *  @brief  Add details to the stitching information block
     *
     *  @param  pNewPfo the address of a new pfo, recreated from an input pfo
     *  @param  larTPC the lar tpc description association with the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void AddStitchingInfo(const pandora::ParticleFlowObject *const pNewPfo, const pandora::LArTPC &larTPC, StitchingInfo &stitchingInfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_OBJECT_CREATION_TOOL_H
