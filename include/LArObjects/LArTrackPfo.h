/**
 *  @file   LArContent/include/LArObjects/LArTrackPfo.h
 *
 *  @brief  Header file for the lar pfo class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_PFO_H
#define LAR_TRACK_PFO_H 1

#include "Objects/TrackState.h"
#include "Objects/CartesianVector.h"
#include "Objects/ParticleFlowObject.h"

#include "Pandora/ObjectFactory.h"

namespace lar_content
{

class LArTrackState : public pandora::TrackState
{
public:    
    /**
     *  @brief  Constructor
     *
     *  @param  position
     *  @param  direction
     *  @param  hitType
     *  @param  dQ
     *  @param  dL
     */
    LArTrackState(const pandora::CartesianVector &position, const pandora::CartesianVector &direction, const pandora::HitType hitType,
        const float dQ, const float dL);

    /**
     *  @brief  Return direction at this trajectory point
     */
    const pandora::CartesianVector &GetDirection() const;

    /**
     *  @brief Return hit type of this trajectory point
     */
    pandora::HitType GetHitType() const;

    /**
     *  @brief  Return dQ at this trajectory point
     */
    float GetdQ() const;
  
    /**
     *  @brief  Return dL at this trajectory point
     */
    float GetdL() const;
  
    /**
     *  @brief Return dQ/dL at this trajectory point 
     */
    float GetdQdL() const;

private:
    pandora::HitType   m_hitType;
    float              m_dQ;
    float              m_dL;
};

typedef std::vector<LArTrackState> LArTrackStateVector; 
typedef std::pair<float, LArTrackState> LArTrackTrajectoryPoint;
typedef std::vector<LArTrackTrajectoryPoint> LArTrackTrajectory;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo parameters
 */
class LArTrackPfoParameters : public PandoraContentApi::ParticleFlowObject::Parameters
{
public:
    LArTrackStateVector   m_trackStateVector;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object
 */
class LArTrackPfo : public pandora::ParticleFlowObject
{
public:
    /**
     *  @brief  Constructor
     */
    LArTrackPfo(const LArTrackPfoParameters &parameters);

    /**
     *  @brief  Get vertex position
     */
    const pandora::CartesianVector &GetVertexPosition() const;
  
    /**
     *  @brief  Get end position
     */
    const pandora::CartesianVector &GetEndPosition() const;
 
    /**
     *  @brief  Get vertex direction 
     */
    const pandora::CartesianVector &GetVertexDirection() const;

    /**
     *  @brief  Get end direction
     */
    const pandora::CartesianVector &GetEndDirection() const;

public:
    const LArTrackStateVector   m_trackStateVector;      ///< The vector of track states

private:
    // OTHER MEMBER VARIABLES GO HERE
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  lar pfo object factory responsible for pfo creation
 */
class LArTrackPfoFactory : public pandora::ObjectFactory<PandoraContentApi::ParticleFlowObject::Parameters, pandora::ParticleFlowObject>
{
public:
    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  fileReader the file reader, used to extract any additional parameters from file
     */
    pandora::StatusCode Read(Parameters &parameters, pandora::FileReader &fileReader) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fileWriter the file writer
     */
    pandora::StatusCode Write(const pandora::ParticleFlowObject *const pObject, pandora::FileWriter &fileWriter) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const PandoraContentApi::ParticleFlowObject::Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArTrackPfoFactory::Parameters *LArTrackPfoFactory::NewParameters() const
{
    return (new LArTrackPfoParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArTrackPfoFactory::Create(const Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const
{
    const LArTrackPfoParameters &larPfoParameters(dynamic_cast<const LArTrackPfoParameters&>(parameters));
    pObject = new LArTrackPfo(larPfoParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArTrackPfoFactory::Read(Parameters&, pandora::FileReader&) const
{
    // TODO: Provide this functionality when necessary

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArTrackPfoFactory::Write(const pandora::ParticleFlowObject*, pandora::FileWriter&) const
{
    // TODO: Provide this functionality when necessary

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_PFO_H
