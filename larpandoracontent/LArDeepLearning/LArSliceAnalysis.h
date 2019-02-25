/**
 *  @file   larpandoracontent/LArDeepLearning/LArSliceAnalysis.h
 *
 *  @brief  Header file for the lar slice analysis class.
 *
 *  $Log: $
 */

#ifndef LAR_SLICE_ANALYSIS_H
#define LAR_SLICE_ANALYSIS_H 1

#include "Objects/CaloHit.h"
#include "Objects/MCParticle.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/StatusCodes.h"

#include <map>
#include <vector>

namespace lar_content
{

typedef std::vector<pandora::CaloHitList> SliceVector;

/**
 *  @brief  SliceAnalysis class
 */
class SliceAnalysis
{
public:
    /**
     *  @brief  Constructor
     */
    SliceAnalysis();

    /**
     *  @brief  Constructor
     *
     *  @param  trainingOutputFile name of input training file
     */
    SliceAnalysis(const std::string &trainingOutputFile);

    /**
     *  @brief  Destructor
     */
    ~SliceAnalysis();

    /**
     *  @brief  Process slices
     *
     *  @param  sliceVector slices to process
     *  @param  eventNumber
     */
    void ProcessSlice(const SliceVector &sliceVector, const int eventNumber);

private:

    /**
     *  @brief  Write hit list to text file
     *
     *  @param  eventNumber event number
     *  @param  sliceCounter slice index
     *  @param  hitType of calo hits being written out
     *  @param  caloHitList to write out
     */
    void WriteHitList(const int eventNumber, const int sliceCounter, const pandora::HitType &hitType,
         const pandora::CaloHitList &caloHitList);

    const bool        m_useTrainingMode;     ///< Should use training mode. If true, training examples will be written to the output file
    const std::string m_trainingOutputFile;  ///< Output file name for training examples, gets appended with view type
    const int         m_minTargetHits;       ///< Min number of hits to be classified as target MCParticle
};

} // namespace lar_content

#endif // #ifndef LAR_SLICE_ANALYSIS_H
