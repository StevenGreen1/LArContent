/**
 *  @file   larpandoracontent/LArVertex/BdtVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the bdt vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_BDT_VERTEX_SELECTION_ALGORITHM_H
#define LAR_BDT_VERTEX_SELECTION_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"

#include "larpandoracontent/LArVertex/MLVertexSelectionBaseAlgorithm.h"

#include <random>

namespace lar_content
{

/**
 *  @brief  BdtVertexSelectionAlgorithm class
 */
class BdtVertexSelectionAlgorithm : public MLVertexSelectionBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    BdtVertexSelectionAlgorithm();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief  Used a binary classifier to compare a set of vertices and pick the best one
     *
     *  @param  vertexVector the vector of vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  eventFeatureList the event feature list
     *  @param  useRPhi whether to include the r/phi feature
     *
     *  @return address of the best vertex
     */
    const pandora::Vertex * CompareVertices(const pandora::VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        const LArMvaHelper::MvaFeatureVector &eventFeatureList, const bool useRPhi, const bool isRegion = false) const;

    std::string           m_filePathEnvironmentVariable;          ///< The environment variable providing a list of paths to bdt files
    std::string           m_bdtFileName;                          ///< The Bdt file name
    std::string           m_regionBdtName;                        ///< The name of the region bdt to find
    std::string           m_vertexBdtName;                        ///< The name of the vertex bdt to find
    AdaBoostDecisionTree  m_bdtRegion;                            ///< The region boosted decision tree
    AdaBoostDecisionTree  m_bdtVertex;                            ///< The vertex boosted decision tree
};

} // namespace lar_content

#endif // #ifndef LAR_BDT_VERTEX_SELECTION_ALGORITHM_H
