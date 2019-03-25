/**
 *  @file   larpandoracontent/LArDeepLearning/LArSliceAnalysis.cc
 *
 *  @brief  Implementation of the lar slice analysis class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArDeepLearning/LArSliceAnalysis.h"
//#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

using namespace pandora;

namespace lar_content
{

SliceAnalysis::SliceAnalysis() :
    m_useTrainingMode(false),
    m_trainingOutputFile(""),
    m_minTargetHits(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceAnalysis::SliceAnalysis(const std::string &trainingOutputFile) :
    m_useTrainingMode(true),
    m_trainingOutputFile(trainingOutputFile),
    m_minTargetHits(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceAnalysis::~SliceAnalysis()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceAnalysis::ProcessSlice(const SliceVector &sliceVector, const int eventNumber)
{
    if (!m_useTrainingMode)
        return;

    int sliceCounter(0);

    for (const CaloHitList &slice : sliceVector)
    {
        sliceCounter++;
        CaloHitList caloHitListU, caloHitListV, caloHitListW, caloHitList3D;

        for (const CaloHit *pCaloHit : slice)
        {
            if (pCaloHit->GetHitType() == TPC_VIEW_U)
            {
                caloHitListU.push_back(pCaloHit);
            }
            else if (pCaloHit->GetHitType() == TPC_VIEW_V)
            {
                caloHitListV.push_back(pCaloHit);
            }
            else if (pCaloHit->GetHitType() == TPC_VIEW_W)
            {
                caloHitListW.push_back(pCaloHit);
            }
            else if (pCaloHit->GetHitType() == TPC_3D)
            {
                caloHitList3D.push_back(pCaloHit);
            }
        }

        this->WriteHitList(eventNumber, sliceCounter, TPC_VIEW_U, caloHitListU);
        this->WriteHitList(eventNumber, sliceCounter, TPC_VIEW_V, caloHitListV);
        this->WriteHitList(eventNumber, sliceCounter, TPC_VIEW_W, caloHitListW);
        this->WriteHitList(eventNumber, sliceCounter, TPC_3D, caloHitList3D);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceAnalysis::WriteHitList(const int eventNumber, const int sliceCounter, const HitType &hitType, const CaloHitList &caloHitList)
{
    LArMvaHelper::MvaFeatureVector featureVector;

    std::string trainingOutputFileName(m_trainingOutputFile);

    if (hitType == TPC_3D)
    {
        trainingOutputFileName += "_3D.txt";
    }
    else if (hitType == TPC_VIEW_U)
    {
        trainingOutputFileName += "_UView.txt";
    }
    else if (hitType == TPC_VIEW_V)
    {
        trainingOutputFileName += "_VView.txt";
    }
    else if (hitType == TPC_VIEW_W)
    {
        trainingOutputFileName += "_WView.txt";
    }
    else
    {
        std::cout << "SliceAnalysis::WriteHitList - Invalid hit type" << std::endl;
        return;
    }

    // EventNumber, SliceNumber, nCaloHits
    featureVector.push_back(static_cast<double>(eventNumber));
    featureVector.push_back(static_cast<double>(sliceCounter));
    featureVector.push_back(static_cast<double>(caloHitList.size()));

    std::map<const MCParticle*, int> mcParticlesToNHitsMap;
    std::map<const MCParticle*, int> mcParticlesToIndex;

    int index(0), globalIndex(0);

    for (const CaloHit *pCaloHit : caloHitList)
    {
        const CaloHit *pCaloHitInMaster(static_cast<const CaloHit*>(pCaloHit->GetParentAddress()));
        int pdg(-1);

        if (hitType == TPC_3D)
            pCaloHitInMaster = static_cast<const CaloHit*>(pCaloHitInMaster->GetParentAddress());

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHitInMaster));
            pdg = pMCParticle->GetParticleId();

            // ATTN: All MCParticles, not just primaries
            if (mcParticlesToNHitsMap.find(pMCParticle) == mcParticlesToNHitsMap.end())
            {
                mcParticlesToNHitsMap.insert(std::make_pair(pMCParticle, 1));
            }
            else
            {
                mcParticlesToNHitsMap.at(pMCParticle)++;
            }

            if (mcParticlesToIndex.find(pMCParticle) == mcParticlesToIndex.end())
            {
                mcParticlesToIndex.insert(std::make_pair(pMCParticle, globalIndex));
                globalIndex++;
            }

            index = mcParticlesToIndex.at(pMCParticle);
        }
        catch(...)
        {
            // ATTN: Only put hits in with an MCParticle
            continue;
        }

        featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetX()));
        featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetY()));
        featureVector.push_back(static_cast<double>(pCaloHit->GetPositionVector().GetZ()));
        featureVector.push_back(static_cast<double>(pdg));
        featureVector.push_back(static_cast<double>(index));
    }

    int nTargetTrks(0), nTargetShws(0);

    for (const auto iter : mcParticlesToNHitsMap)
    {
        if (iter.second < m_minTargetHits)
            continue;

        if (std::fabs(iter.first->GetParticleId()) == 11 || iter.first->GetParticleId() == 22)
        {
            nTargetShws++;
        }
        else
        {
            nTargetTrks++;
        }
    }

    // ATTN: Only put in MC particles with at least a few hits
    featureVector.push_back(static_cast<double>(nTargetTrks));
    featureVector.push_back(static_cast<double>(nTargetShws));

    // ATTN: Bool here is redundant
    LArMvaHelper::ProduceTrainingExample(trainingOutputFileName, true, featureVector);
}

} // namespace lar_content
