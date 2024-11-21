//
//  SampledCladogenesisRootFrequenciesFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 8/11/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <iterator>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "SampledCladogenesisRootFrequenciesFunction.h"
#include "RbException.h"
#include "TransitionProbabilityMatrix.h"
#include "BranchHistory.h"
#include "BranchHistoryDiscrete.h"
#include "CharacterEvent.h"
#include "CharacterEventCompare.h"
#include "CharacterHistoryDiscrete.h"
#include "Distribution.h"
#include "RateGenerator.h"
#include "RbVectorImpl.h"
#include "Simplex.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class MatrixReal; }

using namespace RevBayesCore;

SampledCladogenesisRootFrequenciesFunction::SampledCladogenesisRootFrequenciesFunction(const TypedDagNode<RateGenerator> *rm,
                                                                   const TypedDagNode<MatrixReal> *cp,
                                                                   const TypedDagNode<Tree> *t,
                                                                   const TypedDagNode<double> *r)
: TypedFunction<Simplex>(new Simplex(rm->getValue().getNumberOfStates(), (1.0/rm->getValue().getNumberOfStates()))),
    rateMatrix(rm),
    cladogenesisProbabilities(cp),
    tree(t),
    clock(r),
    numStates(rm->getValue().getNumberOfStates()),
    initialLengthForStationarity(30.0),
    tunedLengthForStationarity(0.0),
    absorptionTolerance(1.0 - 1e-10),
    stationarityDifferenceTolerance(1e-6),
    stationarityPositiveTolerance(1e-6),
    isSampledSpeciationTree(false)
{
    
    // add the rate and frequency parameters as parents
    addParameter( rateMatrix );
    addParameter( cladogenesisProbabilities );
    addParameter( tree );
    addParameter( clock );
    
    tunedLengthForStationarity = initialLengthForStationarity;
    isSampledSpeciationTree = dynamic_cast<const RevBayesCore::AbstractCharacterHistoryBirthDeathProcess* >( &(tree->getDistribution()) ) != NULL;
    
    update();
}


SampledCladogenesisRootFrequenciesFunction::~SampledCladogenesisRootFrequenciesFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



SampledCladogenesisRootFrequenciesFunction* SampledCladogenesisRootFrequenciesFunction::clone( void ) const {
    return new SampledCladogenesisRootFrequenciesFunction( *this );
}


void SampledCladogenesisRootFrequenciesFunction::update( void ) {
    
    // get anagenesis clock
    double rate = clock->getValue();
    
    // get anagenetic rate matrix
    const RateGenerator& rm = rateMatrix->getValue();

    // get root
    const TopologyNode& root = tree->getValue().getRoot();
    
    // compute root freqs
    TransitionProbabilityMatrix tp(numStates);
    TransitionProbabilityMatrix tp_next(numStates);
    
    // this many anagenetic events are expected to occur BEFORE cladogenetic events are applied
    rm.calculateTransitionProbabilities( tunedLengthForStationarity, 0, 1.0, tp );
    
    
//    std::cout << tp << "\n";
    // tune this amount if needed
    bool tune = true;
    bool failed = false;
    size_t tuningAttempts = 0;
    while (tune ) {
        
        tune = false;

        double tp_10 = tp[1][0];
        double tp_11 = tp[1][1];
        double tp_k1 = tp[numStates-1][1];
        double tp_1k = tp[1][numStates-1];
        
        // if the branch length is too std::int64_t, the tp_ij for j!=0 probabilities are too small
        if (tp_10 > absorptionTolerance ||
            tp_1k < 0.0 ||
            tp_k1 < 0.0)
        {
            tunedLengthForStationarity *= 0.80;
            rm.calculateTransitionProbabilities( tunedLengthForStationarity, 0, 1.0, tp );
            tune = true;
            //            std::cout << "Absorption error term: " << tp[1][1] - tp[numStates-1][1] << "\n";
            //            std::cout << "tunedLengthForStationarity: " << tunedLengthForStationarity << "\n";
        }
        
        // if the expected length is too short, then we have not reached stationarity
        else if (fabs(tp_11 - tp_k1) > stationarityDifferenceTolerance ||
            tp_1k < stationarityPositiveTolerance ||
            tp_k1 < stationarityPositiveTolerance )
        {
            tunedLengthForStationarity *= 1.33;
            rm.calculateTransitionProbabilities( tunedLengthForStationarity, 0, 1.0, tp );
            tune = true;
            //            std::cout << "Stationarity error term: " << tp[1][1] - tp[numStates-1][1] << "\n";
            //            std::cout << "tunedLengthForStationarity: " << tunedLengthForStationarity << "\n";

        }
        
        if (tuningAttempts > 100) {
             RbException("SampledCladogenesisRootFrequenciesFunction could not satisfy absorption/stationarity tolerances");
            tune = false;
            failed = true;
        }
        
        
//        if (tuningAttempts > 0)
//        {
//            std::cout << tuningAttempts << " " << tunedLengthForStationarity << " " << 1.0-tp_10 << " " << tp_11-tp_k1 << " " << tp_1k << " " << tp_k1 << "\n";
//        }
        
        tuningAttempts++;
    }
   
    // account for hidden speciation events only if it is a sampled speciation tree
    if (isSampledSpeciationTree)
    {
        // get cladogenetic probabilities
        //cladogenesisProbabilities->getValue();
        //const DeterministicNode<MatrixReal>* cpn = static_cast<const DeterministicNode<MatrixReal>* >( cladogenesisProbabilities );
        //const TypedFunction<MatrixReal>& tf = cpn->getFunction();
        //const AbstractCladogenicStateFunction* csf = dynamic_cast<const AbstractCladogenicStateFunction*>( &tf );
        
        // TODO: MJL fix this later, use new CladogeneticProbabilityMatrix setup
        std::map<std::vector<unsigned>, double> eventMapProbs;
        // eventMapProbs = csf->getEventMap();
        
        // get root history information
        const AbstractCharacterHistoryBirthDeathProcess* dist = dynamic_cast<const AbstractCharacterHistoryBirthDeathProcess* >( &tree->getDistribution() );
        const CharacterHistoryDiscrete& tree_history = static_cast<const CharacterHistoryDiscrete&>(dist->getCharacterHistory());
        const BranchHistory& branch_history = tree_history[root.getIndex()];
        const std::multiset<CharacterEvent*,CharacterEventCompare>& events = branch_history.getHistory();
        
        // for each interval between events, go from present to past
        double brlen = root.getAge();
        double startAge = root.getAge() + brlen;
        double endAge = root.getAge();
        
        double t = 0.0;
        double dt = 0.0;
        double prev_age = startAge;
        double event_age = prev_age;
        std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it;
        for (it = events.rbegin(); it != events.rend(); it++)
        {
            
            prev_age = event_age;
            event_age = (*it)->getAge();
            
            // anagenetic changes occurring between (event_age, event_age-dt)
            
            rm.calculateTransitionProbabilities(prev_age, event_age, rate, tp_next );
            tp *= tp_next;
            
            // first compute clado probs at younger end of branch
            TransitionProbabilityMatrix cp(numStates);
            std::map<std::vector<unsigned>, double>::iterator it;
            for ( it = eventMapProbs.begin(); it != eventMapProbs.end(); ++it)
            {
                // sparse elements from map
                const std::vector<unsigned>& idx = it->first;
                const size_t j = idx[0];
                const size_t k = idx[1];
                const size_t l = idx[2];
                double p_clado = 0.5 * it->second;
                
                cp[j][k] += p_clado;
                cp[j][l] += p_clado;
                
            }
            
            t += dt;
            tp *= cp;
        }
        
        // last interval
        dt = brlen - t;
        
        rm.calculateTransitionProbabilities( event_age, endAge, rate, tp_next );
        tp *= tp_next;
        

    }
    
    // get sum of stationary probabilities
    double sum = 0.0;
    for (size_t i = 0; i < value->size(); i++)
    {
        sum += tp[1][i];
    }
    
    // make sure sum is within certain tolerances
    if ( fabs(1.0 - sum) > 1e-4 || sum < 0.0 || std::isnan(sum) )
    {
        failed = true;
//        std::cout << tp << "\n";
//        RbException("SampledCladogenesisRootFrequenciesFunction simplex does not sum to one");
        
    }
   
    // force 0.0 probability if the tuning heuristic failed
    if (failed) {
        for (size_t i = 0; i < value->size(); i++)
        {
            (*value)[i] = 0.0;
        }
        // reset the tuned branch length
        tunedLengthForStationarity = initialLengthForStationarity;
    }
    else {
        for (size_t i = 0; i < value->size(); i++)
        {
            (*value)[i] = tp[1][i] / sum;
//            if (std::isnan((*value)[i])) {
////                std::cout << tp << "\n";
//                RbException("SampledCladogenesisRootFrequenciesFunction contains nan!");
//            }
//    //        (*value)[i] = (*value)[i] / sum;
        }
    }

}



void SampledCladogenesisRootFrequenciesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    if (oldP == rateMatrix) {
        rateMatrix = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    }
    else if (oldP == cladogenesisProbabilities) {
        cladogenesisProbabilities = static_cast<const TypedDagNode< MatrixReal >* >( newP );
    }
    else if (oldP == tree) {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == clock) {
        clock = static_cast<const TypedDagNode<double>* >( newP );
    }
}
