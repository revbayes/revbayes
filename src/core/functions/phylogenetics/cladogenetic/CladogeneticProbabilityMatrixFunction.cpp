//
//  CladogeneticProbabilityMatrixFunction.cpp
//
//  Created by Will Freyman on 8/1/17.
//


#include "CladogeneticProbabilityMatrixFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

#include <cstddef>

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


//TypedFunction<MatrixReal>( new MatrixReal( mc + 1, (mc + 1) * (mc + 1), 0.0 ) ),
CladogeneticProbabilityMatrixFunction::CladogeneticProbabilityMatrixFunction(const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<long> > >* events, const TypedDagNode<RevBayesCore::RbVector<double> >* probs, int n_states):
TypedFunction<CladogeneticProbabilityMatrix>( new CladogeneticProbabilityMatrix( n_states ) ),
cladogenetic_events( events ),
num_states( n_states ),
probabilities( probs )
{
    addParameter( probabilities );
    
    // since the transition rate matrix will be very large (num_states by num_states^2) but
    // only sparsely filled, the we use an event map instead
    // of the full matrix for efficiency
    update();
}


CladogeneticProbabilityMatrixFunction::~CladogeneticProbabilityMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


CladogeneticProbabilityMatrixFunction* CladogeneticProbabilityMatrixFunction::clone( void ) const
{
    return new CladogeneticProbabilityMatrixFunction( *this );
}


std::map< std::vector<unsigned>, double >  CladogeneticProbabilityMatrixFunction::getEventMap(double t)
{
    return event_map;
}

const std::map< std::vector<unsigned>, double >&  CladogeneticProbabilityMatrixFunction::getEventMap(double t) const
{
    return event_map;
}


void CladogeneticProbabilityMatrixFunction::update( void )
{
    // reset the transition matrix
    delete value;
    value = new CladogeneticProbabilityMatrix( num_states );
    
    // get speciation rates and the clado events
    const std::vector<double>& sr = probabilities->getValue();
    const RbVector<RbVector<long> >& events = cladogenetic_events->getValue();
    
    if (sr.size() != events.size())
    {
        throw RbException("You must enter the same number of cladogenetic events and speciation rates.");
    }
    
    // get number of events
    size_t num_events = events.size();
    
    // get sum of unnormalized probabilities per ancestral state
    std::map<unsigned, double> prob_sum;
    
    // for each clado event type build a map in the structure:
    // pair< [ancestor_state, daughter_1_state, daughter_2_state], speciation_rate >
    
    for (size_t i = 0; i < num_events; i++)
    {
        if (events[i].size() != 3)
        {
            throw RbException("Invalid cladogenetic event type.");
        }
        
        std::vector<unsigned> idx(3);
        idx[0] = (unsigned)events[i][0];
        idx[1] = (unsigned)events[i][1];
        idx[2] = (unsigned)events[i][2];
        event_map[ idx ] = sr[i];
        
        if ( prob_sum.find( idx[0] ) == prob_sum.end() ) {
            prob_sum[ idx[0] ] = 0.0;
        }
        prob_sum[ idx[0] ] += sr[i];
    }
    
    
    // now, divide each probability by each anc. state's normalizing sum
    for (size_t i = 0; i < num_events; i++)
    {
        std::vector<unsigned> idx(3);
        idx[0] = (unsigned)events[i][0];
        idx[1] = (unsigned)events[i][1];
        idx[2] = (unsigned)events[i][2];
        event_map[ idx ] = event_map[ idx ] / prob_sum[ idx[0] ];
    }
    
    // assign normalized event map to function's value
    value->setEventMap(event_map);
}


void CladogeneticProbabilityMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == probabilities)
    {
        probabilities = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
