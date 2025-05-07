//
//  CladogeneticSpeciationRateMatrixFunction.cpp
//
//  Created by Will Freyman on 8/1/17.
//


#include "CladogeneticSpeciationRateMatrixFunction.h"

#include <cstddef>

#include "CladogeneticSpeciationRateMatrix.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


//TypedFunction<MatrixReal>( new MatrixReal( mc + 1, (mc + 1) * (mc + 1), 0.0 ) ),
CladogeneticSpeciationRateMatrixFunction::CladogeneticSpeciationRateMatrixFunction(const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<long> > >* events, const TypedDagNode<RevBayesCore::RbVector<double> >* spec_rates, int n_states):
TypedFunction<CladogeneticSpeciationRateMatrix>( new CladogeneticSpeciationRateMatrix( n_states ) ),
cladogenetic_events( events ), 
num_states( n_states ),
speciation_rates( spec_rates )
{
    addParameter( speciation_rates );

    // since the transition rate matrix will be very large (num_states by num_states^2) but
    // only sparsely filled, the we use an event map instead
    // of the full matrix for efficiency
    update();
}


CladogeneticSpeciationRateMatrixFunction::~CladogeneticSpeciationRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


CladogeneticSpeciationRateMatrixFunction* CladogeneticSpeciationRateMatrixFunction::clone( void ) const
{
    return new CladogeneticSpeciationRateMatrixFunction( *this );
}


std::map< std::vector<unsigned>, double >  CladogeneticSpeciationRateMatrixFunction::getEventMap(double t)
{
    return event_map;
}

const std::map< std::vector<unsigned>, double >&  CladogeneticSpeciationRateMatrixFunction::getEventMap(double t) const
{
    return event_map;
}


void CladogeneticSpeciationRateMatrixFunction::update( void )
{
    // reset the transition matrix
    delete value;
    value = new CladogeneticSpeciationRateMatrix( num_states );
   
    // create temp variables for exiting speciation rates and cladogenetic event probabilities
    std::vector<double> speciation_rate_sum_per_state( num_states, 0.0 );
    CladogeneticProbabilityMatrix cladogenetic_probability_matrix( num_states );

    // get speciation rates and the clado events
    const std::vector<double>& sr = speciation_rates->getValue();
    const RbVector<RbVector<long> >& events = cladogenetic_events->getValue();

    if (sr.size() != events.size())
    {
        throw RbException("You must enter the same number of cladogenetic events and speciation rates.");
    }

    // for each clado event type build a map in the structure:
    // pair< [ancestor_state, daughter_1_state, daughter_2_state], speciation_rate >
    size_t num_events = events.size();
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
        speciation_rate_sum_per_state[ idx[0] ] += sr[i];
    }
    
    // populate TensorPhylo rate/prob structures
    std::map<std::vector<unsigned>, double> clado_prob_event_map = cladogenetic_probability_matrix.getEventMap();
    for (std::map<std::vector<unsigned>, double>::iterator jt = event_map.begin(); jt != event_map.end(); jt++) {
        const std::vector<unsigned>& idx = jt->first;
        clado_prob_event_map[ idx ] = event_map[ idx ] / speciation_rate_sum_per_state[ idx[0] ];
    }
    cladogenetic_probability_matrix.setEventMap(clado_prob_event_map);

    // done!
    value->setEventMap(event_map);
    value->setCladogeneticProbabilityMatrix( cladogenetic_probability_matrix );
    value->setSpeciationRateSumPerState( speciation_rate_sum_per_state );

}


void CladogeneticSpeciationRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == speciation_rates)
    {
        speciation_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
