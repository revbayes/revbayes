//
//  ChromosomesCladogenicBirthDeathFunction.cpp
//
//  Created by Will Freyman on 6/22/16.
//


#include "ChromosomesCladogenicBirthDeathFunction.h"

#include <utility>

#include "CladogeneticSpeciationRateMatrix.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class BranchHistory; }
namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


//TypedFunction<MatrixReal>( new MatrixReal( mc + 1, (mc + 1) * (mc + 1), 0.0 ) ),
ChromosomesCladogenicBirthDeathFunction::ChromosomesCladogenicBirthDeathFunction(const TypedDagNode< RbVector<double> > *sr, unsigned mc ):
TypedFunction<CladogeneticSpeciationRateMatrix>( new CladogeneticSpeciationRateMatrix( mc + 1 ) ),
speciationRates( sr ),
maxChromo(mc),
numEventTypes( (unsigned)sr->getValue().size() ),
use_hidden_rate(false)
{
    addParameter( speciationRates );
    
    // since the transition rate matrix will be very large (maxChromo by maxChromo^2) but
    // only sparsely filled, the we use an event map instead
    // of the full matrix for efficiency
    buildEventMap();
    update();
}


ChromosomesCladogenicBirthDeathFunction::~ChromosomesCladogenicBirthDeathFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


void ChromosomesCladogenicBirthDeathFunction::buildEventMap( void ) {
    
    eventMapCounts.resize(maxChromo + 1, std::vector<unsigned>(numEventTypes, 0));
    
    // for each ancestor state build a map of the possible events in the structure:
    // pair< [ancestor_state, daughter_1_state, daughter_2_state], transition_probability >
    std::vector<unsigned> idx(3);
    for (unsigned i = 1; i <= maxChromo; i++)
    {
        // set ancestor state
        idx[0] = i;
        
        // loop through all possible events
        for (unsigned j = 0; j < numEventTypes; j++)
        {
            // if both daughters have same number chromsomes as ancestor
            if (j == NO_CHANGE)
            {
                idx[1] = i;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(NO_CHANGE));
                eventMapCounts[ i ][ NO_CHANGE ] += 1;
                eventMap[ idx ] = 0.0;
            }
            // if one of the daughters gains a chromosome
            else if (j == FISSION && i + 1 <= maxChromo)
            {
                idx[1] = i + 1;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(FISSION));
                eventMapCounts[ i ][ FISSION ] += 1;
                eventMap[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i + 1;
                eventMapTypes[ idx ].push_back(unsigned(FISSION));
                eventMapCounts[ i ][ FISSION ] += 1;
                eventMap[ idx ] = 0.0;
            }
            // if one of the daughters loses a chromosome
            else if ( j == FUSION && i > 1 )
            {
                idx[1] = i - 1;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(FUSION));
                eventMapCounts[ i ][ FUSION ] += 1;
                eventMap[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i - 1;
                eventMapTypes[ idx ].push_back(unsigned(FUSION));
                eventMapCounts[ i ][ FUSION ] += 1;
                eventMap[ idx ] = 0.0;
            }
            // if one of the daughters undergoes polyploidizations
            else if ( j == POLYPLOIDIZATION && i * 2 <= maxChromo )
            {
                idx[1] = i * 2;
                idx[2] = i;
                eventMapTypes[ idx ].push_back(unsigned(POLYPLOIDIZATION));
                eventMapCounts[ i ][ POLYPLOIDIZATION ] += 1;
                eventMap[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i * 2;
                eventMapTypes[ idx ].push_back(unsigned(POLYPLOIDIZATION));
                eventMapCounts[ i ][ POLYPLOIDIZATION ] += 1;
                eventMap[ idx ] = 0.0;
            }
            // if one of the daughters undergoes demipolyploidizations
            else if ( j == DEMIPOLYPLOIDIZATION && i * 1.5 <= maxChromo && i > 1 )
            {
                if (i % 2 == 0)
                {
                    idx[1] = i * 1.5;
                    idx[2] = i;
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMap[ idx ] = 0.0;
                    
                    idx[1] = i;
                    idx[2] = i * 1.5;
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMap[ idx ] = 0.0;
                }
                else
                {
                    // round down
                    idx[1] = (unsigned)( (double)i * 1.5 - 0.5 );
                    idx[2] = i;
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMap[ idx ] = 0.0;
                    
                    idx[1] = i;
                    idx[2] = (unsigned)( (double)i * 1.5 - 0.5 );
                    eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                    eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                    eventMap[ idx ] = 0.0;
                    
                    if ( i * 1.5 + 0.5 <= maxChromo )
                    {
                        // round up
                        idx[1] = i;
                        idx[2] = (unsigned)( (double)i * 1.5 + 0.5 );
                        eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                        eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                        eventMap[ idx ] = 0.0;
                        
                        idx[1] = (unsigned)( (double)i * 1.5 + 0.5 );
                        idx[2] = i;
                        eventMapTypes[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION));
                        eventMapCounts[ i ][ DEMIPOLYPLOIDIZATION ] += 1;
                        eventMap[ idx ] = 0.0;
                    }
                }
            }
        }
    }
}


ChromosomesCladogenicBirthDeathFunction* ChromosomesCladogenicBirthDeathFunction::clone( void ) const
{
    return new ChromosomesCladogenicBirthDeathFunction( *this );
}


double ChromosomesCladogenicBirthDeathFunction::computeDataAugmentedCladogeneticLnProbability(const std::vector<BranchHistory*>& histories,
                                                                                              size_t node_index,
                                                                                              size_t left_index,
                                                                                              size_t right_index ) const
{
    throw RbException("ChromosomesCladogenicBirthDeathFunction::computeDataAugmentedCladogeneticLnProbability is not currently implemented.");
    double lnP = 0.0;
    return lnP;
    
}



std::map< std::vector<unsigned>, double >  ChromosomesCladogenicBirthDeathFunction::getEventMap(double t)
{
    return eventMap;
}

const std::map< std::vector<unsigned>, double >&  ChromosomesCladogenicBirthDeathFunction::getEventMap(double t) const
{
    return eventMap;
}


void ChromosomesCladogenicBirthDeathFunction::update( void )
{
    // reset the transition matrix
    delete value;
    
    // create temp variables for exiting speciation rates and cladogenetic event probabilities
    std::vector<double> speciation_rate_sum_per_state;
    CladogeneticProbabilityMatrix cladogenetic_probability_matrix;

    // check for a hidden rate category
    if (use_hidden_rate) { 
        value = new CladogeneticSpeciationRateMatrix( (maxChromo + 1) * 2 );
        cladogenetic_probability_matrix = CladogeneticProbabilityMatrix((maxChromo + 1) * 2);
        speciation_rate_sum_per_state = std::vector<double>( (maxChromo + 1) * 2, 0.0 );
    } else { 
        value = new CladogeneticSpeciationRateMatrix( maxChromo + 1 );
        cladogenetic_probability_matrix = CladogeneticProbabilityMatrix((maxChromo + 1));
        speciation_rate_sum_per_state = std::vector<double>( (maxChromo + 1), 0.0 );
    }
    
    const std::vector<double>& sr = speciationRates->getValue();
    
    // assign the correct rate to each event
    for (unsigned i = 1; i <= maxChromo; i++)
    {
        std::map<std::vector<unsigned>, std::vector<unsigned> >::iterator it;
        for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            if (idx[0] == i)
            {
                // reset all rates to 0.0
                eventMap[ idx ] = 0.0;
                if (use_hidden_rate == true) 
                {
                    std::vector<unsigned> idx_hidden(3);
                    idx_hidden[0] = idx[0] + maxChromo + 1;
                    idx_hidden[1] = idx[1] + maxChromo + 1;
                    idx_hidden[2] = idx[2] + maxChromo + 1;
                    eventMap[ idx_hidden ] = 0.0;
                }
            }
        }

        for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            if (idx[0] == i)
            {
                std::vector<unsigned> event_types = it->second;
                for (size_t e = 0; e < event_types.size(); e++)
                {
                    double speciation_rate = 0.0;
                
                    // check for NaN values
                    if (sr[ event_types[e] ] == sr[ event_types[e] ])
                    {
                        speciation_rate = sr[ event_types[e] ];
                    }
                    
                    // normalize for all possible instances of this event type
                    double v = ( speciation_rate / eventMapCounts[ i ][ event_types[e] ] );
                    
                    // save the rate in the event map
                    eventMap[ idx ] += v;
                    speciation_rate_sum_per_state[ idx[0] ] += v;
                    if (use_hidden_rate == true) 
                    {
                        std::vector<unsigned> idx_hidden(3);
                        idx_hidden[0] = idx[0] + maxChromo + 1;
                        idx_hidden[1] = idx[1] + maxChromo + 1;
                        idx_hidden[2] = idx[2] + maxChromo + 1;
                        const std::vector<double>& rate_multipliers = hiddenRateMultipliers->getValue();
                        eventMap[ idx_hidden ] += (v * rate_multipliers[0]);
                        speciation_rate_sum_per_state[ idx_hidden[0] ] += (v * rate_multipliers[0]);
                    }
                }
            }
        }
    }

    // populate TensorPhylo rate/prob structures
    std::map<std::vector<unsigned>, double> clado_prob_event_map = cladogenetic_probability_matrix.getEventMap();

    speciation_rate_sum_per_state[0] = 0.0;
    std::vector<unsigned> idx_null(3);
    idx_null[0] = 0;
    idx_null[1] = 0;
    idx_null[2] = 0;
    clado_prob_event_map[ idx_null ] = 1.0;

    for (std::map<std::vector<unsigned>, double>::iterator jt = eventMap.begin(); jt != eventMap.end(); jt++) {
        const std::vector<unsigned>& idx = jt->first;
        clado_prob_event_map[ idx ] = eventMap[ idx ] / speciation_rate_sum_per_state[ idx[0] ];
    }
    cladogenetic_probability_matrix.setEventMap(clado_prob_event_map);

    // done!
    value->setEventMap(eventMap);
    value->setCladogeneticProbabilityMatrix( cladogenetic_probability_matrix );
    value->setSpeciationRateSumPerState( speciation_rate_sum_per_state );
}


void ChromosomesCladogenicBirthDeathFunction::setRateMultipliers(const TypedDagNode< RbVector< double > >* rm) 
{
    hiddenRateMultipliers = rm;
    addParameter( hiddenRateMultipliers );
    use_hidden_rate = true;
    
    buildEventMap();
    update();
}


void ChromosomesCladogenicBirthDeathFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == speciationRates)
    {
        speciationRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == hiddenRateMultipliers)
    {
        hiddenRateMultipliers = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
