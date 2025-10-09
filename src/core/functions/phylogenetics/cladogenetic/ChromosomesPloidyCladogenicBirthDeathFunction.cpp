//
//  ChromosomesPloidyCladogenicBirthDeathFunction.cpp
//
//  Created by Will Freyman on 5/22/17.
//


#include "ChromosomesPloidyCladogenicBirthDeathFunction.h"

#include <cstddef>
#include <utility>

#include "CladogeneticSpeciationRateMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


ChromosomesPloidyCladogenicBirthDeathFunction::ChromosomesPloidyCladogenicBirthDeathFunction(const TypedDagNode< RbVector<double> > *sr, unsigned mc):
TypedFunction<CladogeneticSpeciationRateMatrix>( new CladogeneticSpeciationRateMatrix( (2 * (mc + 1)) ) ),
speciation_rates( sr ),
max_chromo(mc),
num_event_types( 10 )
{
    addParameter( speciation_rates );
    
    // since the transition rate matrix will be very large (2 * max_chromo by 2 * max_chromo^2) but
    // only sparsely filled, the we use an event map instead
    // of the full matrix for efficiency
    buildEventMap();
    update();
}


ChromosomesPloidyCladogenicBirthDeathFunction::~ChromosomesPloidyCladogenicBirthDeathFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


void ChromosomesPloidyCladogenicBirthDeathFunction::buildEventMap( void ) {
    
    event_map_counts.resize((2 * (max_chromo + 1)), std::vector<unsigned>(num_event_types, 0));
    
    // for each ancestor state build a map of the possible events in the structure:
    // pair< [ancestor_state, daughter_1_state, daughter_2_state], transition_probability >
    std::vector<unsigned> idx(3);

    // make events for diploid lineages
    for (unsigned i = 1; i <= max_chromo; i++)
    {
        // set ancestor state
        idx[0] = i;
        
        // loop through all possible events
        for (unsigned j = 0; j < num_event_types; j++)
        {
            // if both daughters have same number chromsomes as ancestor
            if (j == NO_CHANGE_D)
            {
                idx[1] = i;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(NO_CHANGE_D));
                event_map_counts[ i ][ NO_CHANGE_D ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters gains a chromosome
            else if (j == FISSION_D && i + 1 <= max_chromo)
            {
                idx[1] = i + 1;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(FISSION_D));
                event_map_counts[ i ][ FISSION_D ] += 1;
                event_map[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i + 1;
                event_map_types[ idx ].push_back(unsigned(FISSION_D));
                event_map_counts[ i ][ FISSION_D ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters loses a chromosome
            else if ( j == FUSION_D && i > 1 )
            {
                idx[1] = i - 1;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(FUSION_D));
                event_map_counts[ i ][ FUSION_D ] += 1;
                event_map[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i - 1;
                event_map_types[ idx ].push_back(unsigned(FUSION_D));
                event_map_counts[ i ][ FUSION_D ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters undergoes polyploidizations
            else if ( j == POLYPLOIDIZATION_D && i * 2 <= max_chromo )
            {
                idx[1] = (i * 2) + max_chromo + 1;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(POLYPLOIDIZATION_D));
                event_map_counts[ i ][ POLYPLOIDIZATION_D ] += 1;
                event_map[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = (i * 2) + max_chromo + 1;
                event_map_types[ idx ].push_back(unsigned(POLYPLOIDIZATION_D));
                event_map_counts[ i ][ POLYPLOIDIZATION_D ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters undergoes demipolyploidizations
            else if ( j == DEMIPOLYPLOIDIZATION_D && i * 1.5 <= max_chromo && i > 1 )
            {
                if (i % 2 == 0)
                {
                    idx[1] = (i * 1.5) + max_chromo + 1;
                    idx[2] = i;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_D));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_D ] += 1;
                    event_map[ idx ] = 0.0;
                    
                    idx[1] = i;
                    idx[2] = (i * 1.5) + max_chromo + 1;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_D));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_D ] += 1;
                    event_map[ idx ] = 0.0;
                }
                else
                {
                    // round down
                    idx[1] = (unsigned)( (double)i * 1.5 - 0.5 ) + max_chromo + 1;
                    idx[2] = i;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_D));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_D ] += 1;
                    event_map[ idx ] = 0.0;
                    
                    idx[1] = i;
                    idx[2] = (unsigned)( (double)i * 1.5 - 0.5 ) + max_chromo + 1;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_D));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_D ] += 1;
                    event_map[ idx ] = 0.0;
                    
                    if ( i * 1.5 + 0.5 <= max_chromo )
                    {
                        // round up
                        idx[1] = i;
                        idx[2] = (unsigned)( (double)i * 1.5 + 0.5 ) + max_chromo + 1;
                        event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_D));
                        event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_D ] += 1;
                        event_map[ idx ] = 0.0;
                        
                        idx[1] = (unsigned)( (double)i * 1.5 + 0.5 ) + max_chromo + 1;
                        idx[2] = i;
                        event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_D));
                        event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_D ] += 1;
                        event_map[ idx ] = 0.0;
                    }
                }
            }
        }
    }
    
    // make events for polyploid lineages
    for (unsigned i = (max_chromo + 2); i < (2 * (max_chromo + 1)); i++)
    {
        // set ancestor state
        idx[0] = i;
        
        // loop through all possible events
        for (unsigned j = 0; j < num_event_types; j++)
        {
            // if both daughters have same number chromsomes as ancestor
            if (j == NO_CHANGE_P)
            {
                idx[1] = i;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(NO_CHANGE_P));
                event_map_counts[ i ][ NO_CHANGE_P ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters gains a chromosome
            else if (j == FISSION_P && i + 1 < (2 * (max_chromo + 1)))
            {
                idx[1] = i + 1;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(FISSION_P));
                event_map_counts[ i ][ FISSION_P ] += 1;
                event_map[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i + 1;
                event_map_types[ idx ].push_back(unsigned(FISSION_P));
                event_map_counts[ i ][ FISSION_P ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters loses a chromosome
            else if ( j == FUSION_P && i > (max_chromo + 2) )
            {
                idx[1] = i - 1;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(FUSION_P));
                event_map_counts[ i ][ FUSION_P ] += 1;
                event_map[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = i - 1;
                event_map_types[ idx ].push_back(unsigned(FUSION_P));
                event_map_counts[ i ][ FUSION_P ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters undergoes polyploidizations
            else if ( j == POLYPLOIDIZATION_P && (i - (max_chromo + 1)) * 2 <= max_chromo )
            {
                idx[1] = ((i - (max_chromo + 1)) * 2) + max_chromo + 1;
                idx[2] = i;
                event_map_types[ idx ].push_back(unsigned(POLYPLOIDIZATION_P));
                event_map_counts[ i ][ POLYPLOIDIZATION_P ] += 1;
                event_map[ idx ] = 0.0;
                
                idx[1] = i;
                idx[2] = ((i - (max_chromo + 1)) * 2) + max_chromo + 1;
                event_map_types[ idx ].push_back(unsigned(POLYPLOIDIZATION_P));
                event_map_counts[ i ][ POLYPLOIDIZATION_P ] += 1;
                event_map[ idx ] = 0.0;
            }
            // if one of the daughters undergoes demipolyploidizations
            else if ( j == DEMIPOLYPLOIDIZATION_P && (i - (max_chromo + 1)) * 1.5 <= max_chromo && i > (max_chromo + 2) )
            {
                if ((i - (max_chromo + 1)) % 2 == 0)
                {
                    idx[1] = ((i - (max_chromo + 1)) * 1.5) + max_chromo + 1;
                    idx[2] = i;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_P));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_P ] += 1;
                    event_map[ idx ] = 0.0;
                    
                    idx[1] = i;
                    idx[2] = ((i - (max_chromo + 1)) * 1.5) + max_chromo + 1;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_P));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_P ] += 1;
                    event_map[ idx ] = 0.0;
                }
                else
                {
                    // round down
                    idx[1] = (unsigned)( (double)(i - (max_chromo + 1)) * 1.5 - 0.5 ) + max_chromo + 1;
                    idx[2] = i;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_P));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_P ] += 1;
                    event_map[ idx ] = 0.0;
                    
                    idx[1] = i;
                    idx[2] = (unsigned)( (double)(i - (max_chromo + 1)) * 1.5 - 0.5 ) + max_chromo + 1;
                    event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_P));
                    event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_P ] += 1;
                    event_map[ idx ] = 0.0;
                    
                    if ( (i - (max_chromo + 1)) * 1.5 + 0.5 <= max_chromo )
                    {
                        // round up
                        idx[1] = i;
                        idx[2] = (unsigned)( (double)(i - (max_chromo + 1)) * 1.5 + 0.5 ) + max_chromo + 1;
                        event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_P));
                        event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_P ] += 1;
                        event_map[ idx ] = 0.0;
                        
                        idx[1] = (unsigned)( (double)(i - (max_chromo + 1)) * 1.5 + 0.5 ) + max_chromo + 1;
                        idx[2] = i;
                        event_map_types[ idx ].push_back(unsigned(DEMIPOLYPLOIDIZATION_P));
                        event_map_counts[ i ][ DEMIPOLYPLOIDIZATION_P ] += 1;
                        event_map[ idx ] = 0.0;
                    }
                }
            }
        }
    }
}


ChromosomesPloidyCladogenicBirthDeathFunction* ChromosomesPloidyCladogenicBirthDeathFunction::clone( void ) const
{
    return new ChromosomesPloidyCladogenicBirthDeathFunction( *this );
}


std::map< std::vector<unsigned>, double >  ChromosomesPloidyCladogenicBirthDeathFunction::getEventMap(double t)
{
    return event_map;
}

const std::map< std::vector<unsigned>, double >&  ChromosomesPloidyCladogenicBirthDeathFunction::getEventMap(double t) const
{
    return event_map;
}


void ChromosomesPloidyCladogenicBirthDeathFunction::update( void )
{
    // reset the transition matrix
    delete value;
    value = new CladogeneticSpeciationRateMatrix( (2 * (max_chromo + 1)) );
    
    const std::vector<double>& sr = speciation_rates->getValue();
    
    // assign the correct rate to each event
    for (unsigned i = 1; i <= (2 * (max_chromo + 1)); i++)
    {
        std::map<std::vector<unsigned>, std::vector<unsigned> >::iterator it;
        for (it = event_map_types.begin(); it != event_map_types.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            if (idx[0] == i)
            {
                // reset all rates to 0.0
                event_map[ idx ] = 0.0;
//                (*value)[ idx[0] ][ (2 * (max_chromo + 1)) * idx[1] + idx[2] ] = 0.0;
            }
        }

        for (it = event_map_types.begin(); it != event_map_types.end(); it++)
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
                    double v = ( speciation_rate / event_map_counts[ i ][ event_types[e] ] );
                    
                    // save the rate in the rate matrix
//                    (*value)[ idx[0] ][ (2 * (max_chromo + 1)) * idx[1] + idx[2] ] += v;
                    
                    // save the rate in the event map
                    event_map[ idx ] += v;
                }
            }
        }
    }
    value->setEventMap(event_map);
}


void ChromosomesPloidyCladogenicBirthDeathFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == speciation_rates)
    {
        speciation_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}
