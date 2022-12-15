//
//  GeoSSECladogeneticBirthDeathFunction.cpp
//  revbayes-proj
//
//  Created by Michael R May on 12/14/22.
//  Copyright © 2022 Michael R May. All rights reserved.
//

#include "GeoSSECladogeneticBirthDeathFunction.h"

#include <cmath>
#include <iostream>
#include <utility>

#include "CladogeneticSpeciationRateMatrix.h"
#include "RbException.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class BranchHistory; }
namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;


//TypedFunction<MatrixReal>( new MatrixReal( mc + 1, (mc + 1) * (mc + 1), 0.0 ) ),
GeoSSECladogeneticBirthDeathFunction::GeoSSECladogeneticBirthDeathFunction( const TypedDagNode< RbVector< double > >* sr,
                                                                            const TypedDagNode< RbVector<double> >* ar):
    TypedFunction<CladogeneticSpeciationRateMatrix>( new CladogeneticSpeciationRateMatrix(pow(2, (unsigned)sr->getValue().size()) - 1) ),
    sympatryRates( sr ),
    allopatryRates( ar ),
    hasJumps( false ),
    jumpRates( NULL ),
    numCharacters( (unsigned)sr->getValue().size() ),
    numStates( 2 ),
    numIntStates( pow(2,sr->getValue().size())-1 )
{
    
    // add the parameters
    addParameter( sympatryRates );
    addParameter( allopatryRates );
    
    // make the bits
    buildBits();

    // build the ranges
    buildRanges(ranges, true);
    numRanges = (unsigned)ranges.size();
    
    // build the event map
    buildEventMap();
    // printEventMap();

    // make sure the object is updated
    update();
    
}


GeoSSECladogeneticBirthDeathFunction::~GeoSSECladogeneticBirthDeathFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/*
 * This function returns the value of mask, but complementing the 1-valued bits in base
 * e.g.
 *      mask=001110100
 *      base=001100000
 *      ret =000010100
 */

void GeoSSECladogeneticBirthDeathFunction::setJumpRates(const TypedDagNode< RbVector< RbVector<double> > >* jr)
{
    if (jr != NULL) {
        
        // set the flag to true
        hasJumps = true;

        // set the parameter
        jumpRates = jr;
        addParameter( jumpRates );
        
        buildEventMap();
        update();

    }
}

std::vector<unsigned> GeoSSECladogeneticBirthDeathFunction::bitAllopatryComplement( const std::vector<unsigned>& mask, const std::vector<unsigned>& base )
{
    std::vector<unsigned> ret = mask;
    for (size_t i = 0; i < base.size(); i++)
    {
        if (base[i] == 1)
            ret[i] = 0;
    }
    return ret;
}


/*
 *  This recursive function builds all possible 0/1 bit combinations for array
 */

void GeoSSECladogeneticBirthDeathFunction::bitCombinations(std::vector<std::vector<unsigned> >& comb, std::vector<unsigned> array, int i, std::vector<unsigned> accum)
{
    if (i == array.size()) // end recursion
    {
        unsigned n = sumBits(accum);
        
        if ( n == 0 || n == sumBits(array) )
            ;  // ignore all-0, all-1 vectors
        else
            comb.push_back(accum);
    }
    else
    {
        unsigned b = array[i];
        std::vector<unsigned> tmp0(accum);
        std::vector<unsigned> tmp1(accum);
        tmp0.push_back(0);
        bitCombinations(comb,array,i+1,tmp0);
        if (b == 1)
        {
            tmp1.push_back(1);
            bitCombinations(comb,array,i+1,tmp1);
        }
    }
}

/*
 * This function returns the state associated with a bit vector
 */

unsigned GeoSSECladogeneticBirthDeathFunction::bitsToState( const std::vector<unsigned>& b )
{
    return bitsToStatesByNumOn[b];
}

/*
 *  This function converts a bit vector into a string (mostly for printing)
 */

std::string GeoSSECladogeneticBirthDeathFunction::bitsToString( const std::vector<unsigned>& b )
{
    std::stringstream ss;
    for (size_t i = 0; i < b.size(); i++)
    {
        ss << b[i];
    }
    return ss.str();
}


/*
 * This function generates the interchangeable state <-> bits <-> area-set
 * containers that define the state space.
 */

void GeoSSECladogeneticBirthDeathFunction::buildBits( void )
{
    
    eventTypes.push_back("s");
    eventTypes.push_back("a");
    eventTypes.push_back("j");
    for (size_t i = 0; i < eventTypes.size(); i++) {
        if (eventTypes[i]=="s")
            eventStringToStateMap[ eventTypes[i] ] = SYMPATRY;
        else if (eventTypes[i]=="a")
            eventStringToStateMap[ eventTypes[i] ] = ALLOPATRY;
        else if (eventTypes[i]=="j")
            eventStringToStateMap[ eventTypes[i] ] = JUMP_DISPERSAL;
    }
    
    bitsByNumOn.resize(numCharacters+1);
    statesToBitsByNumOn.resize(numIntStates);
    statesToBitsetsByNumOn.resize(numIntStates);
    bits = std::vector<std::vector<unsigned> >(numIntStates, std::vector<unsigned>(numCharacters, 0));
    for (size_t i = 0; i < numIntStates; i++)
    {
        size_t m = i + 1; // offset by one (no null range)
        for (size_t j = 0; j < numCharacters; j++)
        {
            bits[i][j] = m % 2;
            m /= 2;
            if (m == 0)
                break;
        }
        size_t j = sumBits(bits[i]);
        bitsByNumOn[j].push_back(bits[i]);

    }
    
    for (size_t i = 0; i < numIntStates; i++)
    {
        inverseBits[ bits[i] ] = (unsigned)i;
    }
    
    // assign state to each bit vector, sorted by numOn
    size_t k = 0;
    for (size_t i = 0; i < bitsByNumOn.size(); i++)
    {
        for (size_t j = 0; j < bitsByNumOn[i].size(); j++)
        {
            // assign to presence-absence vector
            statesToBitsByNumOn[k] = bitsByNumOn[i][j];
            
            // assign to set of present areas
            std::set<unsigned> s;
            for (size_t m = 0; m < statesToBitsByNumOn[k].size(); m++)
            {
                if (statesToBitsByNumOn[k][m] == 1) {
                    s.insert( (unsigned)m );
                }
            }
            statesToBitsetsByNumOn[k] = s;
            k++;
        }
    }
    
    for (size_t i = 0; i < statesToBitsByNumOn.size(); i++)
    {
        bitsToStatesByNumOn[ statesToBitsByNumOn[i] ] = (unsigned)i;
    }
    
}

/*
 *  This function populates the eventMap, eventMapTypes, and eventMapCounts
 *  so it may be rapidly filled with values when update() is called.
 */

void GeoSSECladogeneticBirthDeathFunction::buildEventMap( void ) {
    
    // clear events
    eventMapCounts.clear();
    
    // get L,R states per A state
    std::vector<unsigned> idx(3);

    // loop over possible ranges
    for (std::set<unsigned>::iterator its = ranges.begin(); its != ranges.end(); its++)
    {

        // get the state and index
        unsigned i = *its;
        idx[0] = i;

        // std::cout << "State " << i << "\n";
        // std::cout << "Bits  " << bitsToString(statesToBitsByNumOn[i]) << "\n";
        eventMapCounts[i] = std::vector<unsigned>(NUM_CLADO_EVENT_TYPES, 0);

        // get on bits for A
        const std::vector<unsigned>& ba = statesToBitsByNumOn[i];
        std::vector<unsigned> on;
        std::vector<unsigned> off;
        for (unsigned j = 0; j < ba.size(); j++)
        {
            if (ba[j] == 1)
                on.push_back(j);
            else
                off.push_back(j);
        }

        // bits for the left and right descendants
        std::vector<unsigned> bl(numCharacters, 0);
        std::vector<unsigned> br(numCharacters, 0);

        if (sumBits(ba) == 1)
        { // only one ancestral area
            
            // sympatry
            // both descendants get same area
            idx[1] = i;
            idx[2] = i;
            if (ranges.find(i) == ranges.end())
            {
                continue;
            }

            eventMapTypes[ idx ] = SYMPATRY;
            eventMapCounts[ i ][  SYMPATRY ] += 1;
            eventMap[ idx ] = 0.0;
            eventMapArea[ idx ] = i;

            // std::cout << "Narrow sympatry\n";
            // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
            // std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
            // std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";

        }
        else if (sumBits(ba) > 1)
        {            
            // two or more ancestral areas
            idx[1] = i;
            idx[2] = i;
            if (ranges.find(i) == ranges.end())
            {
                continue;
            }
            
            // sympatry
            if (eventStringToStateMap.find("s") != eventStringToStateMap.end()) {
                
                // std::cout << "Subset sympatry (L-trunk, R-bud)\n";
                // get set of possible sympatric events for L-trunk, R-bud
                for (size_t j = 0; j < on.size(); j++)
                {
                    br = std::vector<unsigned>(numCharacters, 0);
                    br[ on[j] ] = 1;
                    unsigned sr = bitsToStatesByNumOn[br];
                    idx[1] = i;
                    idx[2] = sr;
                    
                    if (ranges.find(sr) == ranges.end())
                    {
                        br[ on[j] ] = 0;
                        continue;
                    }
                    
                    eventMapTypes[ idx ] = SYMPATRY;
                    eventMapCounts[ i ][  SYMPATRY ] += 1;
                    eventMap[ idx ] = 0.0;
                    eventMapArea[ idx ] = sr;
                    
                    // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                    // std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                    // std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n\n";
                    
                    br[ on[j] ] = 0;

                }
                
                // std::cout << "Subset sympatry (L-bud, R-trunk)\n";
                // get set of possible sympatric events for R-trunk, L-bud
                for (size_t j = 0; j < on.size(); j++)
                {
                    bl = std::vector<unsigned>(numCharacters, 0);
                    
                    bl[ on[j] ] = 1;
                    //                unsigned sl = bitsToState(bl);
                    unsigned sl = bitsToStatesByNumOn[bl];
                    idx[1] = sl;
                    idx[2] = i;
                    
                    if (ranges.find(sl) == ranges.end())
                    {
                        bl[ on[j] ] = 0;
                        continue;
                    }
                    
                    eventMapTypes[ idx ] =  SYMPATRY;
                    eventMapCounts[ i ][  SYMPATRY ] += 1;
                    eventMap[ idx ] = 0.0;
                    eventMapArea[ idx ] = sl;
                    
                    // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " "<< bitsToString(statesToBitsByNumOn[i]) << "\n";
                    // std::cout << "L " << bitsToState(bl) << " "<< bitsToString(bl) << "\n";
                    // std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";
                    
                    bl[ on[j] ] = 0;

                }

            }
            
            // allopatry
            if (eventStringToStateMap.find("a") != eventStringToStateMap.end()) {

                // std::cout << "Allopatry (L-trunk, R-bud)\n";
                // get set of possible allopatric events for L-trunk, R-bud
                bl = ba;
                for (size_t j = 0; j < on.size(); j++)
                {
                    
                    br = std::vector<unsigned>(numCharacters, 0);
                    
                    bl[ on[j] ] = 0;
                    br[ on[j] ] = 1;

                    unsigned sl = bitsToStatesByNumOn[bl];
                    unsigned sr = bitsToStatesByNumOn[br];
                    idx[1] = sl;
                    idx[2] = sr;
                                        
                    eventMapTypes[ idx ] =  ALLOPATRY;
                    eventMapCounts[ i ][  ALLOPATRY ] += 1;
                    eventMap[ idx ] = 0.0;
                    eventMapArea[ idx ] = sr;
                    
                    // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " "<< bitsToString(statesToBitsByNumOn[i]) << "\n";
                    // std::cout << "L " << bitsToState(bl) << " "<< bitsToString(bl) << "\n";
                    // std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n\n";
                    
                    bl[ on[j] ] = 1;
                    br[ on[j] ] = 0;

                }


                // std::cout << "Allopatry (R-trunk, L-bud)\n";
                // get set of possible allopatric events for R-trunk, L-bud
                br = ba;
                for (size_t j = 0; j < on.size(); j++)
                {
                    
                    bl = std::vector<unsigned>(numCharacters, 0);
                    
                    bl[ on[j] ] = 1;
                    br[ on[j] ] = 0;

                    unsigned sl = bitsToStatesByNumOn[bl];
                    unsigned sr = bitsToStatesByNumOn[br];
                    idx[1] = sl;
                    idx[2] = sr;
                                        
                    eventMapTypes[ idx ] =  ALLOPATRY;
                    eventMapCounts[ i ][  ALLOPATRY ] += 1;
                    eventMap[ idx ] = 0.0;
                    eventMapArea[ idx ] = sl;
                    
                    // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " "<< bitsToString(statesToBitsByNumOn[i]) << "\n";
                    // std::cout << "L " << bitsToState(bl) << " "<< bitsToString(bl) << "\n";
                    // std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n\n";
                    
                    bl[ on[j] ] = 0;
                    br[ on[j] ] = 1;
                    
                }

            }

        }

        if (sumBits(ba) < numCharacters) {

            // jump dispersal
            idx[1] = i;
            idx[2] = i;
            if (ranges.find(i) == ranges.end())
            {
                continue;
            }

            // std::cout << "Jump dispersal (L-trunk, R-bud)\n";
            // get set of possible jump dispersal events for L-trunk, R-bud
            for (size_t j = 0; j < off.size(); j++)
            {
                
                br = std::vector<unsigned>(numCharacters, 0);
                br[ off[j] ] = 1;
                unsigned sr = bitsToStatesByNumOn[br];
                idx[1] = i;
                idx[2] = sr;
                
                if (ranges.find(sr) == ranges.end())
                {
                    br[ off[j] ] = 0;
                    continue;
                }
                
                
                eventMapTypes[ idx ] = JUMP_DISPERSAL;
                eventMapCounts[ i ][  JUMP_DISPERSAL ] += 1;
                eventMap[ idx ] = 0.0;
                eventMapArea[ idx ] = sr;
                
                // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                // std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                // std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n\n";
                
                br[ off[j] ] = 0;

            }
            
            // std::cout << "Jump dispersal (L-bud, R-trunk)\n";
            // get set of possible jump dispersal events for R-trunk, L-bud
            for (size_t j = 0; j < off.size(); j++)
            {
                bl = std::vector<unsigned>(numCharacters, 0);
                
                bl[ off[j] ] = 1;
                unsigned sl = bitsToStatesByNumOn[bl];
                idx[1] = sl;
                idx[2] = i;
                
                if (ranges.find(sl) == ranges.end())
                {
                    bl[ off[j] ] = 0;
                    continue;
                }
                
                eventMapTypes[ idx ] =  JUMP_DISPERSAL;
                eventMapCounts[ i ][  JUMP_DISPERSAL ] += 1;
                eventMap[ idx ] = 0.0;
                eventMapArea[ idx ] = sl;
                
                // std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " "<< bitsToString(statesToBitsByNumOn[i]) << "\n";
                // std::cout << "L " << bitsToState(bl) << " "<< bitsToString(bl) << "\n";
                // std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";
                
                bl[ off[j] ] = 0;

            }

        }

    }
    
}

/*
 * This function builds all defined ranges in the model
 */

void GeoSSECladogeneticBirthDeathFunction::buildRanges(std::set<unsigned>& range_set, bool all)
{
    
    std::set<std::set<unsigned> > r;
    for (size_t i = 0; i < numCharacters; i++)
    {
        std::set<unsigned> s;
        s.insert((unsigned)i);
        r.insert(s);
        buildRangesRecursively(s, r, all);
    }
    
    for (std::set<std::set<unsigned> >::iterator it = r.begin(); it != r.end(); it++)
    {
        std::vector<unsigned> v(numCharacters, 0);
        for (std::set<unsigned>::iterator jt = it->begin(); jt != it->end(); jt++)
        {
            v[*jt] = 1;
        }
        range_set.insert( bitsToStatesByNumOn[v] );
    }

}


/* 
 * This recursive function accumulates areas to build a range
 */

void GeoSSECladogeneticBirthDeathFunction::buildRangesRecursively(std::set<unsigned> s, std::set<std::set<unsigned> >& r, bool all)
{
    
    // add candidate range to list of ranges
    if (s.size() <= numCharacters)
        r.insert(s);
    
    // stop recursing if range equals max size
    if (s.size() == numCharacters)
        return;
    
    // otherwise, recurse along
    for (std::set<unsigned>::iterator it = s.begin(); it != s.end(); it++)
    {
        for (size_t i = 0; i < numCharacters; i++)
        {
            std::set<unsigned> t = s;
            t.insert((unsigned)i);
            if (r.find(t) == r.end())
            {
                buildRangesRecursively(t, r, all);
            }
        }
    }
}



GeoSSECladogeneticBirthDeathFunction* GeoSSECladogeneticBirthDeathFunction::clone( void ) const
{
    return new GeoSSECladogeneticBirthDeathFunction( *this );
}


/*
 * Returns the eventMap container
 */

std::map< std::vector<unsigned>, double >  GeoSSECladogeneticBirthDeathFunction::getEventMap(double t)
{
    return eventMap;
}

/*
 * Returns the eventMap container (const)
 */


const std::map< std::vector<unsigned>, double >&  GeoSSECladogeneticBirthDeathFunction::getEventMap(double t) const
{
    return eventMap;
}

/*
 * Prints the event map -- for debugging mostly
 */

void GeoSSECladogeneticBirthDeathFunction::printEventMap()
{
    std::map< std::vector< unsigned >, double >::iterator it;
    
    for (size_t i = 0; i < 3; i++) {
        
        std::string clado_str = "SYMPATRY";
        if (i == ALLOPATRY) {
            clado_str = "ALLOPATRY";
        } else if (i == JUMP_DISPERSAL) {
            clado_str = "JUMP_DISPERSAL";
        }
            
        std::cout << "Event type : " << clado_str << "\n";
        for (it = eventMap.begin(); it != eventMap.end(); it++)
        {
            std::vector<unsigned> idx = it->first;
            double rate = it->second;
            unsigned event_type = eventMapTypes[ idx ];
            
            if (i == event_type) {
            
                std::vector<unsigned> b0= statesToBitsByNumOn[ idx[0] ];
                std::vector<unsigned> b1= statesToBitsByNumOn[ idx[1] ];
                std::vector<unsigned> b2= statesToBitsByNumOn[ idx[2] ];
                
                std::string s0 = bitsToString( b0 );
                std::string s1 = bitsToString( b1 );
                std::string s2 = bitsToString( b2 );
                
                std::cout << s0 << " -> " << s1 << " | " << s2 << " = " << rate << "\n";
            }
        }
        
        std::cout << "\n";
    }
    
}

/*
 *  Computes the sum of bits (how many bits are set to 1)
 */

unsigned GeoSSECladogeneticBirthDeathFunction::sumBits(const std::vector<unsigned>& b)
{
    unsigned n = 0;
    for (int i = 0; i < b.size(); i++)
        n += b[i];
    return n;
}


/*
 *  Standard swap parameters for moves and monitors
 */

void GeoSSECladogeneticBirthDeathFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == sympatryRates)
    {
        sympatryRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == allopatryRates)
    {
        allopatryRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == jumpRates)
    {
        jumpRates = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
}


/*
 *  Update the rates in eventMap container
 */

void GeoSSECladogeneticBirthDeathFunction::update( void )
{
    // reset the transition matrix
    delete value;

    // create temp variables for exiting speciation rates and cladogenetic event probabilities
    std::vector<double> speciation_rate_sum_per_state;
    CladogeneticProbabilityMatrix cladogenetic_probability_matrix;
    
    value = new CladogeneticSpeciationRateMatrix( numRanges );
    cladogenetic_probability_matrix = CladogeneticProbabilityMatrix(numRanges);
    speciation_rate_sum_per_state = std::vector<double>( numRanges, 0.0 );
    
    // get parameters
    const std::vector<double>& sr = sympatryRates->getValue();
    const std::vector<double>& ar = allopatryRates->getValue();

    // assign the correct rate to each event
    std::map<std::vector<unsigned>, unsigned>::iterator it;
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        const std::vector<unsigned>& idx = it->first;
        eventMap[ idx ] = 0.0;
    }

    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        
        // get the event
        const std::vector<unsigned>& idx = it->first;
        unsigned event_type = it->second;

        // get the speciation rate
        double speciation_rate = 0.0;
        unsigned event_area = eventMapArea[ idx ];
        if (event_type == SYMPATRY)
        {
            speciation_rate = sr.at(event_area);
        }
        else if (event_type == ALLOPATRY)
        {
            speciation_rate = ar.at(event_area);
        }
        else if (hasJumps & (event_type == JUMP_DISPERSAL))
        {
            
            // get the dispersal-rate matrix
            const RbVector<RbVector<double> >& jr = jumpRates->getValue();

            // get the source region
            std::vector<unsigned> source_region = statesToBitsByNumOn[ idx[0] ];

            // accumulate rates per potential source area
            for(size_t i = 0; i < numCharacters; ++i)
            {
                if (source_region[i] == 1)
                {
                    speciation_rate += jr[i][event_area];
                }
            }

        }

        // divide by two if asymmetric event
        double f_asymm = ( idx[1] == idx[2] ? 1.0 : 0.5 );

        // compute the rate
        double clado_rate = speciation_rate * f_asymm;
        
        // save the rate in the event map
        eventMap[ idx ] += clado_rate;
        speciation_rate_sum_per_state[ idx[0] ] += eventMap[ idx ];

    }
    
    // populate TensorPhylo rate/prob structures
    std::map<std::vector<unsigned>, double> clado_prob_event_map = cladogenetic_probability_matrix.getEventMap();
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

