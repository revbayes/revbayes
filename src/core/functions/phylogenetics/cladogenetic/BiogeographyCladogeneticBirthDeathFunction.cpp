//
//  BiogeographyCladogeneticBirthDeathFunction.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 12/15/18.
//  Copyright © 2018 Michael Landis. All rights reserved.
//

#include "BiogeographyCladogeneticBirthDeathFunction.h"

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
BiogeographyCladogeneticBirthDeathFunction::BiogeographyCladogeneticBirthDeathFunction( const TypedDagNode< RbVector< double > >* sr,
                                                                                        TypedDagNode< RbVector<double> >* wf,
                                                                                        TypedDagNode< RbVector< RbVector<double> > >* bf,
                                                                                        unsigned mrs,
                                                                                        unsigned msss,
                                                                                        std::string ct ):
TypedFunction<CladogeneticSpeciationRateMatrix>( new CladogeneticSpeciationRateMatrix(  pow(2,mrs)-1) ),
speciationRates( sr ),
withinRegionFeatures( wf ),
betweenRegionFeatures( bf ),
numCharacters( (unsigned)wf->getValue().size() ),
num_states( 2 ),
numIntStates( pow(2,wf->getValue().size())-1 ),
maxRangeSize(mrs),
maxSubrangeSplitSize(msss),
numEventTypes( (unsigned)sr->getValue().size() ),
use_hidden_rate(false),
use_cutset_mean(true),
connectivityType( ct )
{
    addParameter( speciationRates );
    addParameter( withinRegionFeatures );
    addParameter( betweenRegionFeatures );
    
    if (numCharacters > 8)
    {
        std::cout << "Warning: analyses may be prohibitively slow for >8 regions.\n";
//        throw RbException(">10 characters currently unsupported");
    }
    
    buildBits();
    buildRanges(ranges, betweenRegionFeatures, true);
    
    numRanges = (unsigned)ranges.size();
//    numRanges++; // add one for the null range
    
    buildEventMap();
    if (connectivityType == "none")
    {
        ; // do nothing
    }
    else if (connectivityType == "cutset")
    {
        buildCutsets();
//        buildCutsetFactors();
    }
    buildBuddingRegions();
    buildEventMapFactors();
//    updateEventMapWeights();
    
    update();
    
}


BiogeographyCladogeneticBirthDeathFunction::~BiogeographyCladogeneticBirthDeathFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/*
 * This function returns the value of mask, but complementing the 1-valued bits in base
 * e.g.
 *      mask=001110100
 *      base=001100000
 *      ret =000010100
 */

std::vector<unsigned> BiogeographyCladogeneticBirthDeathFunction::bitAllopatryComplement( const std::vector<unsigned>& mask, const std::vector<unsigned>& base )
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

void BiogeographyCladogeneticBirthDeathFunction::bitCombinations(std::vector<std::vector<unsigned> >& comb, std::vector<unsigned> array, int i, std::vector<unsigned> accum)
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

unsigned BiogeographyCladogeneticBirthDeathFunction::bitsToState( const std::vector<unsigned>& b )
{
    return bitsToStatesByNumOn[b];
}

/*
 *  This function converts a bit vector into a string (mostly for printing)
 */

std::string BiogeographyCladogeneticBirthDeathFunction::bitsToString( const std::vector<unsigned>& b )
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

void BiogeographyCladogeneticBirthDeathFunction::buildBits( void )
{
    
    eventTypes.push_back("s");
    eventTypes.push_back("a");
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
//    bitsByNumOn[0].push_back(bits[0]); // commented out to ignore null range
    for (size_t i = 0; i < numIntStates; i++)
    {
        size_t m = i+1; // offset by one (no null range)
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
 * This function builds the allopatric cutset, which is defined
 * as the set of edges removed in order to create the bipartition
 */
void BiogeographyCladogeneticBirthDeathFunction::buildBuddingRegions( void ) {
    
    std::map< std::vector<unsigned>, unsigned>::iterator it;
    
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event
        std::vector<unsigned> idx = it->first;
        unsigned event_type = it->second;
        
        // for allopatry events
        if (event_type == SYMPATRY)
        {
            
            // get right and left bitsets
            const std::set<unsigned>& s1 = statesToBitsetsByNumOn[ idx[1] ];
            const std::set<unsigned>& s2 = statesToBitsetsByNumOn[ idx[2] ];
            
            // get bud area where new species emerges
            unsigned bud_area = ( s1.size() > 1 ? *s2.begin() : *s1.begin() );
            
            eventMapBuddingRegions[ idx ] = bud_area;
        }
//        eventMapCutsets[ idx ] = cutset;
    }
    
    return;
}

/*
 * This function builds the allopatric cutset, which is defined
 * as the set of edges removed in order to create the bipartition
 */
void BiogeographyCladogeneticBirthDeathFunction::buildCutsets( void ) {
    
    std::map< std::vector<unsigned>, unsigned>::iterator it;
    
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event
        std::vector<unsigned> idx = it->first;
        unsigned event_type = it->second;
        
        // get right and left bitsets
        const std::set<unsigned>& s1 = statesToBitsetsByNumOn[ idx[1] ];
        const std::set<unsigned>& s2 = statesToBitsetsByNumOn[ idx[2] ];
     
        // fill vector with edges to cut
        std::vector< std::vector<unsigned> > cutset;

        // for allopatry events
        if (event_type == ALLOPATRY)
        {
            // find the edges between regions in daughter ranges
            std::set<unsigned>::iterator jt, kt;
            for (jt = s1.begin(); jt != s1.end(); jt++)
            {
                for (kt = s2.begin(); kt != s2.end(); kt++)
                {
                    if ( *jt != *kt )
                    {
                        std::vector<unsigned> edge;
                        edge.push_back( *jt );
                        edge.push_back( *kt );
                        cutset.push_back( edge );
                    }
                }
            }
        }
        eventMapCutsets[ idx ] = cutset;
    }
    
    return;
}


/*
 *  This function populates the eventMap, eventMapTypes, and eventMapCounts
 *  so it may be rapidly filled with values when update() is called.
 */

void BiogeographyCladogeneticBirthDeathFunction::buildEventMap( void ) {
    
    // clear events
    eventMapCounts.clear();
    
    // get L,R states per A state
    std::vector<unsigned> idx(3);

    // loop over possible ranges
    for (std::set<unsigned>::iterator its = ranges.begin(); its != ranges.end(); its++)
    {
        unsigned i = *its;
        idx[0] = i;
        eventMapCounts[i] = std::vector<unsigned>(NUM_CLADO_EVENT_TYPES, 0);
        
#ifdef DEBUG_DEC
        std::cout << "State " << i << "\n";
        std::cout << "Bits  " << bitsToString(statesToBitsByNumOn[i]) << "\n";
#endif
        
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
        
        std::vector<unsigned> bl(numCharacters, 0);
        std::vector<unsigned> br(numCharacters, 0);
        
        // narrow sympatry
        if (sumBits(ba) == 1)
        {
            idx[1] = i;
            idx[2] = i;
            if (ranges.find(i) == ranges.end())
            {
                continue;
            }
            
            eventMapTypes[ idx ] = SYMPATRY;
            eventMapCounts[ i ][  SYMPATRY ] += 1;
            eventMap[ idx ] = 0.0;
            
#ifdef DEBUG_DEC
            std::cout << "Narrow sympatry\n";
            std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
            std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
            std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";
#endif
            
        }
        
        // subset/widespread sympatry
        else if (sumBits(ba) > 1)
        {
            idx[1] = i;
            idx[2] = i;
            if (ranges.find(i) == ranges.end())
            {
                continue;
            }
            
            if (eventStringToStateMap.find("s") != eventStringToStateMap.end()) {
                
#ifdef DEBUG_DEC
                std::cout << "Subset sympatry (L-trunk, R-bud)\n";
#endif
                
                // get set of possible sympatric events for L-trunk, R-bud
                for (size_t j = 0; j < on.size(); j++)
                {
                    br = std::vector<unsigned>(numCharacters, 0);
                    br[ on[j] ] = 1;
                    //                unsigned sr = bitsToState(br);
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
                    
#ifdef DEBUG_DEC
                    std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                    std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                    std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n\n";
#endif
                    
                    br[ on[j] ] = 0;
                }
                
                
#ifdef DEBUG_DEC
                std::cout << "Subset sympatry (L-bud, R-trunk)\n";
#endif
                
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
                    
#ifdef DEBUG_DEC
                    std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " "<< bitsToString(statesToBitsByNumOn[i]) << "\n";
                    std::cout << "L " << bitsToState(bl) << " "<< bitsToString(bl) << "\n";
                    std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";
#endif
                    
                    bl[ on[j] ] = 0;
                }
            }
            
            
            if (eventStringToStateMap.find("a") != eventStringToStateMap.end()) {
                
                // get set of possible allopatry events
                bl = ba;
                std::vector<std::vector<unsigned> > bc;
                bitCombinations(bc, ba, 0, std::vector<unsigned>());
                
#ifdef DEBUG_DEC
                std::cout << "Allopatry combinations\n";
                std::cout << "A " << bitsToState(ba) << " " << bitsToString(ba) << "\n";
#endif
                
                for (size_t j = 0; j < bc.size(); j++)
                {
                    
                    bl = bc[j];
                    br = bitAllopatryComplement(ba, bl);
                    
                    // limit max allopatric split size
                    //if (sumBits(bl)==1 || sumBits(br)==1)
                    if ( sumBits(bl) <= maxSubrangeSplitSize || sumBits(br) <= maxSubrangeSplitSize )
                    {
                        
//                        unsigned sa = bitsToStatesByNumOn[ba];
                        unsigned sl = bitsToStatesByNumOn[bl];
                        unsigned sr = bitsToStatesByNumOn[br];
                        idx[1] = sl;
                        idx[2] = sr;
                        
#ifdef DEBUG_DEC
                        std::cout << "L " << bitsToState(bl) << " " << bitsToString(bl) << "\n";
                        std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n";
#endif
                        
                        
                        eventMapTypes[ idx ] = ALLOPATRY;
                        eventMapCounts[ i ][  ALLOPATRY ] += 1;
                        eventMap[ idx ] = 0.0;
                        
#ifdef DEBUG_DEC
                        std::cout << "\n";
#endif
                    }
                }
            }
        }
        
        // jump dispersal
        if (eventStringToStateMap.find("j") != eventStringToStateMap.end()) {
            
#ifdef DEBUG_DEC
            std::cout << "Jump dispersal (L-trunk, R-bud)\n";
#endif
            
            // get set of possible jump dispersal events for L-trunk, R-bud
            for (size_t j = 0; j < off.size(); j++)
            {
                br = std::vector<unsigned>(numCharacters, 0);
                br[ off[j] ] = 1;
                //                unsigned sr = bitsToState(br);
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
                
#ifdef DEBUG_DEC
                std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n\n";
#endif
                
                br[ off[j] ] = 0;
            }
            
            
#ifdef DEBUG_DEC
            std::cout << "Jump dispersal (L-bud, R-trunk)\n";
#endif
            
            // get set of possible jump dispersal events for R-trunk, L-bud
            for (size_t j = 0; j < off.size(); j++)
            {
                bl = std::vector<unsigned>(numCharacters, 0);
                
                bl[ off[j] ] = 1;
                //                unsigned sl = bitsToState(bl);
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
                
#ifdef DEBUG_DEC
                std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " "<< bitsToString(statesToBitsByNumOn[i]) << "\n";
                std::cout << "L " << bitsToState(bl) << " "<< bitsToString(bl) << "\n";
                std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";
#endif
                
                bl[ off[j] ] = 0;
            }
        }
#ifdef DEBUG_DEC
        std::cout << "\n\n";
#endif
    }
#ifdef DEBUG_DEC
    //    for (size_t i = 0; i < eventMapCounts.size(); i++) {
    //        std::cout << bitsToState(statesToBitsByNumOn[i]) << " " << eventMapCounts[ i ] << "\n";
    //    }
    //
    std::cout << "------\n";
#endif
    
}


/*
 * This function precomputes the base factors for the event map. Modified values
 * of the factors are then applied to the model rates in the update() function.
 */

void BiogeographyCladogeneticBirthDeathFunction::buildEventMapFactors(void)
{
    
    std::vector<double> max_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    std::vector<double> sum_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    std::vector<double> ln_sum_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    std::vector<double> prod_value( NUM_CLADO_EVENT_TYPES, 1.0 );
    std::vector<double> mean_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    std::vector<double> geomean_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    std::vector<unsigned> n_value( NUM_CLADO_EVENT_TYPES, 0 );
    
    // loop over all events and their types
    std::map< std::vector<unsigned>, unsigned >::iterator it;
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event info
        std::vector<unsigned> idx = it->first;
        unsigned event_type = it->second;
        
        // get event score
        double v = 1.0;
        
        if (connectivityType == "none") {
            ; // do nothing
        }
        else if (connectivityType == "cutset") {
            v = computeCutsetScore(idx, event_type);
        }

        eventMapFactors[ idx ] = v;

        // get event map factor statistics for renormalization
        if ( v > max_value[event_type] )
        {
            max_value[event_type] = v;
        }
        if (v > 0.0) {
            sum_value[event_type] += v;
            ln_sum_value[event_type] += std::log(v);
            prod_value[event_type] *= v;
            n_value[event_type] += 1;
        }

    }
    
    for (size_t i = 0; i < mean_value.size(); i++) {
        mean_value[i] = 0.0;
        geomean_value[i] = 0.0;
        if (n_value[i] > 0) {
            mean_value[i] = sum_value[i] / n_value[i];
//            std::cout << i << " " << prod_value[i] << " " << n_value
//            geomean_value[i] = std::pow( prod_value[i], (1.0/n_value[i])); //std::exp( (1.0/n_value[i]) * ln_sum_value[i] );
            geomean_value[i] = std::exp( (1.0/n_value[i]) * ln_sum_value[i] );
        }
    }
    
    // normalize event factors by max factor of event type
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event info
        std::vector<unsigned> idx = it->first;
        unsigned event_type = it->second;

        eventMapFactors[ idx ] = eventMapFactors[ idx ] / geomean_value[ event_type ];
        eventMapWeights[ idx ] = eventMapFactors[ idx ];
    }
        
    return;
}

/*
 * This function builds all defined ranges in the model
 */

void BiogeographyCladogeneticBirthDeathFunction::buildRanges(std::set<unsigned>& range_set, const TypedDagNode< RbVector<RbVector<double> > >* g, bool all)
{
    
    std::set<std::set<unsigned> > r;
    for (size_t i = 0; i < numCharacters; i++)
    {
        std::set<unsigned> s;
        s.insert((unsigned)i);
        r.insert(s);
        buildRangesRecursively(s, r, maxRangeSize, g, all);
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
    
#ifdef DEBUG_DEC
    for (std::set<std::set<unsigned> >::iterator it = r.begin(); it != r.end(); it++)
    {
        for (std::set<unsigned>::iterator jt = it->begin(); jt != it->end(); jt++)
        {
            std::cout << *jt << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
#endif
}


/* 
 * This recursive function accumulates areas to build a range
 */

void BiogeographyCladogeneticBirthDeathFunction::buildRangesRecursively(std::set<unsigned> s, std::set<std::set<unsigned> >& r, size_t k, const TypedDagNode< RbVector<RbVector<double> > >* g, bool all)
{
    
    // add candidate range to list of ranges
    if (s.size() <= k)
        r.insert(s);
    
    // stop recursing if range equals max size, k
    if (s.size() == k)
        return;
    
    
    // otherwise, recurse along
    for (std::set<unsigned>::iterator it = s.begin(); it != s.end(); it++)
    {
        for (size_t i = 0; i < numCharacters; i++)
        {
            if (g->getValue()[*it][i] > 0 || all) {
                std::set<unsigned> t = s;
                t.insert((unsigned)i);
                if (r.find(t) == r.end())
                {
                    buildRangesRecursively(t, r, k, g, all);
                }
            }
        }
    }
}



BiogeographyCladogeneticBirthDeathFunction* BiogeographyCladogeneticBirthDeathFunction::clone( void ) const
{
    return new BiogeographyCladogeneticBirthDeathFunction( *this );
}


double BiogeographyCladogeneticBirthDeathFunction::computeDataAugmentedCladogeneticLnProbability(const std::vector<BranchHistory*>& histories,
                                                                                              size_t node_index,
                                                                                              size_t left_index,
                                                                                              size_t right_index ) const
{
    throw RbException("BiogeographyCladogeneticBirthDeathFunction::computeDataAugmentedCladogeneticLnProbability is not currently implemented.");
    double lnP = 0.0;
    return lnP;
    
}

/*
 * This function computes the cutset score for a cladogenetic outcome (optionally, divided by number of cut edges)
 */

double BiogeographyCladogeneticBirthDeathFunction::computeCutsetScore( std::vector<unsigned> idx, unsigned event_type)
{
    double cost = 0.0;
    
    // get value for connectivity mtx
    const RbVector<double>& wf = withinRegionFeatures->getValue();
    const RbVector<RbVector<double> >& bf = betweenRegionFeatures->getValue();
    
    // compute modularity score depending on event type
    if (event_type == SYMPATRY)
    {
        // sympatry depends on matrix diagonal value
        unsigned i = eventMapBuddingRegions[idx];
        cost = wf[i];
    }
    else if (event_type == ALLOPATRY)
    {
        
        std::string s0 = bitsToString( statesToBitsByNumOn[ idx[0] ] );
        std::string s1 = bitsToString( statesToBitsByNumOn[ idx[1] ] );
        std::string s2 = bitsToString( statesToBitsByNumOn[ idx[2] ] );
        
        
//        std::cout << "COST " << s0 << " -> " << s1 << " | " << s2 << "\n";
        // allopatry depends on inverse sum of cutset cost of edge weights
        const std::vector<std::vector<unsigned> >& cutset = eventMapCutsets[idx];
        for (size_t i = 0; i < cutset.size(); i++) {
            size_t v1 = cutset[i][0];
            size_t v2 = cutset[i][1];
            cost += (1.0 / bf[v1][v2]);
//            std::cout << "\t" << v1 << " -- " << v2 << " : " << bf[v1][v2] << "\n";
        }
        
        // take the inverse sum of costs
        cost = 1.0 / cost;
//        std::cout << " = " << cost << "\n\n";
    }
    return cost;
}
//
///*
// * This function computes the modularity score for a cladogenetic outcome
// */
//
//double BiogeographyCladogeneticBirthDeathFunction::computeModularityScore(std::vector<unsigned> idx, unsigned event_type)
//{
//    // return value
//    double Q = 0.0;
//    
//    // get value for connectivity mtx
//    const RbVector<RbVector<double> >& mtx = betweenRegionFeatures->getValue();
//    size_t n = mtx.size();
//    
//    // get left/right area sets
//    std::vector< std::set<unsigned> > daughter_ranges;
//    daughter_ranges.push_back( statesToBitsetsByNumOn[ idx[0] ] );
//    daughter_ranges.push_back( statesToBitsetsByNumOn[ idx[1] ] );
//    
//    // compute modularity score depending on event type
//    if (event_type == SYMPATRY)
//    {
//        
//        if (daughter_ranges[0].size() == 1 && daughter_ranges[1].size() == 1)
//        {
//            return 0.0;
//        }
//        
//        std::vector<double> z(n, 0.0);
//        
//        // get trunk range (larger range)
//        const std::set<unsigned>& s = ( daughter_ranges[0].size() > daughter_ranges[1].size() ? daughter_ranges[0] : daughter_ranges[1] );
//        
//        // compute z, the sum of ranges for the trunk range
//        std::set<unsigned>::iterator it1, it2;
//        for (it1 = s.begin(); it1 != s.end(); it1++)
//        {
//            for (it2 = s.begin(); it2 != s.end(); it2++)
//            {
//                if ( (*it1) != (*it2) ) {
//                    z[ *it1 ] += mtx[ *it1 ][ *it2 ];
//                }
//            }
//        }
//        
//        double z_sum = 0.0;
//        for (size_t i = 0; i < z.size(); i++) {
//            z_sum += z[i];
//        }
//        
//        for (size_t i = 0; i < n; i++)
//        {
//            for (size_t j = 0; j < n; j++)
//            {
//                if (i != j) {
//                    Q += mtx[i][j] - (z[i] * z[j] / (2 * z_sum));
//                }
//            }
//        }
//        
//    }
//    else if (event_type == ALLOPATRY)
//    {
//        std::vector<double> z(n, 0.0);
//        
//        // compute z, the sum of edges across daughter ranges
//        for (size_t i = 0; i < daughter_ranges.size(); i++) {
//            
//            // get one daughter range
//            const std::set<unsigned>& s = daughter_ranges[i];
//            
//            // sum connectivity scores
//            std::set<unsigned>::iterator it1, it2;
//            for (it1 = s.begin(); it1 != s.end(); it1++)
//            {
//                for (it2 = s.begin(); it2 != s.end(); it2++)
//                {
//                    if ( (*it1) != (*it2) ) {
//                        z[ *it1 ] += mtx[ *it1 ][ *it2 ];
//                    }
//                }
//            }
//        }
//        
//        double z_sum = 0.0;
//        for (size_t i = 0; i < z.size(); i++) {
//            z_sum += z[i];
//        }
//        
//        for (size_t i = 0; i < n; i++)
//        {
//            for (size_t j = 0; j < n; j++)
//            {
//                if (i != j) {
//                    Q += mtx[i][j] - (z[i] * z[j] / (2 * z_sum));
//                }
//            }
//        }
//        
//    }
//    
//    return Q;
//}

/*
 * Returns the eventMap container
 */

std::map< std::vector<unsigned>, double >  BiogeographyCladogeneticBirthDeathFunction::getEventMap(double t)
{
    return eventMap;
}

/*
 * Returns the eventMap container (const)
 */


const std::map< std::vector<unsigned>, double >&  BiogeographyCladogeneticBirthDeathFunction::getEventMap(double t) const
{
    return eventMap;
}

/*
 * Prints the event map -- for debugging mostly
 */

void BiogeographyCladogeneticBirthDeathFunction::printEventMap(std::map< std::vector< unsigned >, double > x)
{
    std::map< std::vector< unsigned >, double >::iterator it;
    
    for (size_t i = 0; i < 2; i++) {
        
        std::string clado_str = "WITHIN_SPECIATION";
        if (i == ALLOPATRY) {
            clado_str = "BETWEEN_SPECIATION";
        }
            
        std::cout << "Event type : " << clado_str << "\n";
        for (it = x.begin(); it != x.end(); it++)
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
 * Sets the hidden rate multipliers
 */

void BiogeographyCladogeneticBirthDeathFunction::setRateMultipliers(const TypedDagNode< RbVector< double > >* rm)
{
    
    if (rm != NULL) {
        hiddenRateMultipliers = rm;
        addParameter( hiddenRateMultipliers );
        use_hidden_rate = true;
        
        buildEventMap();
        update();
    }
}


/*
 *  Computes the sum of bits (how many bits are set to 1)
 */

unsigned BiogeographyCladogeneticBirthDeathFunction::sumBits(const std::vector<unsigned>& b)
{
    unsigned n = 0;
    for (int i = 0; i < b.size(); i++)
        n += b[i];
    return n;
}


/*
 *  Standard swap parameters for moves and monitors
 */

void BiogeographyCladogeneticBirthDeathFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == speciationRates)
    {
        speciationRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == hiddenRateMultipliers)
    {
        hiddenRateMultipliers = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    if (oldP == withinRegionFeatures)
    {
        withinRegionFeatures = static_cast<const TypedDagNode<RbVector<double> >* >( newP );
    }
    if (oldP == betweenRegionFeatures)
    {
        betweenRegionFeatures = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
}


/*
 *  Update the rates in eventMap container
 */

void BiogeographyCladogeneticBirthDeathFunction::update( void )
{
    // reset the transition matrix
    delete value;
    
    // create temp variables for exiting speciation rates and cladogenetic event probabilities
    std::vector<double> speciation_rate_sum_per_state;
    CladogeneticProbabilityMatrix cladogenetic_probability_matrix;
    
    // check for a hidden rate category
    if (use_hidden_rate) {
        value = new CladogeneticSpeciationRateMatrix( numRanges * 2 );
        cladogenetic_probability_matrix = CladogeneticProbabilityMatrix(numRanges * 2);
        speciation_rate_sum_per_state = std::vector<double>( numRanges * 2, 0.0 );
    } else {
        value = new CladogeneticSpeciationRateMatrix( numRanges );
        cladogenetic_probability_matrix = CladogeneticProbabilityMatrix(numRanges);
        speciation_rate_sum_per_state = std::vector<double>( numRanges, 0.0 );
    }
    
    // update modularity score
    buildEventMapFactors();

    // get speciation rates across cladogenetic events
    const std::vector<double>& sr = speciationRates->getValue();
    
    // assign the correct rate to each event
    std::map<std::vector<unsigned>, unsigned>::iterator it;
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        const std::vector<unsigned>& idx = it->first;
        eventMap[ idx ] = 0.0;
        if (use_hidden_rate == true)
        {
            std::vector<unsigned> idx_hidden(3);
            idx_hidden[0] = idx[0] + numRanges + 1;
            idx_hidden[1] = idx[1] + numRanges + 1;
            idx_hidden[2] = idx[2] + numRanges + 1;
            eventMap[ idx_hidden ] = 0.0;
        }
    }
    
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        const std::vector<unsigned>& idx = it->first;
        unsigned event_type = it->second;
        double speciation_rate = 0.0;
        
        // check for NaN values
        if (sr[ event_type ] == sr[ event_type ])
        {
            speciation_rate = sr[ event_type ];
        }
        
        // divide by two if asymmetric event
        double f_asymm = ( idx[1] == idx[2] ? 1.0 : 0.5 );
        
        // rescale by connectivity weight
        double c_weight = eventMapWeights[ idx ];
        
        // compute the cladogenetic event rate
        double clado_rate = speciation_rate * f_asymm * c_weight;

        // save the rate in the event map
        eventMap[ idx ] += clado_rate;
        speciation_rate_sum_per_state[ idx[0] ] += eventMap[ idx ];
        if (use_hidden_rate == true)
        {
            std::vector<unsigned> idx_hidden(3);
            idx_hidden[0] = idx[0] + numRanges + 1;
            idx_hidden[1] = idx[1] + numRanges + 1;
            idx_hidden[2] = idx[2] + numRanges + 1;
            const std::vector<double>& rate_multipliers = hiddenRateMultipliers->getValue();
            eventMap[ idx_hidden ] += (clado_rate * rate_multipliers[0]);
            speciation_rate_sum_per_state[ idx_hidden[0] ] += eventMap[ idx_hidden ];
        }
    }
    
    // populate TensorPhylo rate/prob structures
    std::map<std::vector<unsigned>, double> clado_prob_event_map = cladogenetic_probability_matrix.getEventMap();
    for (std::map<std::vector<unsigned>, double>::iterator jt = eventMap.begin(); jt != eventMap.end(); jt++) {
        const std::vector<unsigned>& idx = jt->first;
        // if the speciation rate for the state is zero
        // and does not involve state change, set outcome prob
        // equal to 1.0 (no change only)
        if (speciation_rate_sum_per_state[ idx[0] ] == 0 && (idx[0] == idx[1]) && (idx[0]==idx[2])) {
            clado_prob_event_map[ idx ] = 1.0;
        }
        // otherwise if speciation rate is zero but
        // does involve state change, set outcome prob
        // equal to 0.0 (change impossible)
        else if (speciation_rate_sum_per_state[ idx[0] ] == 0) {
            clado_prob_event_map[ idx ] = 0.0;
        }
        // lastly, if speciation rate for state is non-zero
        // compute the probability as the rate of the particular
        // event divided by the sum of rates departing the state
        else {
            clado_prob_event_map[ idx ] = eventMap[ idx ] / speciation_rate_sum_per_state[ idx[0] ];
        }
    }
    cladogenetic_probability_matrix.setEventMap(clado_prob_event_map);
    
    // done!
    value->setEventMap(eventMap);
    value->setCladogeneticProbabilityMatrix( cladogenetic_probability_matrix );
    value->setSpeciationRateSumPerState( speciation_rate_sum_per_state );
}

