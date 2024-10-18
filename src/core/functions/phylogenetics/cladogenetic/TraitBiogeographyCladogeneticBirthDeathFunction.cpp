//
//  TraitBiogeographyCladogeneticBirthDeathFunction.cpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/27/24.
//

#define DEBUG_TRAITFIG 1

#include "TraitBiogeographyCladogeneticBirthDeathFunction.h"

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

//TraitBiogeographyCladogeneticBirthDeathFunction::TraitBiogeographyCladogeneticBirthDeathFunction(
//    const TypedDagNode< RbVector< double > >* sr,
//    TypedDagNode< RbVector<double> >* m_w,
//    TypedDagNode< RbVector< RbVector<double> > >* bf,
//    unsigned mrs,
//    unsigned msss,
//    bool nss,
//    std::string ct) :
TraitBiogeographyCladogeneticBirthDeathFunction::TraitBiogeographyCladogeneticBirthDeathFunction(
    const TypedDagNode<RbVector<RbVector<double> > >* rw,
    const TypedDagNode<RbVector<RbVector<double> > >* rb,
    const TypedDagNode<RbVector<RbVector<RbVector<double> > > >* mw,
    const TypedDagNode<RbVector<RbVector<RbVector<RbVector<double> > > > >* mb,
    unsigned mrs, unsigned msss, bool nss, std::string ct) :
TypedFunction<CladogeneticSpeciationRateMatrix>( new CladogeneticSpeciationRateMatrix(  pow(2,mrs)-1) ),
rho_w( rw ),
rho_b( rb ),
m_w( mw ),
m_b( mb ),
numTraits( (unsigned)m_w->getValue().size() ),
numRegions( (unsigned)m_w->getValue()[0].size() ),
numTraitValues( rho_w->getValue()[0].size() ),
numIntStates( 0 ),
maxRangeSize(mrs),
maxSubrangeSplitSize(msss),
numEventTypes( 2 ),     // 2 types: within-region, between-region
use_cutset_mean(true),
normalize_split_scores(nss),
connectivityType( ct )
{
    addParameter( rho_w );
    addParameter( rho_b );
    addParameter( m_w );
    addParameter( m_b );
    
//    if (numRegions > 8) {
//        std::cout << "Warning: analyses may be prohibitively slow for >8 regions.\n";
//    }
    
    // total possible number of integer-valued states
    numTraitValues = 2;
    numTraitSets = pow(numTraitValues,numTraits);
    numIntStates = (pow(2,numRegions)-1) * numTraitSets;
    
    if (numTraits + numRegions >= 15) {
        std::cout << "Warning: TraitBiogeographyCladogeneticBirthDeathFunction not designed for models with numTraits + numRegions >= 15\n";
    }
    
    buildStateSpace();
    
    numCompositeStates = compositeBitsToStatesByNumOn.size();
    numRanges = numCompositeStates / numTraitSets;
    
    buildEventMap();
    
    if (connectivityType == "none")
    {
        ; // do nothing
    }
    else if (connectivityType == "cutset")
    {
        buildCutsets();
    }
    buildBuddingRegions();
    buildEventMapFactors();
    
    update();
    
}


TraitBiogeographyCladogeneticBirthDeathFunction::~TraitBiogeographyCladogeneticBirthDeathFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

/*
 * This function returns the value of mask, but complementing the 1-valued rangeBits in base
 * e.g.
 *      mask=001110100
 *      base=001100000
 *      ret =000010100
 */

std::vector<unsigned> TraitBiogeographyCladogeneticBirthDeathFunction::rangeBitComplement( const RangeBits& mask, const RangeBits& base )
{
    RangeBits ret = mask;
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

void TraitBiogeographyCladogeneticBirthDeathFunction::rangeBitCombinations(std::vector<RangeBits>& comb, RangeBits array, int i, RangeBits accum)
{
    if (i == array.size()) // end recursion
    {
        unsigned n = sumRangeBits(accum);
        
        if ( n == 0 || n == sumRangeBits(array) )
            ;  // ignore all-0, all-1 vectors
        else
            comb.push_back(accum);
    }
    else
    {
        unsigned b = array[i];
        RangeBits tmp0(accum);
        RangeBits tmp1(accum);
        tmp0.push_back(0);
        rangeBitCombinations(comb,array,i+1,tmp0);
        if (b == 1)
        {
            tmp1.push_back(1);
            rangeBitCombinations(comb,array,i+1,tmp1);
        }
    }
}

/*
 * This function returns the state associated with a bit vector
 */

unsigned TraitBiogeographyCladogeneticBirthDeathFunction::compositeBitsToState( const CompositeBits& b )
{
    return compositeBitsToStatesByNumOn[b];
}

/*
 *  This function converts a bit vector into a string (mostly for printing)
 */

std::string TraitBiogeographyCladogeneticBirthDeathFunction::compositeBitsToString( const CompositeBits& b )
{
    std::stringstream ss;
    for (size_t i = 0; i < b.first.size(); i++) {
        ss << b.first[i];
    }
    ss << ":";
    for (size_t j = 0; j < b.second.size(); j++) {
        ss << b.second[j];
    }
    return ss.str();
}


/*
 * This function generates the interchangeable state <-> rangeBits <-> region-trait
 * containers that define the state space.
 *
 * Example with 2 regions (AB), 2 binary traits (XY)
 *
 *      integers    regions     trait       state-vector
 *      0           A                       10,00
 *      1           A           X           10,10
 *      2           A            Y          10,01
 *      3           A           XY          10,11
 *      4            B                      01,00
 *      5            B          X           01,10
 *      6            B           Y          01,01
 *      7            B          XY          01,11
 *      8           AB                      11,00
 *      9           AB          X           11,10
 *      10          AB           Y          11,01
 *      11          AB          XY          11,11
 */

void TraitBiogeographyCladogeneticBirthDeathFunction::buildStateSpace( void )
{
    
    // determine event types
    eventTypes.push_back("s");
    eventTypes.push_back("a");
    for (size_t i = 0; i < eventTypes.size(); i++) {
        if (eventTypes[i]=="s")
            eventStringToStateMap[ eventTypes[i] ] = WITHIN_SPECIATION;
        else if (eventTypes[i]=="a")
            eventStringToStateMap[ eventTypes[i] ] = BETWEEN_SPECIATION;
        else if (eventTypes[i]=="j")
            eventStringToStateMap[ eventTypes[i] ] = FOUNDER_SPECIATION;
    }
    
    // find all region-sets (ranges)
    // e.g. rangeBitsByNumOn[2] contains all ranges w/ 2 regions occupied
    rangeBitsByNumOn.resize(numRegions+1);
    
    // integer state as key
    statesToRangeBitsByNumOn.resize(numIntStates);
    statesToRangeBitsetsByNumOn.resize(numIntStates);
    statesToTraitBitsByNumOn.resize(numIntStates);
    statesToTraitBitsetsByNumOn.resize(numIntStates);
    statesToCompositeBitsByNumOn.resize(numIntStates);
    statesToCompositeBitsetsByNumOn.resize(numIntStates);
    
    // generate all trait-bit patterns
    size_t numAllTraits = pow(2, numTraits);
    traitBits = std::vector<TraitBits>(numAllTraits, TraitBits(numTraits, 0));
    for (size_t i = 0; i < numAllTraits; i++)
    {
        size_t m = i; // do not offset by one (absence of all traits)
        for (size_t j = 0; j < numTraits; j++)
        {
            traitBits[i][j] = m % 2;
            m /= 2;
            if (m == 0)
                break;
        }
    }
    
    // generate all range-bit patterns
    size_t numAllRanges = pow(2, numRegions) - 1;
    rangeBits = std::vector<RangeBits>(numAllRanges, RangeBits(numRegions, 0));

    for (size_t i = 0; i < numAllRanges; i++)
    {
        size_t m = i+1; // offset by one (no null range)
        for (size_t j = 0; j < numRegions; j++)
        {
            rangeBits[i][j] = m % 2;
            m /= 2;
            if (m == 0)
                break;
        }
        size_t j = sumRangeBits(rangeBits[i]);
        
        if (j <= maxRangeSize) {
            rangeBitsByNumOn[j].push_back(rangeBits[i]);
        }
    }

    
    // assign state to each bit vector, sorted by numOn
    size_t idx = 0;
    // each set of ranges with i regions occupied
    for (size_t i = 0; i < rangeBitsByNumOn.size(); i++)
    {
        // each range with i regions occupied
        for (size_t j = 0; j < rangeBitsByNumOn[i].size(); j++)
        {
            // each trait combination for range with i regions occupied
            for (size_t k = 0; k < traitBits.size(); k++) {
                // assign to presence-absence vector
                statesToRangeBitsByNumOn[idx] = rangeBitsByNumOn[i][j];
                statesToTraitBitsByNumOn[idx] = traitBits[k];
                
                // assign to set of present regions
                RangeBitset s_range;
                for (size_t m = 0; m < statesToRangeBitsByNumOn[idx].size(); m++)
                {
                    if (statesToRangeBitsByNumOn[idx][m] == 1) {
                        s_range.insert( (unsigned)m );
                    }
                }
                statesToRangeBitsetsByNumOn[idx] = s_range;
                
                // assign to set of present traits
                TraitBitset s_trait;
                for (size_t m = 0; m < statesToTraitBitsByNumOn[idx].size(); m++)
                {
                    if (statesToTraitBitsByNumOn[idx][m] == 1) {
                        s_trait.insert( (unsigned)m );
                    }
                }
                statesToTraitBitsetsByNumOn[idx] = s_trait;
                
                // increment counter for integer state index
                idx++;
            }
        }
    }
    
    // create inverse look-up table
    for (size_t i = 0; i < statesToRangeBitsByNumOn.size(); i++)
    {
        CompositeBits p1( statesToRangeBitsByNumOn[i], statesToTraitBitsByNumOn[i] );
        compositeBitsToStatesByNumOn[p1] = (unsigned)i;
        statesToCompositeBitsByNumOn[i] = p1;
//        CompositeBitset p2( statesToRangeBitsetByNumOn[i], statesToTraitBitsetsByNumOn[i] );
        
    }
    
    
    return;
    
}

/*
 * This function builds the allopatric cutset, which is defined
 * as the set of edges removed in order to create the bipartition
 */
void TraitBiogeographyCladogeneticBirthDeathFunction::buildBuddingRegions( void ) {
    
    std::map<StateTriplet, unsigned>::iterator it;
    
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event
        StateTriplet idx = it->first;
        unsigned event_type = it->second;
        
        // for allopatry events
        if (event_type == WITHIN_SPECIATION)
        {
            
            // get right and left rangeBitsets
            const RangeBitset& s1 = statesToRangeBitsetsByNumOn[ idx[1] ];
            const RangeBitset& s2 = statesToRangeBitsetsByNumOn[ idx[2] ];
            
            // get bud area where new species emerges
            unsigned bud_area = ( s1.size() > 1 ? *s2.begin() : *s1.begin() );
            
            eventMapBuddingRegions[ idx ] = bud_area;
        }
    }
    
    return;
}

/*
 * This function builds the allopatric cutset, which is defined
 * as the set of edges removed in order to create the bipartition
 */
void TraitBiogeographyCladogeneticBirthDeathFunction::buildCutsets( void ) {
    
    std::map<StateTriplet, unsigned>::iterator it;
    
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event
        StateTriplet idx = it->first;
        unsigned event_type = it->second;
        
        // get right and left rangeBitsets
        const RangeBitset& s1 = statesToRangeBitsetsByNumOn[ idx[1] ];
        const RangeBitset& s2 = statesToRangeBitsetsByNumOn[ idx[2] ];
     
        // fill vector with edges to cut
        std::vector< std::vector<unsigned> > cutset;

        // for allopatry events
        if (event_type == BETWEEN_SPECIATION)
        {
            // find the edges between regions in daughter ranges
            RangeBitset::iterator jt, kt;
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

void TraitBiogeographyCladogeneticBirthDeathFunction::buildEventMap( void ) {
   
    // clear events
    eventMapCounts.clear();
    
    // get L,R states per A state
    StateTriplet idx(3);

    // loop over possible comp
    for (unsigned i = 0; i < compositeBitsToStatesByNumOn.size(); i++)
    {
        idx[0] = i;
        eventMapCounts[i] = std::vector<unsigned>(NUM_CLADO_EVENT_TYPES, 0);
        
#ifdef DEBUG_TRAITFIG
        std::cout << "State " << i << "\n";
        std::cout << "Bits  " << compositeBitsToString(statesToCompositeBitsByNumOn[i]) << "\n";
#endif
        
        // get on rangeBits for A
        const RangeBits& ba = statesToCompositeBitsByNumOn[i].first;
        const TraitBits& ta = statesToCompositeBitsByNumOn[i].second;
        
        // find regions that are on and off
        std::vector<unsigned> on_region;
        std::vector<unsigned> off_region;
        for (unsigned j = 0; j < ba.size(); j++)
        {
            if (ba[j] == 1)
                on_region.push_back(j);
            else
                off_region.push_back(j);
        }
        
        RangeBits bl(numRegions, 0);
        RangeBits br(numRegions, 0);
        
        // narrow sympatry
        if (sumRangeBits(ba) == 1)
        {
            idx[1] = i;
            idx[2] = i;
            
            eventMapTypes[ idx ] = WITHIN_SPECIATION;
            eventMapCounts[ i ][  WITHIN_SPECIATION ] += 1;
            eventMap[ idx ] = 0.0;
//            eventMapByTraits[ ta ][ idx ] = 0.0;
            
//#ifdef DEBUG_TRAITFIG
//            std::cout << "Narrow sympatry\n";
//            std::cout << "A " << rangeBitsToState(statesToCompositeBitsByNumOn[i]) << " " << compositeBitsToString(statesToCompositeBitsByNumOn[i]) << "\n";
//            std::cout << "L " << rangeBitsToState(statesToCompositeBitsByNumOn[i]) << " " << compositeBitsToString(statesToCompositeBitsByNumOn[i]) << "\n";
//            std::cout << "R " << rangeBitsToState(statesToCompositeBitsByNumOn[i]) << " " << compositeBitsToString(statesToCompositeBitsByNumOn[i]) << "\n\n";
//#endif
            
        }
        
        
        // subset/widespread sympatry
        else if (sumRangeBits(ba) > 1)
        {
            idx[1] = i;
            idx[2] = i;
            //            if (ranges.find(i) == ranges.end())
            //            {
            //                continue;
            //            }
            
            
#ifdef DEBUG_TRAITFIG
            std::cout << "Subset sympatry (L-trunk, R-bud)\n";
#endif
            
            // get set of possible sympatric events for L-trunk, R-bud
            for (size_t j = 0; j < on_region.size(); j++)
            {
                br = RangeBits(numRegions, 0);
                br[ on_region[j] ] = 1;
                CompositeBits p(br, ta);
                unsigned sr = compositeBitsToStatesByNumOn[p];
                idx[1] = (unsigned)i;
                idx[2] = sr;
                
                //                if (ranges.find(sr) == ranges.end())
                //                {
                //                    br[ on_region[j] ] = 0;
                //                    continue;
                //                }
                
                eventMapTypes[ idx ] = WITHIN_SPECIATION;
                eventMapCounts[ i ][  WITHIN_SPECIATION ] += 1;
                eventMap[ idx ] = 0.0;
//                eventMapByTraits[ ta ][ idx ] = 0.0;
                
                //#ifdef DEBUG_TRAITFIG
                //                std::cout << "A " << rangeBitsToState(statesToRangeBitsByNumOn[i]) << " " << compositeBitsToString(statesToRangeBitsByNumOn[i]) << "\n";
                //                std::cout << "L " << rangeBitsToState(statesToRangeBitsByNumOn[i]) << " " << compositeBitsToString(statesToRangeBitsByNumOn[i]) << "\n";
                //                std::cout << "R " << rangeBitsToState(br) << " " << compositeBitsToString(br) << "\n\n";
                //#endif
                
                br[ on_region[j] ] = 0;
            }
            
            
#ifdef DEBUG_TRAITFIG
            std::cout << "Subset sympatry (L-bud, R-trunk)\n";
#endif
            
            // get set of possible sympatric events for R-trunk, L-bud
            for (size_t j = 0; j < on_region.size(); j++)
            {
                bl = RangeBits(numRegions, 0);
                bl[ on_region[j] ] = 1;
                CompositeBits p(bl, ta);
                unsigned sl = compositeBitsToStatesByNumOn[p];
                idx[1] = sl;
                idx[2] = i;
                
                //                if (ranges.find(sl) == ranges.end())
                //                {
                //                    bl[ on_region[j] ] = 0;
                //                    continue;
                //                }
                
                eventMapTypes[ idx ] =  WITHIN_SPECIATION;
                eventMapCounts[ i ][  WITHIN_SPECIATION ] += 1;
                eventMap[ idx ] = 0.0;
//                eventMapByTraits[ ta ][ idx ] = 0.0;
                
                //#ifdef DEBUG_TRAITFIG
                //                std::cout << "A " << rangeBitsToState(statesToRangeBitsByNumOn[i]) << " "<< compositeBitsToString(statesToRangeBitsByNumOn[i]) << "\n";
                //                std::cout << "L " << rangeBitsToState(bl) << " "<< compositeBitsToString(bl) << "\n";
                //                std::cout << "R " << rangeBitsToState(statesToRangeBitsByNumOn[i]) << " " << compositeBitsToString(statesToRangeBitsByNumOn[i]) << "\n\n";
                //#endif
                
                bl[ on_region[j] ] = 0;
            }
            
            
            // get set of possible allopatry events
            bl = ba;
            std::vector<RangeBits> bc;
            rangeBitCombinations(bc, ba, 0, RangeBits());
            
            //#ifdef DEBUG_TRAITFIG
            //            std::cout << "Allopatry combinations\n";
            //            std::cout << "A " << rangeBitsToState(ba) << " " << compositeBitsToString(ba) << "\n";
            //#endif
            
            for (size_t j = 0; j < bc.size(); j++)
            {
                
                bl = bc[j];
                br = rangeBitComplement(ba, bl);
                
                // limit max allopatric split size
                if ( sumRangeBits(bl) <= maxSubrangeSplitSize || sumRangeBits(br) <= maxSubrangeSplitSize )
                {
                    // unsigned sa = compositeBitsToStatesByNumOn[ba];
                    CompositeBits pl(bl, ta);
                    CompositeBits pr(br, ta);
                    unsigned sl = compositeBitsToStatesByNumOn[pl];
                    unsigned sr = compositeBitsToStatesByNumOn[pr];
                    idx[1] = sl;
                    idx[2] = sr;
                    
                    //#ifdef DEBUG_TRAITFIG
                    //                    std::cout << "L " << rangeBitsToState(bl) << " " << compositeBitsToString(bl) << "\n";
                    //                    std::cout << "R " << rangeBitsToState(br) << " " << compositeBitsToString(br) << "\n";
                    //#endif
                    
                    
                    eventMapTypes[ idx ] = BETWEEN_SPECIATION;
                    eventMapCounts[ i ][  BETWEEN_SPECIATION ] += 1;
                    eventMap[ idx ] = 0.0;
//                    eventMapByTraits[ ta ][ idx ] = 0.0;
                    
//#ifdef DEBUG_TRAITFIG
//                    std::cout << "\n";
//#endif
                }
            }
        }
        
//#ifdef DEBUG_TRAITFIG
//        std::cout << "\n\n";
//#endif
    }
//#ifdef DEBUG_TRAITFIG
//    //    for (size_t i = 0; i < eventMapCounts.size(); i++) {
//    //        std::cout << rangeBitsToState(statesToRangeBitsByNumOn[i]) << " " << eventMapCounts[ i ] << "\n";
//    //    }
//    //
//    std::cout << "------\n";
//#endif
    
    return;
    
}


/*
 * This function precomputes the base factors for the event map. Modified values
 * of the factors are then applied to the model rates in the update() function.
 */

void TraitBiogeographyCladogeneticBirthDeathFunction::buildEventMapFactors(void)
{
    
    // containers for normalizing factors, but assumes no traits
    //    std::vector<double> max_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    //    std::vector<double> sum_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    //    std::vector<double> prod_value( NUM_CLADO_EVENT_TYPES, 1.0 );
    //    std::vector<double> mean_value( NUM_CLADO_EVENT_TYPES, 0.0 );
//    std::vector<std::vector<std::vector<double> > > ln_sum_value( NUM_CLADO_EVENT_TYPES, 0.0 );
    //    std::vector<std::vector<std::vector<double> > > geomean_value; // ( NUM_CLADO_EVENT_TYPES, 0.0 );
    //    std::vector<std::vector<std::vector<unsigned> > > n_value; // ( NUM_CLADO_EVENT_TYPES, 0 );

    std::vector<std::vector<std::vector<double> > > ln_sum_value( numTraits, std::vector<std::vector<double> >(numTraitValues, std::vector<double>(NUM_CLADO_EVENT_TYPES, 0.0 )));
    std::vector<std::vector<std::vector<double> > > geomean_value( numTraits, std::vector<std::vector<double> >(numTraitValues, std::vector<double>(NUM_CLADO_EVENT_TYPES, 0.0 )));
    std::vector<std::vector<std::vector<unsigned> > > n_value( numTraits, std::vector<std::vector<unsigned> >(numTraitValues, std::vector<unsigned>(NUM_CLADO_EVENT_TYPES, 0 )));
    
    // initialize containers to be right sizes
    eventMapFactors = std::vector<std::vector<std::map<StateTriplet, double> > >(numTraits, std::vector<std::map<StateTriplet, double> >(numTraitValues));
    eventMapWeights = std::vector<std::vector<std::map<StateTriplet, double> > >(numTraits, std::vector<std::map<StateTriplet, double> >(numTraitValues));
    
    
    // normalization factors for each trait combination
    std::map<TraitBits, double> trait_factors;
    
    // loop over all events and their types
    std::map<StateTriplet, unsigned >::iterator it;
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event info
        StateTriplet idx = it->first;
        unsigned event_type = it->second;
        
        // get TraitBits for ancestral state
        const TraitBits& tb = statesToTraitBitsByNumOn[idx[0]];
        
        // get event score
        double v = 1.0;
        
        for (size_t trait_idx = 0; trait_idx < numTraits; trait_idx++) {
            for (size_t val_idx = 0; val_idx < numTraitValues; val_idx++) {
                
                if (connectivityType == "none") {
                    ; // do nothing
                }
                else if (connectivityType == "cutset") {
                    v = computeCutsetScore(idx, event_type, trait_idx, val_idx);
                }
                
                eventMapFactors[ trait_idx ][ val_idx ][ idx ] = v;
                
                if (v > 0.0) {
        //            sum_value[event_type] += v;
        //            prod_value[event_type] *= v;
                    ln_sum_value[trait_idx][val_idx][event_type] += std::log(v);
                    n_value[trait_idx][val_idx][event_type] += 1;
                }
            }
        }

    
        // get event map factor statistics for renormalization
//        if ( v > max_value[event_type] )
//        {
//            max_value[event_type] = v;
//        }
//        if (v > 0.0) {
////            sum_value[event_type] += v;
////            prod_value[event_type] *= v;
//            ln_sum_value[tb][event_type] += std::log(v);
//            n_value[tb][event_type] += 1;
//        }

    }
    
//    for (auto it = geomean_value.begin(); it != geomean_value.end(); it++) {
//        const TraitBits& tb = it->first;
    for (size_t trait_idx = 0; trait_idx < numTraits; trait_idx++) {
        for (size_t val_idx = 0; val_idx < numTraitValues; val_idx++) {
            for (size_t i = 0; i < geomean_value[trait_idx][val_idx].size(); i++) {
    //            geomean_value[*it][i] = 0.0;
                if (n_value[trait_idx][val_idx][i] > 0) {
                    geomean_value[trait_idx][val_idx][i] = std::exp( (1.0/n_value[trait_idx][val_idx][i]) * ln_sum_value[trait_idx][val_idx][i] );
                }
            }
        }
        
    }
    
    // normalize event factors by max factor of event type
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        // get event info
        std::vector<unsigned> idx = it->first;
        unsigned event_type = it->second;
        
        for (size_t trait_idx = 0; trait_idx < numTraits; trait_idx++) {
            for (size_t val_idx = 0; val_idx < numTraitValues; val_idx++) {
                if (normalize_split_scores) {
                    eventMapFactors[trait_idx][val_idx][idx] = eventMapFactors[trait_idx][val_idx][idx] / geomean_value[trait_idx][val_idx][ event_type ];
                }
                eventMapWeights[trait_idx][val_idx][idx] = eventMapFactors[trait_idx][val_idx][idx];
            }
        }
        // get trait vector
//        const TraitBits& tb = statesToTraitBitsByNumOn[idx[0]];
    }
        
    return;
}


TraitBiogeographyCladogeneticBirthDeathFunction* TraitBiogeographyCladogeneticBirthDeathFunction::clone( void ) const
{
    return new TraitBiogeographyCladogeneticBirthDeathFunction( *this );
}


double TraitBiogeographyCladogeneticBirthDeathFunction::computeDataAugmentedCladogeneticLnProbability(const std::vector<BranchHistory*>& histories,
                                                                                              size_t node_index,
                                                                                              size_t left_index,
                                                                                              size_t right_index ) const
{
    throw RbException("TraitBiogeographyCladogeneticBirthDeathFunction::computeDataAugmentedCladogeneticLnProbability is not currently implemented.");
    double lnP = 0.0;
    return lnP;
    
}

/*
 * This function computes the cutset score for a cladogenetic outcome (optionally, divided by number of cut edges)
 */

double TraitBiogeographyCladogeneticBirthDeathFunction::computeCutsetScore(StateTriplet idx, unsigned event_type, size_t trait_idx, size_t val_idx)
{
    double cost = 0.0;
    
    CompositeBits cb = statesToCompositeBitsByNumOn[ idx[0] ];
    RangeBits rb = cb.first;
    TraitBits tb = cb.second;
    
    // compute score depending on event type
    if (event_type == WITHIN_SPECIATION)
    {
        const RbVector<RbVector<RbVector<double> > >& m_w_values = m_w->getValue();
        unsigned region_idx = eventMapBuddingRegions[idx];
        cost = 1.0;
//        for (size_t trait_idx = 0; trait_idx < m_w_values.size(); trait_idx++) {
//            size_t trait_val = tb[trait_idx];
        cost *= m_w_values[trait_idx][val_idx][region_idx];
//        }
    }
    else if (event_type == BETWEEN_SPECIATION)
    {
        const RbVector<RbVector<RbVector<RbVector<double> > > >& m_b_values = m_b->getValue();
        
        // allopatry depends on inverse sum of cutset cost of edge weights
        const std::vector<std::vector<unsigned> >& cutset = eventMapCutsets[idx];
        for (size_t i = 0; i < cutset.size(); i++) {
            // we want the weight for the edge of trait_idx that connects region v1 to v2
            size_t v1 = cutset[i][0];
            size_t v2 = cutset[i][1];
            
            double mb12 = m_b_values[trait_idx][val_idx][v1][v2];
            double mb21 = m_b_values[trait_idx][val_idx][v2][v1];
            double mean_m_b_values = (mb12 + mb21) / 2.0;
            cost += (1.0 / mean_m_b_values);
            
//            std::cout << "\t" << v1 << " -- " << v2 << " : " << bf[v1][v2] << "\n";
        }
        
        // take the inverse sum of costs
        cost = 1.0 / cost;
//        std::cout << " = " << cost << "\n\n";
    }
    return cost;
}


/*
 * Returns the eventMap container
 */

std::map< std::vector<unsigned>, double >  TraitBiogeographyCladogeneticBirthDeathFunction::getEventMap(double t)
{
    return eventMap;
}

/*
 * Returns the eventMap container (const)
 */


const std::map< std::vector<unsigned>, double >&  TraitBiogeographyCladogeneticBirthDeathFunction::getEventMap(double t) const
{
    return eventMap;
}

/*
 * Prints the event map -- for debugging mostly
 */

void TraitBiogeographyCladogeneticBirthDeathFunction::printEventMap(std::map<StateTriplet, double> x)
{
    std::map< std::vector< unsigned >, double >::iterator it;
    
    for (size_t i = 0; i < 2; i++) {
        
        std::string clado_str = "WITHIN_SPECIATION";
        if (i == BETWEEN_SPECIATION) {
            clado_str = "BETWEEN_SPECIATION";
        }
            
        std::cout << "Event type : " << clado_str << "\n";
        for (it = x.begin(); it != x.end(); it++)
        {
            std::vector<unsigned> idx = it->first;
            double rate = it->second;
            unsigned event_type = eventMapTypes[ idx ];
            
            if (i == event_type) {
                CompositeBits b0 = statesToCompositeBitsByNumOn[ idx[0] ];
                CompositeBits b1 = statesToCompositeBitsByNumOn[ idx[1] ];
                CompositeBits b2 = statesToCompositeBitsByNumOn[ idx[2] ];

                std::string s0 = compositeBitsToString( b0 );
                std::string s1 = compositeBitsToString( b1 );
                std::string s2 = compositeBitsToString( b2 );
                
                std::cout << s0 << " -> " << s1 << " | " << s2 << " = " << rate << "\n";
            }
        }
        
        std::cout << "\n";
    }
    
}



/*
 *  Computes the sum of rangeBits (how many rangeBits are set to 1)
 */

unsigned TraitBiogeographyCladogeneticBirthDeathFunction::sumRangeBits(const std::vector<unsigned>& b)
{
    unsigned n = 0;
    for (int i = 0; i < b.size(); i++)
        n += b[i];
    return n;
}


/*
 *  Standard swap parameters for moves and monitors
 */

void TraitBiogeographyCladogeneticBirthDeathFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == rho_w)
    {
        rho_w = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    if (oldP == rho_b)
    {
        rho_b = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >( newP );
    }
    if (oldP == m_w)
    {
        m_w = static_cast<const TypedDagNode<RbVector<RbVector<RbVector<double> > > >* >( newP );
    }
    if (oldP == m_b)
    {
        m_b = static_cast<const TypedDagNode<RbVector<RbVector<RbVector<RbVector<double> > > > >* >( newP );
    }
}


/*
 *  Update the rates in eventMap container
 */

void TraitBiogeographyCladogeneticBirthDeathFunction::update( void )
{
    // reset the transition matrix
    delete value;
    
    // create temp variables for exiting speciation rates and cladogenetic event probabilities
    std::vector<double> speciation_rate_sum_per_state;
    CladogeneticProbabilityMatrix cladogenetic_probability_matrix;
    
    // make new speciation rate matrix container
    value = new CladogeneticSpeciationRateMatrix( numCompositeStates );
    cladogenetic_probability_matrix = CladogeneticProbabilityMatrix( numCompositeStates );
    speciation_rate_sum_per_state = std::vector<double>( numCompositeStates, 0.0 );

    // update modularity score
    buildEventMapFactors();

    // get speciation rates across cladogenetic events
    const RbVector<RbVector<double> >& rho_w_values = rho_w->getValue();
    const RbVector<RbVector<double> >& rho_b_values = rho_b->getValue();
    
    // assign the correct rate to each event
    std::map<std::vector<unsigned>, unsigned>::iterator it;
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        const std::vector<unsigned>& idx = it->first;
        eventMap[ idx ] = 0.0;
    }
    
    for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
    {
        const StateTriplet& idx = it->first;
        const RangeBits& rb = statesToRangeBitsByNumOn[ idx[0] ];
        const TraitBits& tb = statesToTraitBitsByNumOn[ idx[0] ];
        unsigned event_type = it->second;
        double speciation_rate = 0.0;
        
        // divide by two if asymmetric event
        double f_asymm = ( idx[1] == idx[2] ? 1.0 : 0.5 );
        
        // add the absolute rate contributed by each trait-value pair
        for (size_t trait_idx = 0; trait_idx < tb.size(); trait_idx++) {
            size_t val_idx = tb[trait_idx];
            double base_rate = 0.0;
            
            // rescale by relative rate factor for (anc -> left, right)
            double rel_rate = eventMapWeights[trait_idx][val_idx][idx];
            
            if (event_type == WITHIN_SPECIATION) {
                base_rate = rho_w_values[trait_idx][val_idx];
            } else if (event_type == BETWEEN_SPECIATION) {
                base_rate = rho_b_values[trait_idx][val_idx];
            }
            
            // compute the cladogenetic event rate
            double clado_rate = base_rate * rel_rate * f_asymm;
            speciation_rate += clado_rate;
//            speciation_rate += speciation_rate_part * f_asymm * c_weight;
        }
        
        // save the rate in the event map
        eventMap[ idx ] += speciation_rate;
        speciation_rate_sum_per_state[ idx[0] ] += eventMap[ idx ];
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
    
    // print
    printEventMap( this->eventMap );
}

