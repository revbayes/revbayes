//#define DEBUG_DEC

#include "DECCladogeneticStateFunction.h"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

#include "BiogeographicCladoEvent.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "BranchHistory.h"
#include "CharacterEventDiscrete.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Simplex.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class CharacterEvent; }
namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

DECCladogeneticStateFunction::DECCladogeneticStateFunction(const TypedDagNode< Simplex > *ep,
                                                           const TypedDagNode<RbVector<RbVector<double> > >* cg,
                                                           const TypedDagNode<RbVector<RbVector<double> > >* vg,
                                                           unsigned nc,
                                                           unsigned ns,
                                                           std::vector<std::string> et,
                                                           bool epawa,
                                                           bool wa,
                                                           bool uv,
                                                           unsigned mrs):
//    TypedFunction<MatrixReal>( new MatrixReal( pow(ns,nc), pow(ns,nc*2), 0.0) ),
    TypedFunction<CladogeneticProbabilityMatrix>( new CladogeneticProbabilityMatrix( computeNumStates(nc, mrs) )),
    eventProbs( ep ),
    connectivityGraph( cg ),
    vicarianceGraph( vg ),
    numCharacters(nc),
    num_states(2),
    numIntStates(pow(num_states,nc)),
    numEventTypes( BiogeographicCladoEvent::NUM_STATES ),
    maxRangeSize( (mrs < 1 || mrs > nc ? nc : mrs)),
    eventTypes( et ),
    eventProbsAsWeightedAverages(epawa),
    wideAllopatry(wa),
    useVicariance(uv)
{
    // add the lambda parameter as a parent
    addParameter( eventProbs );
    addParameter( connectivityGraph );
    addParameter( vicarianceGraph );
    
    if (numCharacters <= MAX_NUM_AREAS)
    {
        buildBits();
        buildRanges(beforeRanges, connectivityGraph, false);
        buildRanges(afterRanges, vicarianceGraph, false);
        
        numRanges = (unsigned)beforeRanges.size();
        numRanges++; // add one for the null range
        
        buildEventMap();
    } else {
        std::stringstream ss;
        ss << "DECCladogeneticStateFunction only supports up to " << MAX_NUM_AREAS << " areas\n";
        throw RbException(ss.str());
    }
    
    update();
   
    value->setEventTypes( eventTypes );
}

DECCladogeneticStateFunction::~DECCladogeneticStateFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}

std::vector<unsigned> DECCladogeneticStateFunction::bitAllopatryComplement( const std::vector<unsigned>& mask, const std::vector<unsigned>& base )
{
    std::vector<unsigned> ret = mask;
    for (size_t i = 0; i < base.size(); i++)
    {
        if (base[i] == 1)
            ret[i] = 0;
    }
    return ret;
}

void DECCladogeneticStateFunction::bitCombinations(std::vector<std::vector<unsigned> >& comb, std::vector<unsigned> array, int i, std::vector<unsigned> accum)
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

unsigned DECCladogeneticStateFunction::sumBits(const std::vector<unsigned>& b)
{
    unsigned n = 0;
    for (int i = 0; i < b.size(); i++)
        n += b[i];
    return n;
}

unsigned DECCladogeneticStateFunction::bitsToState( const std::vector<unsigned>& b )
{
    return bitsToStatesByNumOn[b];
}

std::string DECCladogeneticStateFunction::bitsToString( const std::vector<unsigned>& b )
{
    std::stringstream ss;
    for (size_t i = 0; i < b.size(); i++)
    {
        ss << b[i];
    }
    return ss.str();
}

void DECCladogeneticStateFunction::buildBits( void )
{
    
    for (size_t i = 0; i < eventTypes.size(); i++) {
        if (eventTypes[i]=="s")
            eventStringToStateMap[ eventTypes[i] ] = BiogeographicCladoEvent::SYMPATRY_SUBSET;
        else if (eventTypes[i]=="a")
            eventStringToStateMap[ eventTypes[i] ] = BiogeographicCladoEvent::ALLOPATRY;
        else if (eventTypes[i]=="j")
            eventStringToStateMap[ eventTypes[i] ] = BiogeographicCladoEvent::JUMP_DISPERSAL;
        else if (eventTypes[i]=="f")
            eventStringToStateMap[ eventTypes[i] ] = BiogeographicCladoEvent::SYMPATRY_WIDESPREAD;
    }
    
    bitsByNumOn.resize(numCharacters+1);
    statesToBitsByNumOn.resize(numIntStates);
    bits = std::vector<std::vector<unsigned> >(numIntStates, std::vector<unsigned>(numCharacters, 0));
    bitsByNumOn[0].push_back(bits[0]);
    for (size_t i = 1; i < numIntStates; i++)
    {
        size_t m = i;
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
            statesToBitsByNumOn[k++] = bitsByNumOn[i][j];
        }
    }
    
    for (size_t i = 0; i < statesToBitsByNumOn.size(); i++)
    {
        bitsToStatesByNumOn[ statesToBitsByNumOn[i] ] = (unsigned)i;
    }

}

void DECCladogeneticStateFunction::buildEventMap( void ) {
    
//    eventMapCounts.resize(numRanges, std::vector<unsigned>(numEventTypes, 0));
    eventMapCounts.clear();
    std::map<std::vector<unsigned>, double> eventMapProbs;
    
//    statesToBitsByNumOn = bits;
//    bitsToStatesByNumOn = inverseBits;
    
    // get L,R states per A state
    std::vector<unsigned> idx(3);
//    for (unsigned i = 0; i < numRanges; i++) {
    for (std::set<unsigned>::iterator its = beforeRanges.begin(); its != beforeRanges.end(); its++)
    {
        unsigned i = *its;
        idx[0] = i;
        eventMapCounts[i] = std::vector<unsigned>(BiogeographicCladoEvent::NUM_STATES, 0);
        
        // only include supported ranges
//        if (beforeRanges.find(i) == beforeRanges.end())
//            continue;
   
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
            if (beforeRanges.find(i) == beforeRanges.end())
            {
                continue;
            }
            
//            eventMapTypes[ idx ] = BiogeographicCladoEvent::SYMPATRY_NARROW;
//            eventMapCounts[ i ][  BiogeographicCladoEvent::SYMPATRY_NARROW ] += 1;
            eventMapTypes[ idx ] = BiogeographicCladoEvent::SYMPATRY_SUBSET;
            eventMapCounts[ i ][  BiogeographicCladoEvent::SYMPATRY_SUBSET ] += 1;
            eventMapProbs[ idx ] = 0.0;
 
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
            if (beforeRanges.find(i) == beforeRanges.end())
            {
                continue;
            }
            
            
            if (eventStringToStateMap.find("f") != eventStringToStateMap.end()) {
                eventMapTypes[ idx ] = BiogeographicCladoEvent::SYMPATRY_WIDESPREAD;
                eventMapCounts[ i ][  BiogeographicCladoEvent::SYMPATRY_WIDESPREAD ] += 1;
                eventMapProbs[ idx ] = 0.0;
            
#ifdef DEBUG_DEC
                std::cout << "Widespread sympatry\n";
                std::cout << "A " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                std::cout << "L " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n";
                std::cout << "R " << bitsToState(statesToBitsByNumOn[i]) << " " << bitsToString(statesToBitsByNumOn[i]) << "\n\n";
#endif
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
                    
                    if (beforeRanges.find(sr) == beforeRanges.end())
                    {
                        br[ on[j] ] = 0;
                        continue;
                    }
                    
                    
                    eventMapTypes[ idx ] = BiogeographicCladoEvent::SYMPATRY_SUBSET;
                    eventMapCounts[ i ][  BiogeographicCladoEvent::SYMPATRY_SUBSET ] += 1;
                    eventMapProbs[ idx ] = 0.0;
                
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
                    
                    if (beforeRanges.find(sl) == beforeRanges.end())
                    {
                        bl[ on[j] ] = 0;
                        continue;
                    }
                    
                    eventMapTypes[ idx ] =  BiogeographicCladoEvent::SYMPATRY_SUBSET;
                    eventMapCounts[ i ][  BiogeographicCladoEvent::SYMPATRY_SUBSET ] += 1;
                    eventMapProbs[ idx ] = 0.0;
             
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
                    
                    if (sumBits(bl)==1 || sumBits(br)==1 || wideAllopatry)
                    {
                        
    //                    unsigned sl = bitsToState(bl);
    //                    unsigned sr = bitsToState(br);
                        unsigned sa = bitsToStatesByNumOn[ba];
                        unsigned sl = bitsToStatesByNumOn[bl];
                        unsigned sr = bitsToStatesByNumOn[br];
                        idx[1] = sl;
                        idx[2] = sr;
                        
    //                    if (beforeRanges.find(sl) == beforeRanges.end() || beforeRanges.find(sr) == beforeRanges.end())
    //                    {
    //                        continue;
    //                    }
                        
//                        std::cout << sa << " -> " << sl << " | " << sr << "\n";
//                        std::cout << idx[0] << " -> " << idx[1] << " | " << idx[2] << "\n";
                        
                        // one daughter lineage must be contained within the vicariance range set
                        if (useVicariance)
                        {
//                            if (afterRanges.find(sl) != afterRanges.end() && afterRanges.find(sr) != afterRanges.end()) {
    //                            std::cout << sl << " " << sr << "\n";
    //                            for (std::set<unsigned>::iterator it_after = afterRanges.begin(); it_after != afterRanges.end(); it_after++)
    //                                std::cout << *it_after << "\n";
                            if (afterRanges.find(sa) != afterRanges.end())
//                                (afterRanges.find(sr) == afterRanges.end() ||
//                                 afterRanges.find(sl) == afterRanges.end()))
                            {
                                continue;
                            }
                        }

#ifdef DEBUG_DEC
                        std::cout << "L " << bitsToState(bl) << " " << bitsToString(bl) << "\n";
                        std::cout << "R " << bitsToState(br) << " " << bitsToString(br) << "\n";
#endif

                    
                        eventMapTypes[ idx ] = BiogeographicCladoEvent::ALLOPATRY;
                        eventMapCounts[ i ][  BiogeographicCladoEvent::ALLOPATRY ] += 1;
                        eventMapProbs[ idx ] = 0.0;

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
                
                if (beforeRanges.find(sr) == beforeRanges.end())
                {
                    br[ off[j] ] = 0;
                    continue;
                }
                
                
                eventMapTypes[ idx ] = BiogeographicCladoEvent::JUMP_DISPERSAL;
                eventMapCounts[ i ][  BiogeographicCladoEvent::JUMP_DISPERSAL ] += 1;
                eventMapProbs[ idx ] = 0.0;
                
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
                
                if (beforeRanges.find(sl) == beforeRanges.end())
                {
                    bl[ off[j] ] = 0;
                    continue;
                }
                
                eventMapTypes[ idx ] =  BiogeographicCladoEvent::JUMP_DISPERSAL;
                eventMapCounts[ i ][  BiogeographicCladoEvent::JUMP_DISPERSAL ] += 1;
                eventMapProbs[ idx ] = 0.0;
                
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

void DECCladogeneticStateFunction::buildRanges(std::set<unsigned>& range_set, const TypedDagNode< RbVector<RbVector<double> > >* g, bool all)
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

void DECCladogeneticStateFunction::buildRangesRecursively(std::set<unsigned> s, std::set<std::set<unsigned> >& r, size_t k, const TypedDagNode< RbVector<RbVector<double> > >* g, bool all)
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

DECCladogeneticStateFunction* DECCladogeneticStateFunction::clone( void ) const
{
    return new DECCladogeneticStateFunction( *this );
}

double DECCladogeneticStateFunction::computeDataAugmentedCladogeneticLnProbability( const std::vector<BranchHistory*>& histories,
                                                             size_t node_index,
                                                             size_t left_index,
                                                             size_t right_index ) const
{
    double lnP = 0.0;
    double p = 1.0;
    
    std::vector<CharacterEvent*>& nodeChildState = histories[ node_index ]->getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState = histories[ left_index ]->getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[ right_index ]->getParentCharacters();

//    std::cout << "nodeChildState   : ";
//    for (size_t i = 0; i < nodeChildState.size(); i++)
//    {
//        std::cout << static_cast<CharacterEventDiscrete*>(nodeChildState[i])->getState();
//    }
//    std::cout << "\n";
//    
//    std::cout << "leftParentState  : ";
//    for (size_t i = 0; i < leftParentState.size(); i++)
//    {
//        std::cout << static_cast<CharacterEventDiscrete*>(leftParentState[i])->getState();
//    }
//    std::cout << "\n";
//    
//    
//    std::cout << "rightParentState : ";
//    for (size_t i = 0; i < rightParentState.size(); i++)
//    {
//        std::cout << static_cast<CharacterEventDiscrete*>(rightParentState[i])->getState();
//    }
//    std::cout << "\n";
//    
//    histories[ node_index ]->print();
    
    // determine what type of cladogenetic event it is
    // determine the cladogenetic state based on un/shared areas
    size_t n_n_on = 0;
    size_t n_l_on = 0;
    size_t n_r_on = 0;
    size_t n_nlr_off = 0;
    size_t n_lr_mismatch = 0;
    size_t n_nlr_on = 0;
    size_t n_jump_mismatch = 0;
    for (size_t i = 0; i < nodeChildState.size(); i++)
    {
        size_t sn = static_cast<CharacterEventDiscrete*>(nodeChildState[i])->getState();
        size_t sl = static_cast<CharacterEventDiscrete*>(leftParentState[i])->getState();
        size_t sr = static_cast<CharacterEventDiscrete*>(rightParentState[i])->getState();
        
        n_n_on += sn;
        n_l_on += sl;
        n_r_on += sr;
        
        // copy sympatry
        if (sn == 1 && sl == 1 && sr == 1)
        {
            n_nlr_on += 1;
        }
        
        // jump dispersal
        if ( (sn==0&&sr==1&&sl==0) || (sn==0&&sr==0&&sl==1) )
        {
            n_jump_mismatch += 1;
        }
        
        // allopatry
        if ( (sn==1&&sl==1&&sr==0) || (sn==1&&sl==0&&sr==1) )
        {
            n_lr_mismatch += 1;
        }
        
        if (sn==0&&sr==0&&sl==0) {
            n_nlr_off += 1;
        }
        
    }
    
    
    if (n_nlr_on > 1) {
        
        throw RbException("unknown cladogenetic state");
    }
    
    
    std::string clado_type = "";
    if (n_n_on == 0) {
        clado_type = "n";
    }
    else if (n_nlr_on == 1 && n_nlr_off == (numCharacters-1))
    {
        clado_type = "s";
    }
    else if (n_n_on == n_lr_mismatch)
    {
        clado_type = "a";
    }
    else if (n_jump_mismatch == 1)
    {
        clado_type = "j";
    }
    else if (n_nlr_on == 1 && n_lr_mismatch==(n_n_on-n_nlr_on))
    {
        clado_type = "s";
    }
    else
    {
        
        throw RbException("Unknown cladogenetic event type!");
    }
    
    // get the information from the arguments for reading the file
    const std::vector<double>& ep = eventProbs->getValue();
    
    // get the probability for each clado event type
    std::map<std::string, double> probs;
    for (size_t i = 0; i < eventTypes.size(); i++)
    {
        probs[ eventTypes[i] ] = ep[ i ];
    }
    
    
    // the proposal prob
    if ( clado_type == "n" )
    {
        p = 1.0;
    }
//    else if ( clado_type == "s")
//    {
//        p = 1.0; // probs["s"];
//    }
    else if ( clado_type == "s" && n_n_on == 1)
    {
        p = 1.0;
    }
    else if ( clado_type == "s" && n_n_on > 1)
    {
        // Any single ancestral bit may be set across the two daughter ranges
        p = probs["s"] * (1.0 / (2 * n_n_on));
    }
    else if ( clado_type == "a" )
    {
        // Any combination of bits for one range (with the sister having its complement)
        // excluding the all-zero range and the all-one range (hence, -2)
        p = probs["a"] * (1.0 / (pow(2, n_n_on) - 2));
    }
    else if ( clado_type == "j" )
    {
        // Any single non-ancestral bit may be set across the two daughter ranges
        p = probs["j"] * (1.0 / (2 * (numCharacters - n_n_on)));
    }
    
    lnP = log(p);
    // compute the probability of that event type
    
    return lnP;
    
}


std::string DECCladogeneticStateFunction::simulateDataAugmentedCladogeneticState(std::vector<BranchHistory*>& histories,
                                                                          size_t node_index, size_t left_index, size_t right_index) const

{
    
    // what we will return
    std::string clado_state = "";
    
    // get containers for probabilities
    const std::vector<double>& ep = eventProbs->getValue();
    std::vector<double> probs(numEventTypes, 0.0);

    // sample cladogenetic state
    double u = GLOBAL_RNG->uniform01();
    std::string event_type = "";
    for (size_t i = 0; i < eventTypes.size(); i++)
    {
        event_type = eventTypes[i];
        std::map<std::string, unsigned>::const_iterator it = eventStringToStateMap.find( event_type );
        if (it == eventStringToStateMap.end()) {
            throw RbException("can't find clado event type");
        }
        unsigned clado_idx = it->second;
        probs[clado_idx] = ep[i];
        u -= probs[clado_idx];
        if (u < 0.0)
        {
            break;
        }
    }
    
    // get ranges
    std::vector<CharacterEvent*>& nodeChildState = histories[ node_index ]->getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState = histories[ left_index ]->getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[ right_index ]->getParentCharacters();
    
    // get indices that are on/off for nodeChildState
    std::vector<size_t> node_child_on;
    std::vector<size_t> node_child_off;
    
    // initialize daughter states to zero
    for (size_t site_index = 0; site_index < numCharacters; ++site_index)
    {
        static_cast<CharacterEventDiscrete*>( leftParentState[site_index] )->setState(0);
        static_cast<CharacterEventDiscrete*>( rightParentState[site_index] )->setState(0);
        size_t s = static_cast<CharacterEventDiscrete*>( nodeChildState[site_index] )->getState();
        if (s == 0) {
            node_child_off.push_back(site_index);
        }
        else {
            node_child_on.push_back(site_index);
        }
    }
    
    
    // sample cladogenetic state
    if (node_child_on.size() == 0) {
        clado_state = "n";

        ; // do nothing
    }
    else if ( node_child_on.size() == 1) {
        size_t s = node_child_on[0];
        static_cast<CharacterEventDiscrete*>( leftParentState[s] )->setState(1);
        static_cast<CharacterEventDiscrete*>( rightParentState[s] )->setState(1);
        clado_state = "s";
    }
    else
    {
        if (event_type == "s")
        {
            // choose the area inherited by the right node
            deprecated::random_shuffle( node_child_on.begin(), node_child_on.end());
            size_t area_idx = node_child_on[0];
            
            if (GLOBAL_RNG->uniform01() < 0.5)
            {
                
                // one area inherited by the right node
                static_cast<CharacterEventDiscrete*>( rightParentState[ area_idx ] )->setState(1);
                // all areas inherited by the left node
                for (std::vector<size_t>::iterator it = node_child_on.begin(); it != node_child_on.end(); it++)
                {
                    static_cast<CharacterEventDiscrete*>( leftParentState[ *it ] )->setState(1);
                }
            }
            else
            {
                // one area inherited by the left node
                static_cast<CharacterEventDiscrete*>( leftParentState[ area_idx ] )->setState(1);
                // all areas inherited by the left node
                for (std::vector<size_t>::iterator it = node_child_on.begin(); it != node_child_on.end(); it++)
                {
                    static_cast<CharacterEventDiscrete*>( rightParentState[ *it ] )->setState(1);
                }
            }
            
            clado_state = "s";
            
        }
        else if (event_type == "a")
        {
            // choose the area inherited by the right node
            deprecated::random_shuffle( node_child_on.begin(), node_child_on.end());
            size_t area_left_idx = node_child_on[0];
            size_t area_right_idx = node_child_on[1];
            
            static_cast<CharacterEventDiscrete*>( rightParentState[ area_left_idx ] )->setState(1);
            static_cast<CharacterEventDiscrete*>( leftParentState[ area_left_idx ] )->setState(0);
            static_cast<CharacterEventDiscrete*>( rightParentState[ area_right_idx ] )->setState(0);
            static_cast<CharacterEventDiscrete*>( leftParentState[ area_right_idx ] )->setState(1);
            
            for (size_t i = 2; i < node_child_on.size(); i++)
            {
                size_t j = node_child_on[i];
                if (GLOBAL_RNG->uniform01() < 0.5) {
                    static_cast<CharacterEventDiscrete*>( rightParentState[ j ] )->setState(0);
                    static_cast<CharacterEventDiscrete*>( leftParentState[ j ] )->setState(1);
                }
                else {
                    static_cast<CharacterEventDiscrete*>( leftParentState[ j ] )->setState(0);
                    static_cast<CharacterEventDiscrete*>( rightParentState[ j ] )->setState(1);
                }
            }
            clado_state = "a";
        }
        else if (event_type == "j")
        {
            // choose the area inherited by the right node
            deprecated::random_shuffle( node_child_off.begin(), node_child_off.end());
            size_t area_idx = node_child_off[0];
            
            if (GLOBAL_RNG->uniform01() < 0.5)
            {
                
                // one area inherited by the right node
                static_cast<CharacterEventDiscrete*>( rightParentState[ area_idx ] )->setState(1);
                // all areas inherited by the left node
                for (std::vector<size_t>::iterator it = node_child_on.begin(); it != node_child_on.end(); it++)
                {
                    static_cast<CharacterEventDiscrete*>( leftParentState[ *it ] )->setState(1);
                }
            }
            else
            {
                // one area inherited by the left node
                static_cast<CharacterEventDiscrete*>( leftParentState[ area_idx ] )->setState(1);
                // all areas inherited by the left node
                for (std::vector<size_t>::iterator it = node_child_on.begin(); it != node_child_on.end(); it++)
                {
                    static_cast<CharacterEventDiscrete*>( rightParentState[ *it ] )->setState(1);
                }
            }
            clado_state = "j";

        }
        else
        {
            
            throw RbException("Unknown cladogenetic event type sampled");
        }
    }
    
    return clado_state;
    
}


size_t DECCladogeneticStateFunction::computeNumStates(size_t numAreas, size_t maxRangeSize)
{
    size_t numStates = 1;
    for (size_t i = 1; i <= maxRangeSize; i++)
    {
        numStates += RbMath::choose(numAreas, i);
    }
    
    return numStates;
}

std::map< std::vector<unsigned>, double > DECCladogeneticStateFunction::getEventMap(double t)
{
    return this->getValue().getEventMap(0);
}

const std::map< std::vector<unsigned>, double >& DECCladogeneticStateFunction::getEventMap(double t) const
{
    return this->getValue().getEventMap(0);
}

const std::vector<std::string>& DECCladogeneticStateFunction::getEventTypes(void) const
{
    return eventTypes;
}


void DECCladogeneticStateFunction::update( void )
{
    
    // tmp
    std::map<std::vector<unsigned>, unsigned>::iterator it;
    std::map<std::vector<unsigned>, double> eventMapProbs;
    
    // get the information from the arguments for reading the file
    const std::vector<double>& ep = eventProbs->getValue();
    
    // get the probability for each clado event type
    std::vector<double> probs(numEventTypes, 0.0);
    for (size_t i = 0; i < eventTypes.size(); i++)
    {
        probs[ eventStringToStateMap[eventTypes[i]] ] = ep[i];
    }

    if (numCharacters > MAX_NUM_AREAS) return;
   
    if (eventProbsAsWeightedAverages)
    {
        // get sum of transition weights per starting state
        std::vector<double> z_pattern_probs( numRanges, 0.0 );
        for (size_t i = 0; i < eventTypes.size(); i++)
        {
            size_t k = eventStringToStateMap[ eventTypes[i] ];
            for (std::set<unsigned>::iterator its = beforeRanges.begin(); its != beforeRanges.end(); its++)
            {
                unsigned j = *its;
                z_pattern_probs[j] += probs[k] * eventMapCounts[j][k];
            }
        }
        for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            eventMapProbs[ idx ] = probs[ it->second ] / z_pattern_probs[ idx[0] ];
        }

    }
    else
    {
        // get transition probabilities per event class
        std::vector<double> z_class_probs( numRanges, 0 );
        for (size_t i = 0; i < eventTypes.size(); i++)
        {
            size_t k = eventStringToStateMap[ eventTypes[i] ];
            for (std::set<unsigned>::iterator its = beforeRanges.begin(); its != beforeRanges.end(); its++)
            {
                unsigned j = *its;
                z_class_probs[j] += (eventMapCounts[j][k] > 0 ? probs[k] : 0);
            }
        }
       
        for (it = eventMapTypes.begin(); it != eventMapTypes.end(); it++)
        {
            const std::vector<unsigned>& idx = it->first;
            eventMapProbs[ idx ] = probs[ it->second ] / eventMapCounts[ idx[0] ][ it->second ] / z_class_probs[ idx[0] ];
        }
    }
    
    value->setEventMap(eventMapProbs);
    
    return;
}



void DECCladogeneticStateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == eventProbs)
    {
        eventProbs = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == connectivityGraph)
    {
        connectivityGraph = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >(newP);
    }
    else if (oldP == vicarianceGraph)
    {
        vicarianceGraph = static_cast<const TypedDagNode<RbVector<RbVector<double> > >* >(newP);
    }
    
}


