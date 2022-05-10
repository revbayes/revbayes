#include "PoMoState.h"

#include <boost/lexical_cast.hpp>
#include <math.h>
#include <iostream>
#include <cstddef>
#include <string>

#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "StringUtilities.h"
#include "Cloneable.h"

using namespace RevBayesCore;



/** Constructor that sets the observation and the other fields */
PoMoState::PoMoState(size_t n, size_t vps, const std::string &s, const std::string &chr, size_t pos, const std::vector<double> &w) : DiscreteCharacterState( n + size_t(RbMath::kchoose2(n))*(vps-1) ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    virtual_population_size( vps ),
    num_raw_states( n ),
    num_pomo_states( n + size_t(RbMath::kchoose2(n))*(virtual_population_size-1) ),
    state( num_pomo_states ),
    chromosome( chr ),
    position( pos ),
    string_value(s)
{
    weights = w;
    
    setWeighted( weights.size() > 0 );
    
    if ( s != "" )
    {
        setState(s);
    }
    
}


void PoMoState::addState(const std::string &symbol)
{
    ++num_observed_states;
    
    std::string labels = getStateLabels();
    size_t pos = labels.find(symbol);
    
    state.set( pos );
    index_single_state = pos;
}



PoMoState* PoMoState::clone( void ) const
{
    
    return new PoMoState( *this );
}


const std::string& PoMoState::getChromosome( void ) const
{
    return chromosome;
}


std::string PoMoState::getDataType( void ) const
{
    
    return "PoMo";
}


RbBitSet PoMoState::getState(void) const
{
    
    return state;
}


std::string PoMoState::getStateLabels( void ) const
{
    
    static std::string labels = "A C G T ";
    std::string acgt( "ACGT" );
    std::vector< size_t > frequencies;
    int stepSize = 100 / virtual_population_size;
    for (size_t i = 1; i < virtual_population_size; ++i)
    {
        frequencies.push_back(i*stepSize);
    }
    for ( size_t k = 0; k < acgt.size(); ++k )
    {
        char ch = acgt[k];
        
        for ( size_t j = k+1; j < acgt.size(); ++j )
        {
            char ch2 = acgt[j];
            for (size_t i = 0; i < virtual_population_size-1; ++i)
            {
                labels += ch + boost::lexical_cast<std::string>(frequencies[i]) + ch2 + boost::lexical_cast<std::string>(frequencies[virtual_population_size - 2 - i]) + " ";
            }
        }
    }
    
    return labels;
}


size_t PoMoState::getPosition( void ) const
{
    return position;
}


std::string PoMoState::getStringValue(void) const
{
    
    if ( isMissingState() )
    {
        return "?";
    }
    
    if ( isGapState() )
    {
        return "-";
    }
    
    return string_value;
}


bool PoMoState::isGapState( void ) const
{
    return is_gap;
}


bool PoMoState::isMissingState( void ) const
{
    return is_missing;
}


const std::string& PoMoState::nexusSeparator(void) const
{
    static std::string sep = " ";
    
    return sep;
}


void PoMoState::populateWeightedStatesForMonoallelicState(size_t ind_first, int sum)
{
    // We do PoMo state averaging.
    double n = (double)sum;
    double p = 1.0;
    
    weights[ind_first] = 1.0;
    state.set( ind_first );
    int vps_minus_1 = virtual_population_size - 1;
    
    
    std::vector<double> nd (vps_minus_1, 0.0);
    std::vector<double> id (vps_minus_1, 0.0);
    for (size_t i = 0; i < vps_minus_1; ++i)
    {
        nd[i] = (double)(i + 1) / (double)virtual_population_size;
        id[i] = (double)(vps_minus_1 - i) / (double)virtual_population_size;
    }
    
    if (ind_first < num_raw_states )
    {
        
        for (size_t ind_second=0; ind_second<num_raw_states; ++ind_second)
        {
            if ( ind_first != ind_second )
            {
                size_t first  = (ind_first < ind_second ? ind_first : ind_second);
                size_t second = (ind_first > ind_second ? ind_first : ind_second);
                
                size_t index = num_raw_states*first;
                for (size_t i=0; i<first; ++i)
                {
                    ind_first -= i;
                }
                index *= vps_minus_1;
                index += num_raw_states;
                
                for (size_t offset = 0; offset< vps_minus_1; ++offset)
                {
                    weights[ index + offset ] = pow(nd[offset], (double)sum);//RbStatistics::Binomial::pdf(n, p, (double)(numid1));
                    state.set( index + offset );
                }
            }
        }
        
    }
    else
    {
        throw RbException( "PoMo string state not correct. We found "+ StringUtilities::to_string(ind_first)  );
    }
    
    for (size_t i =0; i < weights.size(); ++i)
    {
        if (weights[i] < 1e-8)
        {
            weights[i] = 1e-8;
        }
    }
    
    
    return;
}


void PoMoState::setChromosome(std::string chr)
{
    chromosome = chr;
}


void PoMoState::setGapState( bool tf )
{
    is_gap = tf;
}


void PoMoState::setMissingState( bool tf )
{
    is_missing = tf;
}


void PoMoState::setPosition(size_t pos)
{
    position = pos;
}



/**
 * Setting the PoMo state from  counts.
 *
 * Example with only ten states and DNA:
 * A C G T A10C90 A20C80 A30C70...A90C10 A10G90 A20G80...A10T90...C10G90...C10T90...G10T90
 * The preferred format is that of counts: e.g.:
 *    0,1,4,0
 * meaning 0 A, 1 C, 4 G, 0 T were sampled at that position.
 */
void PoMoState::setState(const std::string &symbol)
{
    // store for internal use
    string_value = symbol;
        
    std::vector<std::string> counts_string;
    StringUtilities::stringSplit(symbol, ",", counts_string);
    
    if (counts_string.size() != num_raw_states)
    {
        throw RbException( "PoMo string state not correctly formatted. We found "+ symbol +", but the preferred format is that of counts, e.g. 0,1,4,0 meaning 0 A, 1 C, 4 G, 0 T were sampled at that position." );
    }
    
    std::vector<size_t> counts = std::vector<size_t>(num_raw_states, 0);
    for (size_t i = 0; i<num_raw_states; ++i)
    {
        counts[i] = StringUtilities::asIntegerNumber( counts_string[i] );
    }
    
    setState( counts );
}



/**
 * Setting the PoMo state from  counts.
 *
 * Example with only ten states and DNA:
 * A C G T A10C90 A20C80 A30C70...A90C10 A10G90 A20G80...A10T90...C10G90...C10T90...G10T90
 * The preferred format is that of counts: e.g.:
 *    0,1,4,0
 * meaning 0 A, 1 C, 4 G, 0 T were sampled at that position.
 */
void PoMoState::setState(const std::vector<size_t> &counts)
{
    size_t virt_pop_size_minus_1 = virtual_population_size-1;
    
    if (counts.size() != num_raw_states)
    {
        throw RbException( "PoMo string state not correctly formatted. We expected " + StringUtilities::toString(num_raw_states) + " counts but received " + StringUtilities::toString(counts.size()) + "." );
    }
    
    // We have the counts, now we create the state.
    size_t total_num_samples        = 0;
    size_t count_observed_alleles   = 0;
    size_t index_first_allele       = -1;
    size_t index_second_allele      = -1;
    size_t count_first_allele       = 0;
    // Sum over elements and count non-zero elements.
    for (size_t i = 0; i < num_raw_states; ++i)
    {
        size_t allele_count = counts[i];
        if ( allele_count > 0 )
        {
            // determines the nucleotide or allele type.
            if (index_first_allele == -1)
            {
                index_first_allele = i;
                count_first_allele = allele_count;
            }
            else if (index_second_allele == -1)
            {
                index_second_allele = i;
            }
            else
            {
                throw RbException("We current only support biallelic states in the PoMo framework.");
            }
            ++count_observed_alleles;
            total_num_samples += allele_count;
        }
    }
    
    size_t state_index = 0;
    if ( count_observed_alleles == 1 ) // monoallelic state
    {
        state_index = index_first_allele;
        state.clear();
        state.set(index_first_allele);
//        populateWeightedStatesForMonoallelicState(index_first_allele, total_num_samples);
    }
    else if ( count_observed_alleles == 2 ) //biallelic state
    {
        size_t basic_index = num_raw_states * index_first_allele;
        for (size_t i=0; i<index_first_allele; ++i)
        {
            basic_index -= i;
        }
        basic_index *= virt_pop_size_minus_1;
        basic_index += num_raw_states;
        
        if ( index_second_allele > num_raw_states )
        {
            throw RbException( "PoMo string state not correct. The second allele could not be determined." );
        }
        
        state.clear();
        // index corresponds to the closest place where the pomo state is.
        // In case the virtual population size is smaller to the counts in the state, or the reverse,
        // we have to do some maths.
        // We have to get the closest numbers to the observed counts.
        // Basically, the lowest count has to be >=1, and as close as possible
        // to the observed count.
        if (total_num_samples != virtual_population_size)
        {
            double obs_proportion_second_allele = (double) (total_num_samples - count_first_allele) / (double)total_num_samples ;
            size_t virtual_pop_sample_count = (size_t)round((1.0-obs_proportion_second_allele)*virtual_population_size);
            if ( virtual_pop_sample_count == virtual_population_size )
            {
                state_index = index_first_allele;
            }
            else if ( virtual_pop_sample_count == 0 )
            {
                state_index = index_second_allele;
            }
            else
            {
                state_index = basic_index + virtual_population_size - virtual_pop_sample_count - 1;
            }
        }
        else
        {
            state_index = basic_index + virtual_population_size - count_first_allele - 1;
        }
        state.set(state_index);
        
        
//        // Let's try PoMo state averaging.
//        // Basically all cells in the weight matrix that contain only id1, only id2, or a combination of both should have a non-zero weight.
//        double n = (double)total_num_samples;
//        double p = (double)count_first_allele/(double)total_num_samples;
//
//        std::vector<double> prob (virtual_population_size );
//        for (size_t j =0; j <= virt_pop_size_minus_1; ++j)
//        {
//            prob[j] = RbMath::choose(total_num_samples, count_first_allele) *  pow( ( (double)j/double(virtual_population_size)) , count_first_allele) * pow( ( (double)(virtual_population_size - j)/double(virtual_population_size) ) , (double)(total_num_samples-count_first_allele)) ;
//            if (prob[j] < 1e-8)
//            {
//                prob[j] = 0.0;
//            }
//            //RbStatistics::Binomial::pdf(n, (double)j/double(virtual_population_size), (double)(virtual_population_size));
//        }
//
//        for (size_t j=0; j < virt_pop_size_minus_1; ++j)
//        {
//            weights[j+basic_index] = prob[j+1];
//            state.set( j+num_raw_states );
//        }
        
    }
    else if ( count_observed_alleles == 0 )
    {
        setGapState( true ); // We say we have a gap
    }
    
    index_single_state = state_index;
    num_observed_states = 1;
    
}


void PoMoState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.clear();
    state.set( index );
    
}


void PoMoState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.clear();
    state.set( 0 );
}


void PoMoState::setVirtualPopulationSize(size_t ps)
{
    if (ps >= 100)
    {
        throw RbException( "The virtual population size should be < 100 and should be a divisor of 100." );
    }
    if (100 % ps != 0)
    {
        throw RbException( "The virtual population size should be a divisor of 100." );
    }
    virtual_population_size = ps;
}
