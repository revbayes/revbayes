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
#include "DistributionBinomial.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"

using namespace RevBayesCore;



/** Constructor that sets the observation and the other fields */
PoMoState::PoMoState(size_t n, size_t vps, const std::string &s, const std::string &chr, size_t pos) : DiscreteCharacterState( n + size_t(RbMath::kchoose2(int(n)))*(vps-1) ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    virtual_population_size( vps ),
    num_raw_states( n ),
    num_pomo_states( n + size_t(RbMath::kchoose2(int(n)))*(virtual_population_size-1) ),
    state( num_pomo_states ),
    chromosome( chr ),
    position( pos ),
    string_value(s)
{
    
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


/**
 * Compute the index from which this biallelic state will start.
 * Remember that the first k indices are used for the monomorphic raw states/
 * Then, we sort the indices by the first allele index in ascending order of the frequency of the first allele.
 *
 * Example with only ten states and DNA:
 * A C G T A10C90 A20C80 A30C70...A90C10 A10G90 A20G80...A10T90...C10G90...C10T90...G10T90
 */
size_t PoMoState::computeIndexBiallelic(size_t index_first_allele, size_t index_second_allele) const
{
    size_t basic_index = 0;
    
    // first, we move ahead the raw states
    basic_index += num_raw_states;
    
    // second, we need to move over all the biallelic states with a smaller first index (f_i)
    size_t first_allele_offset = num_raw_states * (num_raw_states-1) / 2;
    first_allele_offset -= (num_raw_states-index_first_allele) * (num_raw_states-index_first_allele-1) / 2;
    basic_index += first_allele_offset * (virtual_population_size-1);
    
    // third, we need to move over all biallelic states with a smaller second index
    basic_index += (index_second_allele - index_first_allele - 1) * (virtual_population_size-1);
    
    return basic_index;
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


bool PoMoState::isStateIncludedInAscertainmentBiasCorrection( void ) const
{    
    size_t index = getStateIndex();
    return index < num_raw_states;
}


const std::string& PoMoState::nexusSeparator(void) const
{
    static std::string sep = " ";
    
    return sep;
}


//void PoMoState::setStateBinomialForMonomorphic(size_t total_samples, size_t index_first_allele)
//{
//    // We do PoMo state averaging.
//    double n = (double)total_samples;
//    double p = 1.0;
//
//    weights[index_first_allele] = 1.0;
//    state.set( index_first_allele );
//
//    int vps_minus_1 = virtual_population_size - 1;
//
//
//    std::vector<double> nd (vps_minus_1, 0.0);
//    std::vector<double> id (vps_minus_1, 0.0);
//    for (size_t i = 0; i < vps_minus_1; ++i)
//    {
//        nd[i] = (double)(i + 1) / (double)virtual_population_size;
//        id[i] = (double)(vps_minus_1 - i) / (double)virtual_population_size;
//    }
//
//    if (index_first_allele < num_raw_states )
//    {
//
//        for (size_t ind_second=0; ind_second<num_raw_states; ++ind_second)
//        {
//            if ( index_first_allele != ind_second )
//            {
//                size_t first  = (index_first_allele < ind_second ? index_first_allele : ind_second);
//                size_t second = (index_first_allele > ind_second ? index_first_allele : ind_second);
//
//                size_t index = num_raw_states*first;
//                for (size_t i=0; i<first; ++i)
//                {
//                    index_first_allele -= i;
//                }
//                index *= vps_minus_1;
//                index += num_raw_states;
//
//                for (size_t offset = 0; offset< vps_minus_1; ++offset)
//                {
//                    weights[ index + offset ] = pow(nd[offset], (double)total_samples);//RbStatistics::Binomial::pdf(n, p, (double)(numid1));
//                    state.set( index + offset );
//                }
//            }
//        }
//
//    }
//    else
//    {
//        throw RbException( "PoMo string state not correct. We found "+ StringUtilities::to_string(index_first_allele)  );
//    }
//
//    for (size_t i =0; i < weights.size(); ++i)
//    {
//        if (weights[i] < 1e-8)
//        {
//            weights[i] = 1e-8;
//        }
//    }
//
//}


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
 * Setting the PoMo state from counts string.
 * This function converts the counts string into a counts vector and delegates the transformation.
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
    
    if ( symbol == "-" || symbol == "?" )
    {
        setMissingState( true );
    }
    else
    {
        
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
}



/**
 * Setting the PoMo state from counts vector.
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
    size_t count_second_allele      = 0;
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
                count_second_allele = allele_count;
            }
            else
            {
                if ( count_first_allele < count_second_allele )
                {
                    if ( allele_count > count_first_allele )
                    {
                        index_first_allele = i;
                        count_first_allele = allele_count;
                    }
                }
                else
                {
                    if ( allele_count > count_second_allele )
                    {
                        index_second_allele = i;
                        count_second_allele = allele_count;
                    }
                    
                }
//                throw RbException("We current only support biallelic states in the PoMo framework.");
            }
            ++count_observed_alleles;
            total_num_samples += allele_count;
        }
    }
    
    size_t state_index = 0;
    if ( count_observed_alleles == 1 ) // monoallelic state
    {
        state.clear();
        
        if ( weighting == FIXED )
        {
            state_index = index_first_allele;
        }
        else if ( weighting == BINOMIAL )
        {
            setStateBinomialForMonomorphic(total_num_samples, index_first_allele);
        }
        else if ( weighting == SAMPLED )
        {
//            setStateSampled();
            state_index = index_first_allele;
        }
        else
        {
            throw RbException() << "Unknown weighting type in PoMo state.";
        }
    }
    else if ( count_observed_alleles == 2 ) // biallelic state
    {
        
        if ( index_second_allele >= num_raw_states )
        {
            throw RbException( "PoMo string state not correct. The second allele could not be determined." );
        }
        
        size_t basic_index = computeIndexBiallelic(index_first_allele, index_second_allele);
        
        state.clear();
        
        
        if ( weighting == FIXED )
        {
            setStateFixed(total_num_samples, count_first_allele, basic_index);
        }
        else if ( weighting == BINOMIAL )
        {
            setStateBinomial(total_num_samples, count_first_allele, basic_index);
        }
        else if ( weighting == SAMPLED )
        {
            setStateSampled(total_num_samples, count_first_allele, basic_index);
        }
        else
        {
            throw RbException() << "Unknown weighting type in PoMo state.";
        }
        
    }
    else if ( count_observed_alleles == 0 )
    {
        setGapState( true ); // We say we have a gap
    }
    
}

/**
 * Computing the PoMo state as weights from a binomial distribution.
 *
 */
void PoMoState::setStateBinomial(size_t total_num_samples, size_t count_first_allele, size_t basic_index)
{
    // clear the weights
    weights.clear();
    
    size_t virt_pop_size_minus_1 = virtual_population_size-1;
        
    // Let's try PoMo state averaging.
    // Basically all cells in the weight matrix that contain only id1, only id2, or a combination of both should have a non-zero weight.
    double n = (double)total_num_samples;

    for (size_t j=1; j <= virt_pop_size_minus_1; ++j)
    {
        double prob = RbStatistics::Binomial::pdf(n, double(j)/double(virtual_population_size), count_first_allele);
        if (prob < 1e-8)
        {
            prob = 0.0;
        }
        weights[j+basic_index-1] = prob;
        state.set( j+basic_index-1 );
    }
    
    index_single_state = -1;
    num_observed_states = 1;
    
}

/**
 * Computing the PoMo state as weights from a binomial distribution for an observed monomorphic states.
 * That means we could either have a monomorphic frequency, or
 * any possible biallelic state where the observe allele is part of the biallelic combination.
 *
 * For example, the true frequency could be A80C20 but we sampled 2 A, 0 C, 0 G, 0 T (all being the allele A).
 *
 */
void PoMoState::setStateBinomialForMonomorphic(size_t total_samples, size_t index_first_allele)
{
    
    // We do PoMo state averaging.
    double n = (double)total_samples;
    double p = 1.0;
    
    // we can directly set the weight for the monomorphic case
    weights[index_first_allele] = 1.0;
    state.set( index_first_allele );
    
    size_t virt_pop_size_minus_1 = virtual_population_size-1;

    // iterate over all possible second alleles
    for (size_t index_second_allele=0; index_second_allele<num_raw_states; ++index_second_allele)
    {
        
        // first check that we don't have the same alleles
        if ( index_first_allele == index_second_allele )
        {
            continue;
        }
        
        // now get allele indices in a sorted order
        size_t index_a = std::min(index_first_allele, index_second_allele);
        size_t index_b = std::max(index_first_allele, index_second_allele);
        
        // compute the basic index (offset) for this biallelic combination
        size_t basic_index = computeIndexBiallelic(index_a, index_b);
        
        // now iterate over all biallelic frequencies
        for (size_t j=1; j < virtual_population_size; ++j)
        {
            // compute the PoMo frequency
            // this depends on whether this observed allele has the lower index
            double freq = 0.0;
            if ( index_first_allele < index_second_allele )
            {
                freq = double(j) / virtual_population_size;
            }
            else
            {
                freq = 1.0 - double(j) / virtual_population_size;
            }
            
            // compute the probability of binomially sampling the monomorphic state
            double prob = pow(freq, n);
            
            // cap the probabilities (for numerical stability???)
            if (prob < 1e-8)
            {
                prob = 0.0;
            }
            
            // set the current weight
            weights[j+basic_index-1] = prob;
            state.set( j+basic_index-1 );
            
        } // end loop over all biallelic frequencies
        
    } // end loop over all biallelic combinations
    
    index_single_state = -1;
    num_observed_states = 1;
    
}



/**
 * Setting the PoMo state from counts as the most closely matching frequency.
 * Here we specifically try to match the frequency if possible by rounding but make sure that the internal PoMo state is not a monomorphic state.
 * Internal monomorphic PoMo state could never generate biallelic observed states.
 */
void PoMoState::setStateFixed(size_t total_num_samples, size_t count_first_allele, size_t basic_index)
{
    
    size_t state_index = 0;
    
    // index corresponds to the closest place where the PoMo state is.
    // We have to get the closest numbers to the observed counts.
    double obs_proportion_first_allele = std::fmax(1.0, double(count_first_allele) / double(total_num_samples) );
    size_t virtual_pop_sample_count = (size_t)round(obs_proportion_first_allele*virtual_population_size);
    // we want to make sure that we don't assume a monomorphic state due to mapping
    // because a monomorphic PoMo frequency can never produce a biallelic observed state
    if ( virtual_pop_sample_count == virtual_population_size )
    {
        --virtual_pop_sample_count;
    }
    else if ( virtual_pop_sample_count == 0 )
    {
        ++virtual_pop_sample_count;
    }
    else
    {
        state_index = basic_index + virtual_population_size - virtual_pop_sample_count - 1;
    }
       
    // set the internal values
    index_single_state = state_index;
    num_observed_states = 1;
    state.set(state_index);


}


/**
 * Computing the PoMo state as being sampled from a binomial distribution.
 * That is, we compute the weight as the binomial probability of sampling the number of first alleles given the PoMo frequencies with the total number of samples as the number of trials.
 * Then we randomly sample according to the weights.
 *
 */
void PoMoState::setStateSampled(size_t total_num_samples, size_t count_first_allele, size_t basic_index)
{
            
    size_t virt_pop_size_minus_1 = virtual_population_size-1;
        
    // Let's try PoMo state averaging.
    // Basically all cells in the weight matrix that contain only id1, only id2, or a combination of both should have a non-zero weight.
    double n = (double)total_num_samples;

    std::vector<double> prob (virtual_population_size );
    double max_prob = 0;
    for (size_t j=1; j <= virt_pop_size_minus_1; ++j)
    {
        prob[j] = RbStatistics::Binomial::pdf(n, (double)j/double(virtual_population_size), (double)(count_first_allele));
        if ( prob[j] > max_prob )
        {
            max_prob = prob[j];
        }
    }
    
    // normalize the probabilities
    for (size_t j=1; j <= virt_pop_size_minus_1; ++j)
    {
        prob[j] /= max_prob;
    }
    
    // sample the index
    size_t sample_index = 1;
    double u = GLOBAL_RNG->uniform01();
    while (sample_index < virtual_population_size && u > prob[sample_index])
    {
        u -= prob[sample_index];
        ++sample_index;
    }
    size_t state_index = basic_index + virtual_population_size - sample_index - 1;

    
    index_single_state = state_index;
    num_observed_states = 1;
    state.set(state_index);

}


void PoMoState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.reset();
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

const std::vector<double>& PoMoState::getWeights( void ) const
{
    return weights;
}

bool PoMoState::isWeighted( void ) const
{
    return weighted;
}

void PoMoState::setWeighted( bool tf )
{
    weighted = tf;
}
void PoMoState::setWeighting( PoMoState::WEIGHTING weight_type )
{
    weighting = weight_type;
}


