#include "PoMoState.h"

#include <boost/lexical_cast.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <math.h>
#include <iostream>
#include <cstddef>
#include <string>
#include <algorithm>

#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "StringUtilities.h"
#include "Cloneable.h"
#include "DistributionBinomial.h"
#include "RandomNumberGenerator.h"
#include "RandomNumberFactory.h"

using namespace RevBayesCore;



/** Constructor that sets the observation and the other fields */
PoMoState::PoMoState(   size_t n, 
                        size_t vps, 
                        const std::string &s, 
                        const std::string &chr, 
                        size_t pos,
                        WEIGHTING w , 
                        const long eps ) : DiscreteCharacterState( n + size_t(RbMath::kchoose2(int(n)))*(vps-1) ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    virtual_population_size( vps ),
    n_alleles( n ),
    n_pomo_states( n + size_t(RbMath::kchoose2(int(n)))*(virtual_population_size-1) ),
    state( n_pomo_states ),
    chromosome( chr ),
    position( pos ),
    string_value(s),
    weighting( w )
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
 * Then, we sort the indices by the first allele index in descending order of the frequency of the first allele.
 *
 * Example with only ten states and the four nucleotide basis:
 * A C G T A9C1 A8C2 ... A1C9 A9G1 A8G2 ... A9T1 
 */
size_t PoMoState::computeEdgeFirstState(size_t index_first_allele, size_t index_second_allele) const
{
    /*
    RUI
    The following code is an overcomplicated way of compute the edge first state
    I am commenting it out and add my code below.

    size_t edge_first_state = 0;
    
    // first, we move ahead the raw states
    edge_first_state += n_alleles;
    
    // second, we need to move over all the biallelic states with a smaller first index (f_i)
    size_t first_allele_offset = n_alleles * (n_alleles-1) / 2;
    first_allele_offset -= (n_alleles-index_first_allele) * (n_alleles-index_first_allele-1) / 2;
    edge_first_state += first_allele_offset * (virtual_population_size-1);
    
    // third, we need to move over all biallelic states with a smaller second index
    edge_first_state += (index_second_allele - index_first_allele - 1) * (virtual_population_size-1);
    */

    size_t edge             = n_alleles*index_first_allele - ((index_first_allele+2)*(index_first_allele+2)-(index_first_allele+2))/2 + index_second_allele;
    size_t edge_first_state = n_alleles+edge*(virtual_population_size-1);

    return edge_first_state;
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


/*
SEBASTIAN CHECK:
I am not sure about the utility of this function.
I think PoMoState is not implemented for a general K.
*/
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
    return index < n_alleles;
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
//    if (index_first_allele < n_alleles )
//    {
//
//        for (size_t ind_second=0; ind_second<n_alleles; ++ind_second)
//        {
//            if ( index_first_allele != ind_second )
//            {
//                size_t first  = (index_first_allele < ind_second ? index_first_allele : ind_second);
//                size_t second = (index_first_allele > ind_second ? index_first_allele : ind_second);
//
//                size_t index = n_alleles*first;
//                for (size_t i=0; i<first; ++i)
//                {
//                    index_first_allele -= i;
//                }
//                index *= vps_minus_1;
//                index += n_alleles;
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
    
        if (counts_string.size() != n_alleles)
        {
            throw RbException( "PoMo string state not correctly formatted. We found "+ symbol +", but the preferred format is that of counts, e.g. 0,1,4,0 meaning 0 A, 1 C, 4 G, 0 T were sampled at that position." );
        }
    
        std::vector<size_t> counts = std::vector<size_t>(n_alleles, 0);
        for (size_t i = 0; i<n_alleles; ++i)
        {
            counts[i] = StringUtilities::asIntegerNumber( counts_string[i] );
        }
    
        setState( counts );
    }
}



/*
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
    
    if (counts.size() != n_alleles)
    {
        throw RbException( "PoMo string state not correctly formatted. We expected " + StringUtilities::toString(n_alleles) + " counts but received " + StringUtilities::toString(counts.size()) + "." );
    }
    

    /*
    // We have the counts, now we create the state.
    size_t total_count        = 0;
    size_t n_observed_alleles   = 0;
    size_t index_first_allele       = -1;
    size_t index_second_allele      = -1;
    size_t count_first_allele       = 0;
    size_t count_second_allele      = 0;

    // Sum over elements and count non-zero elements.
    for (size_t i = 0; i < n_alleles; ++i)
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
            ++n_observed_alleles;
            total_count += allele_count;
        }
    }

    RUI
    This code can be dangerous because admits triallelic sites by masking the allele with the lowest frequency
    while this seems a reasonable thing to do (e.g., sequencing errors are usually the lowest counts), 
    the user should decide what to do with these counts
    */

    // getting the counts and indexes for the observed alleles
    // throw error when more than 2 alleles are observed
    size_t index_first_allele, index_second_allele , count_first_allele, count_second_allele;

    size_t total_count        = 0;
    size_t n_observed_alleles = 0;

    for (size_t i = 0; i < n_alleles; ++i) 
    {
        size_t allele_count = counts[i];
        if (allele_count > 0) 
        {
            if (n_observed_alleles == 0) 
            {
                index_first_allele = i;
                count_first_allele = allele_count;
            } 
            else if (n_observed_alleles == 1) 
            {
                index_second_allele = i;
                count_second_allele = allele_count;
            } 
            else 
            {
                throw RbException("A 3- or 4-allelic count was found. We current only support biallelic counts in the PoMo framework.");
            }
            total_count += allele_count;
            n_observed_alleles++;
        }
    }


    size_t state_index = 0;

    // The observed count is monomorphic
    if ( n_observed_alleles == 1 ) 
    {
        //state.clear();
        
        if ( weighting == FIXED )
        {   
            state.reset();
            state_index = index_first_allele;
            state.set(state_index);
            std::cout << "  Index: " << state_index << "\n";


        }
        else if ( weighting == BINOMIAL )
        {
            /*
            SEBASTIAN CHECK:
            I was having and segentation fault error because of the weights vector. 
            Added this and it stoped, but not sure this is correct. 
            */
            setWeighted(true);
            weights.assign( n_pomo_states , 0.0);

            setStateBinomialForMonomorphic(total_count, index_first_allele);

        }
        else if ( weighting == SAMPLED )
        {
            //setStateSampled();
            /*
            SEBASTIAN CHECK:
            This is not the original sampled method. Monomorphic observed counts should still be able to produce polymorphic counts.
            This is how we recue variation that is not observed due to the sampling process.
            This should be simply sampling a state from the binomial weights directly.
            */
            state.reset();
            state_index = index_first_allele;
            state.set(state_index);
            std::cout << "  Index: " << state_index << "\n";

        }
        else if ( weighting == HYPERGEOMETRIC )
        {
             /*
            SEBASTIAN CHECK:
            I was having and segentation fault error because of the weights vector. 
            Added this and it stoped, but not sure this is correct. 
            */
            setWeighted(true);
            weights.assign( n_pomo_states , 0.0);
            setStateHypergeometricForMonomorphic(total_count, index_first_allele);
        }
        else
        {
            throw RbException() << "Unknown weighting method in PoMo state.";
        }
    }

    // The observed count is polymorphic and biallelic 
    else if ( n_observed_alleles == 2 ) 
    {
        
        if ( index_second_allele >= n_alleles )
        {
            throw RbException( "PoMo string state not correct. The second allele could not be determined." );
        }
        
        size_t edge_first_state = computeEdgeFirstState(index_first_allele, index_second_allele);
        
        if ( weighting == FIXED )
        {
             /*
            SEBASTIAN CHECK:
            Rev was claiming 0 states were observed by the time we clamp
            Adding state.reset() solved the problem, please make sure this makes sense
            borrow it from PoMoState4
            */
            state.reset();
            setStateFixed(total_count, count_first_allele, edge_first_state);
        }
        else if ( weighting == BINOMIAL )
        {
             /*
            SEBASTIAN CHECK:
            I was having and segentation fault error because of the weights vector. 
            Added this and it stoped, but not sure this is correct. 
            */
            setWeighted(true);
            weights.assign( n_pomo_states , 0.0);
            setStateBinomialForPolymorphic(total_count, count_first_allele, edge_first_state);
        }
        else if ( weighting == SAMPLED )
        {
            setStateSampled(total_count, count_first_allele, edge_first_state);
        }
        else if ( weighting == HYPERGEOMETRIC )
        {
             /*
            SEBASTIAN CHECK:
            I was having and segentation fault error because of the weights vector. 
            Added this and it stoped, but not sure this is correct. 
            */
            setWeighted(true);
            weights.assign( n_pomo_states , 0.0);
            setStateHypergeometricForPolymorphic(total_count, count_first_allele, edge_first_state);
        }
        else
        {
            throw RbException() << "Unknown weighting method in PoMo state.";
        }
        
    }
    // The count is empty
    else if ( n_observed_alleles == 0 )
    {
        setGapState( true ); // We say we have a gap
    }
    
}

/**
 * Computing the PoMo state as weights from a binomial distribution.
 *
 */
void PoMoState::setStateBinomialForPolymorphic(size_t total_count, size_t count_first_allele, size_t edge_first_state)
{
    // clear the weights
    //weights.clear();
            
    // PoMo state averaging.
    // Basically all cells in the weight matrix that contain a combination of both should have a non-zero weight.
    double n = (double)total_count;

    for (size_t j=1; j < virtual_population_size; ++j)
    {
        /*
        Note: The allele that increases in frequency in the edge is the second allele. 
        See description of PoMoKN. 
        */
        double prob = RbStatistics::Binomial::pdf(n, double(j)/double(virtual_population_size), total_count - count_first_allele);
        if (prob < 1e-8)
        {
            prob = 0.0;
        }
        weights[j+edge_first_state-1] = prob;
        state.set( j+edge_first_state-1 );
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
    for (size_t index_second_allele=0; index_second_allele<n_alleles; ++index_second_allele)
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
        size_t edge_first_state = computeEdgeFirstState(index_a, index_b);
        
        // now iterate over all biallelic frequencies
        for (size_t j=1; j < virtual_population_size; ++j)
        {
            // compute the PoMo frequency
            // this depends on whether this observed allele has the lower index
            double freq = 0.0;
            if ( index_first_allele < index_second_allele )
            {
                freq = 1.0 - double(j) / virtual_population_size;
            }
            else
            {
                freq = double(j) / virtual_population_size;
            }
            
            // compute the probability of binomially sampling the monomorphic state
            double prob = pow(freq, n);
            
            // cap the probabilities (for numerical stability???)
            if (prob < 1e-8)
            {
                prob = 0.0;
            }
            
            // set the current weight
            weights[j+edge_first_state-1] = prob;
            state.set( j+edge_first_state-1 );
            
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
void PoMoState::setStateFixed(size_t total_count, size_t count_first_allele, size_t edge_first_state)
{
    
    size_t state_index = 0;
    
    // index corresponds to the closest place where the PoMo state is.
    // We have to get the closest numbers to the observed counts.
    
    /*
    RUI 
    I corrected the indexing here
    The first allele descreases in frequency with the edge state indexing. See PoMoKN.

    We could as well have inputed directly the count of the second allele, which we already have
    double obs_proportion_first_allele = std::fmax(1.0, double(count_second_allele) / double(total_count) );
    */

    //double obs_proportion_first_allele = std::fmax(1.0, double(count_first_allele) / double(total_count) );
    double obs_proportion_second_allele = double(total_count-count_first_allele) / double(total_count);
    size_t virtual_pop_sample_count    = (size_t)round(obs_proportion_second_allele*virtual_population_size);

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
    
    state_index = edge_first_state + virtual_pop_sample_count - 1;
    
    // set the internal values
    index_single_state = state_index;
    num_observed_states = 1;
    state.set(state_index);
    std::cout << "  Index: " << state_index << "\n";

}


/**
 * Computing the PoMo state as being sampled from a binomial distribution.
 * That is, we compute the weight as the binomial probability of sampling the number of first alleles given the PoMo frequencies with the total number of samples as the number of trials.
 * Then we randomly sample according to the weights.
 *
 */
void PoMoState::setStateSampled(size_t total_count, size_t count_first_allele, size_t edge_first_state)
{
            
    size_t virt_pop_size_minus_1 = virtual_population_size-1;
        
    // Let's try PoMo state averaging.
    // Basically all cells in the weight matrix that contain only id1, only id2, or a combination of both should have a non-zero weight.
    double n = (double)total_count;

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
    size_t state_index = edge_first_state + virtual_population_size - sample_index - 1;


    state.reset();
    index_single_state = state_index;
    num_observed_states = 1;
    state.set(state_index);

    std::cout << "  Index: " << state_index << "\n";

}


/*
 * Computing the PoMo state as weights from the hypergeometric distribution given an effective population size of N
*/
void PoMoState::setStateHypergeometricForPolymorphic(size_t total_count, size_t count_first_allele, size_t edge_first_state)
{
    // clear the weights
    //weights.clear();

    // some inportant variables
    long   N = 10000;  // this will eventually be an input parameter with 10000 as default value
    size_t M = virtual_population_size;
    long   C = total_count;
    long   c = C - count_first_allele;  // this makes sence because edge states increase with aj and not ai

    // normalization constant
    // p(m!c)=p(c!m)*p(m)/p(c) where p(c!m) is a hypergeometric sampling; the other are discrete uniform
    // nc = p(m)/p(c)
    double ns = n_pomo_states;


    // weighting the polymorphic states

    // ni0 and nf0 are the effecive frequencies we want to sum up
    // there are V-1 ranges, where V is the virtual population size
    // ni and nf is just to guaretee that we do not sum off the reasonable ranfe of effective frequencies
    // i.e.: ni>=c and nf <= N-(C-c)
    long ni, nf, ni0, nf0;

    // some imporant quantities
    double hyper1, hyper2, weights_sum;
    //double total_sum = 0.0;

    // loops over the virtual frequencies
    for (int m=0; m<(M-1); ++m){

        // creates the ni:nf interval for a given m frequency
        ni0 = std::floor(1.0*m*(N-1)/(M-1)) + 1;
        nf0 = std::floor(1.0*(m+1)*(N-1)/(M-1)) + 1;
        ni  = std::max(c,ni0);
        nf  = std::min(N-C+c,nf0);

        // sums up the hypergeometric probabilities in the range ni:nf
        /*
        for (int n=ni; n<nf; ++n){
            weights_matrix(M-m-1,i) += hypergeometric_pdf(c, C, n, N);
        }
        
        This is the orginal method but was far too slow.
        Apparently it is very time consuming to call the hypergeometric probability in boost.
        
        I built an alternative method, which is not thouroughly tested but 
        is several orders of magnitude faster than using the hypergeometric of boost.
        */
        hyper1      = hypergeometric_pdf(c, C, ni, N);
        weights_sum = hyper1;

        for (int n=(ni+1); n<nf; ++n){

            // this is to avoid calling the hypergeometric function too many times
            // the racio is DETERMINED based on the properties of binomial coefficients
            hyper2       = hyper1*(1.0*n)*(N-n-C+c+1.0)/((1.0*n-c)*(N-n+1.0));
            weights_sum += hyper2;
            hyper1       = hyper2;
        }

        // the normalization constant for the polymorphic states has to acount for the fact 
        // that several states (nf0-ni0) are being summed up
        weights_sum /= (ns*(nf0-ni0));

        weights[edge_first_state + m] = weights_sum;
        state.set( edge_first_state+m );

        //total_sum += weights_sum;

        //if (weights_matrix(M-m-1,i) < 1.0e-8){
        //    weights_matrix(M-m-1,i) = 1.0e-8;
        //}

    }

    // normalizing the weights
    // I have tried without normalization and let to the same results
    // I don't think this step is necessary
    //for (int m=0; i < (M-1); i++)
    //{
    //    weights[edge_first_state + m] /= total_sum;
    //    std::cout << weights[i] << ' ';
    //}

    // still not sure why this is important
    index_single_state = -1;
    num_observed_states = 1;
}

void PoMoState::setStateHypergeometricForMonomorphic(size_t total_samples, size_t index_first_allele)
{
    // clear the weights
    //weights.clear();

    // some inportant variables
    long   N = 10000;  // this will eventually be an input parameter with 10000 as default value
    size_t M = virtual_population_size;
    long   C = total_samples;
    size_t K = n_alleles;

    // normalizytion constant
    // p(m!c)=p(c!m)*p(m)/p(c) where p(c!m) is a hypergeometric sampling; the other are discrete uniform
    // nc = p(m)/p(c)
    double ns = n_pomo_states;

    // weighting the monomorphic state 
    weights[index_first_allele] = hypergeometric_pdf(C, C, N, N)/ns; // this is 1.0
    state.set( index_first_allele );
    //double total_sum = weights[index_first_allele];

    // weighting the polymorphic states 

    // let us first determine all the alleles interacting with index_first_allele
    // we will be weighting polymorphic states concurently taking adavantage of the symetry of the hypergeometric sampling
    // for that we need the index of the initial or last edge state and the direction of weighting
    std::vector<size_t> edge_state(K-1,0);
    std::vector<int> direction(K-1,1);
    size_t edge;
    size_t index = 0;

    // all the interacting alleles with index smaller than index_first_allele
    // the weighting is in terms of the allele with larger index thus in this case
    // we will start at the first polymoprhic state and weight forwards (direction is positive)
    for (int i=0; i<index_first_allele; ++i) 
    {
        edge = K*i - ((i+2)*(i+2)-(i+2))/2 + index_first_allele;
        edge_state[index] = K+edge*(M-1);
        index += 1;
    }

    // all the interacting alleles with index larger than index_first_allele
    // the weighting is in terms of the allele with larger index thus in this case
    // we will start at the last polymoprhic state and weight backwards (direction is negative)    
    for (int i=(index_first_allele+1); i<K; ++i) 
    {
        edge = K*index_first_allele - ((index_first_allele+2)*(index_first_allele+2)-(index_first_allele+2))/2 + i;
        edge_state[index] = K+(edge+1)*(M-1)-1;
        direction[index] = -1;
        index += 1;
    }


    // ni0 and nf0 are the effecive frequencies we want to sum up
    // there are V-1 ranges, where V is the virtual population size
    // ni and nf is just to guaretee that we do not sum off the reasonable ranfe of effective frequencies
    // i.e.: ni>=c and nf <= N-(C-c)
    long ni, nf, ni0, nf0;

    // some imporant quantities
    double hyper1, hyper2, weights_sum;

    for (int m=0; m<(M-1); ++m)
    {
        // creates the ni:nf interval for a given m frequency
        ni0 = std::floor(1.0*m*(N-1)/(M-1)) + 1;
        nf0 = std::floor(1.0*(m+1)*(N-1)/(M-1)) + 1;
        ni  = std::max(C,ni0);
        nf  = std::min(N,nf0);

        // sums up the hypergeometric probabilities in the range ni:nf
        /*
        for (int n=ni; n<nf; ++n){
            weights_matrix(M-m-1,i) += hypergeometric_pdf(c, C, n, N);
        }
        
        This is the orginal method but was far too slow.
        Apparently it is very time consuming to call the hypergeometric probability in boost.
        
        I built an alternative method, which is not thouroughly tested but 
        is several orders of magnitude faster than using the hypergeometric of boost.
        */
        hyper1      = hypergeometric_pdf(C, C, ni, N);
        weights_sum = hyper1;

        for (int n=(ni+1); n<nf; ++n)
        {
            // this is to avoid calling the hypergeometric function too many times
            // based on the properties of binomial coefficients
            hyper2 = hyper1*(1.0*n)*(N-n+1.0)/((1.0*n-C)*(N-n+1.0));
            weights_sum += hyper2;
            //std::cout << "  " << hyper1 << "," << hyper2 << "\n";
            hyper1 = hyper2;
        }

        // the normalizytion constant for the polymorphic states has to acount for the fact 
        // that several states (nf-ni) are being summed up
        weights_sum /= (ns*(nf0-ni0));

        // updating total sum
        //total_sum += (K-1)*weights_sum;

        // updating the weiths for all the edges 
        for (int i=0; i<(K-1); ++i) 
        {
            weights[edge_state[i] + m*direction[i]] = weights_sum;
            state.set( edge_state[i] + m*direction[i] );
        }

        //if (weights_matrix(M-m-1,i) < 1.0e-8){
        //    weights_matrix(M-m-1,i) = 1.0e-8;
        //}

    }

    index_single_state = -1;
    num_observed_states = 1;

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
    /* 
    SEBASTIAN CHECK:
    I dont see the rational fot the following if statement. I am commenting it out. 
    
    if (ps >= 100)
    {
        throw RbException( "The virtual population size should be < 100 and should be a divisor of 100." );
    }
    if (100 % ps != 0)
    {
        throw RbException( "The virtual population size should be a divisor of 100." );
    }
    */
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


double PoMoState::hypergeometric_pdf(int c, int C, int n, int N)
{
    if (c > n || C-c > N-n){
        return 0.0;
    } else {
        boost::math::hypergeometric_distribution<double> hg_dist(n, C, N);
        return boost::math::pdf<double>(hg_dist, c);  
    }
}