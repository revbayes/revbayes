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



PoMoState* PoMoState::clone( void ) const
{
    
    return new PoMoState( *this );
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
    else  {
        throw RbException() << "PoMo string state not correct. We found " << ind_first;
    }
    
    for (size_t i =0; i < weights.size(); ++i)
    {
        if (weights[i] < 1e-8)
        {
            weights[i] = 1e-8;
        }
        //    std::cout << weights[i] << " ";
    }
    
    
    return;
}


void PoMoState::setState(const std::string &symbol)
{
    /* Example with only ten states:
     A C G T A10C90 A20C80 A30C70...A90C10 A10G90 A20G80...A10T90...C10G90...C10T90...G10T90
     The preferred format is that of counts: e.g.:
     0,1,4,0
     meaning 0 A, 1 C, 4 G, 0 T were sampled at that position.
     */
    size_t index = 0;
    size_t virtMinus1 = virtual_population_size-1;
    
    std::vector<std::string> counts;
    StringUtilities::stringSplit(symbol, ",", counts);
    
    if (counts.size() != num_raw_states)
        throw RbException() << "PoMo string state not correctly formatted. We found " << symbol << ", but the preferred format is that of counts, e.g. 0,1,4,0 meaning 0 A, 1 C, 4 G, 0 T were sampled at that position." ; 
    
    
    //    //Checking if we have the preferred format, i.e. counts.
    //    if ( symbol.find(",") != std::string::npos ) {
    //        std::stringstream ss(symbol);
    //        std::vector<int> vect;
    //        int j;
    //        while (ss >> j)
    //        {
    //            vect.push_back(j);
    //            if (ss.peek() == ',' || ss.peek() == ' ')
    //                ss.ignore();
    //        }
    //    }
    
    // We have the counts, now we create the state.
    // code borrowed in part from IQ-Tree
    int sum = 0;
    int count = 0;
    size_t id1 = -1;
    size_t id2 = -1;
    int num = 0;
    // Sum over elements and count non-zero elements.
    for (size_t i = 0; i < num_raw_states; ++i)
    {
        int allele_count = StringUtilities::asIntegerNumber( counts[i] );
        if ( allele_count > 0 )
        {
            // determines the nucleotide or allele type.
            if (id1 == -1)
            {
                id1 = i;
                num = allele_count;
            }
            else if (id2 == -1)
            {
                id2 = i;
            }
            else
            {
                throw RbException("We current only support biallelic states in the PoMo framework.");
            }
            count++;
            sum += allele_count;
        }
    }
    
    if ( count == 1 ) // monoallelic state
    {
        index = id1 ; //+ 1;
        state.reset();
        state.set(id1);
        populateWeightedStatesForMonoallelicState(id1, sum);
    }
    else if ( count == 2 ) //biallelic state
    {
        std::cout << "BIALLELIC STATE" <<std::endl;
        int basicIndex = num_raw_states * id1;
        for (size_t i=0; i<id1; ++i) {
            basicIndex -= i;
        }
        basicIndex *= virtMinus1;
        basicIndex += num_raw_states;
        
        if ( id2 > num_raw_states )
        {
            throw RbException() << "PoMo string state not correct. We found " << symbol  ; 
        }
        
        std::cout << "basicIndex: "<<basicIndex <<std::endl;
        
        state.reset();
        // index corresponds to the closest place where the pomo state is.
        // In case the virtual population size is inferior to the counts in the state, or the reverse,
        // we have to do some maths.
        // We have to get the closest numbers to the observed counts.
        // Basically, the lowest count has to be >=1, and as close as possible
        // to the observed count.
        if (sum != virtual_population_size)
        {
            double obs_proportion = (double) (sum - num) / (double)sum ;
            size_t corrected_num = 0;
            std::cout<< "obs_proportion*virtual_population_size: " << obs_proportion*virtual_population_size <<std::endl;
            // A minima we decide there must be at least
            // one sample from the least frequent allele
            if  (obs_proportion*virtual_population_size <= 1.0)
            {
                corrected_num = virtual_population_size - 1;
            }
            else
            {
                corrected_num = virtual_population_size - (size_t)round(obs_proportion*virtual_population_size);
            }
            std::cout << "corrected_num: " << corrected_num << " vs num: " << num << std::endl;
            index = basicIndex + virtual_population_size - corrected_num - 1;
        }
        else
        {
            index = basicIndex + virtual_population_size - num - 1;
        }
        
        std::cout << "index: "<<index <<std::endl;
        
        // Let's try PoMo state averaging.
        // Basically all cells in the weight matrix that contain only id1, only id2, or a combination of both should have a non-zero weight.
        double n = (double)sum;
        double p = (double)num/(double)sum;
        
        std::vector<double> prob (virtual_population_size );
        for (size_t j =0; j <= virtMinus1; ++j)
        {
            prob[j] = RbMath::choose(sum, num) *  pow( ( (double)j/double(virtual_population_size)) , num) * pow( ( (double)(virtual_population_size - j)/double(virtual_population_size) ) , (double)(sum-num)) ;
            if (prob[j] < 1e-8)
            {
                prob[j] = 0.0;
            }
            //RbStatistics::Binomial::pdf(n, (double)j/double(virtual_population_size), (double)(virtual_population_size));
        }
        
        for (size_t j=0; j < virtMinus1; ++j)
        {
            weights[j+basicIndex] = prob[j+1];
            state.set( j+num_raw_states );
        }
        
    }
    else if ( count == 0 )
    {
        setGapState( true ); //We say we have a gap
    }
    
    index_single_state = index;
    num_observed_states = 1;
    
    std::cout << "stringvalue: " << string_value <<std::endl;
    for (size_t i =0; i < weights.size(); ++i)
    {
        std::cout << weights[i] << " ";
    }
    std::cout << "\n\n "<<std::endl;
}


std::string PoMoState::getDataType( void ) const
{
    
    return "PoMo";
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

void PoMoState::setChromosome(std::string chr)
{
    chromosome = chr;
}

void PoMoState::setPosition(size_t pos)
{
    position = pos;
}

const std::string& PoMoState::getChromosome( void ) const
{
    return chromosome;
}

size_t PoMoState::getPosition( void ) const
{
    return position;
}


void PoMoState::addState(const std::string &symbol)
{
    ++num_observed_states;
    
    std::string labels = getStateLabels();
    size_t pos = labels.find(symbol);
    
    state.set( pos );
    index_single_state = pos;
}


RbBitSet PoMoState::getState(void) const
{
    // we need to clear the bits first
    //    RbBitSet bs = RbBitSet( state.size());
    //  bs.set(index_single_state);
    
    return state;
    //    return bs;
}


bool PoMoState::isGapState( void ) const
{
    return is_gap;
}


bool PoMoState::isMissingState( void ) const
{
    return is_missing;
}


void PoMoState::setGapState( bool tf )
{
    is_gap = tf;
}


void PoMoState::setMissingState( bool tf )
{
    is_missing = tf;
}


void PoMoState::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.reset();
    state.set( 0 );
}


void PoMoState::setStateByIndex(size_t index)
{
    
    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
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
