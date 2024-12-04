#include "PoMoState4.h"

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

/** Default constructor */
PoMoState4::PoMoState4(void) : DiscreteCharacterState( 4 + 6 * (10 - 1) ),
    is_gap( false ),
    is_missing( false ),
    index_single_state( 0 ),
    num_observed_states( 0 ),
    state(4 + 6 * (10 - 1)),
    chromosome_ (""),
    position_(0),
    virtualPopulationSize_ ( 10 ),
    stringValue_ ("")
{
    setWeighted(true);
    weights.assign(4 + 6 * (virtualPopulationSize_ - 1), 0.0);

}

/** Constructor with virtual population size */
PoMoState4::PoMoState4(size_t vps): DiscreteCharacterState( 4 + 6 * (vps - 1) ),
is_gap( false ),
is_missing( false ),
index_single_state( 0 ),
num_observed_states( 0 ),
state(4 + 6 * (10 - 1)),
chromosome_ (""),
position_(0),
virtualPopulationSize_( vps ),
stringValue_ ("")
{

    setWeighted(true);
    weights.assign(4 + 6 * (virtualPopulationSize_ - 1), 0.0);

}


/** Constructor that sets the observation */
PoMoState4::PoMoState4(const std::string &s) : DiscreteCharacterState( 4 + 6 * (10 - 1) ),
is_gap( false ),
is_missing( false ),
index_single_state( 0 ),
num_observed_states( 0 ),
state(4 + 6 * (10 - 1)),
chromosome_( "" ),
position_( 0 ),
virtualPopulationSize_( 10 ),
stringValue_ (s)
{

    //assert( s <= 15 );
    setWeighted(true);
    weights.assign(4 + 6 * (virtualPopulationSize_ - 1), 0.0);
    setState(s);
}



/** Constructor that sets the observation and the other fields */
PoMoState4::PoMoState4(const std::string &s, const std::string chromosome, const size_t position, const size_t virtualPopulationSize, std::vector<double> weights ) : DiscreteCharacterState( 4 + 6 * (virtualPopulationSize - 1) ),
is_gap( false ),
is_missing( false ),
index_single_state( 0 ),
num_observed_states( 0 ),
state(4 + 6 * (virtualPopulationSize - 1)),
chromosome_( chromosome ),
position_( position ),
virtualPopulationSize_( virtualPopulationSize ),
stringValue_ (s)
{

    setWeighted(true);
    weights.assign(4 + 6 * (virtualPopulationSize_ - 1), 0.0);
    setState(s);

}


/** Constructor that sets the observation and the other fields */
PoMoState4::PoMoState4(const std::string &s, const std::string chromosome, const size_t position, const size_t virtualPopulationSize) : DiscreteCharacterState( 4 + 6 * (virtualPopulationSize - 1) ),
is_gap( false ),
is_missing( false ),
index_single_state( 0 ),
num_observed_states( 0 ),
state(4 + 6 * (virtualPopulationSize - 1)),
chromosome_( chromosome ),
position_( position ),
virtualPopulationSize_( virtualPopulationSize ),
stringValue_ (s)
{

    setWeighted(true);
    weights.assign(4 + 6 * (virtualPopulationSize_ - 1), 0.0);
    setState(s);

}


/* Copy constructor */
PoMoState4::PoMoState4(const PoMoState4& t) : DiscreteCharacterState( 4 + 6 * (t.virtualPopulationSize_ - 1) ), is_gap(t.is_gap), is_missing(t.is_missing), index_single_state(t.index_single_state), num_observed_states(t.num_observed_states), state(t.state), chromosome_(t.chromosome_), position_(t.position_), virtualPopulationSize_(t.virtualPopulationSize_), stringValue_(t.stringValue_)
 {

   setWeighted ( t.weighted );
   weights.clear();
   for (size_t i = 0; i < t.weights.size() ; ++i) {
     weights.push_back(t.weights.at(i));
   }
   state = t.state;

 }



PoMoState4* PoMoState4::clone( void ) const
{

    return new PoMoState4( *this );
}


void PoMoState4::populateWeightedStatesForMonoallelicState(int id1, int sum) {
  // We do PoMo state averaging.
  double n = (double)sum;
  double p = 1.0;

  weights[id1] = 1.0;
  state.set( id1 );
  int virtMinus1 = virtualPopulationSize_ - 1;


  std::vector<double> nd (virtMinus1, 0.0);
  std::vector<double> id (virtMinus1, 0.0);
  for (size_t i = 0; i < virtMinus1; ++i) {
    nd[i] = (double)(i + 1) / (double)virtualPopulationSize_;
    id[i] = (double)(virtMinus1 - i) / (double)virtualPopulationSize_;
  }

  if (id1 == 0) {
    // A, AC, AG, AT weights will be non 0.
    for (size_t numid1 = 0; numid1< virtMinus1; ++numid1)
    {
      //A1C9 and subsequent
      weights[ 4 + numid1 ] = pow(nd[numid1], (double)sum);//RbStatistics::Binomial::pdf(n, p, (double)(numid1));
      state.set( 4 + numid1 );
      //A1G9 and subsequent
      weights[ 4 + virtMinus1 + numid1 ] = pow(nd[numid1], (double)sum);//RbStatistics::Binomial::pdf(n, p, (double)(numid1));
      state.set( 4 + virtMinus1 + numid1 );
      //A1T9 and subsequent
      weights[ 4 + 2 * virtMinus1 + numid1 ] = pow(nd[numid1], (double)sum);//RbStatistics::Binomial::pdf(n, p, (double)(numid1));
      state.set( 4 + 2 * virtMinus1 + numid1 );
    }
  }
  else if (id1 == 1) {
    // C, AC, CG, CT weights will be non 0.
    for (size_t numid1 = 0; numid1< virtMinus1; ++numid1)
    {
      //A.C.
      weights[ 4 + numid1 ] = pow(id[numid1], (double)sum);
      state.set( 4 + numid1 );
      //C.G.
      weights[ 4 + 3 * virtMinus1 + numid1 ] = pow(nd[numid1], (double)sum);
      state.set( 4 + 3 * virtMinus1 + numid1 );
      //C.T.
      weights[ 4 + 4 * virtMinus1 + numid1 ] = pow(nd[numid1], (double)sum);
      state.set( 4 + 4 * virtMinus1 + numid1 );
    }
  }
  else if (id1 == 2) {
    // G, AG, CG, GT weights will be non 0.
    for (size_t numid1 = 0; numid1< virtMinus1; ++numid1)
    {
      //A.G.
      weights[ 4 + virtMinus1 + numid1 ] = pow(id[numid1], (double)sum);
      state.set( 4 + virtMinus1 + numid1 );
      //C.G.
      weights[ 4 + 3 * virtMinus1 + numid1 ] = pow(id[numid1], (double)sum);
      state.set( 4 + 3 * virtMinus1 + numid1 );
      //G.T.
      weights[ 4 + 5 * virtMinus1 + numid1 ] = pow(nd[numid1], (double)sum);
      state.set( 4 + 5 * virtMinus1 + numid1 );
    }
  }
  else if (id1 == 3) {
    // T, AT, CT, GT weights will be non 0.
    for (size_t numid1 = 0; numid1< virtMinus1; ++numid1)
    {
      //A.T.
      weights[ 4 + 2 * virtMinus1 + numid1 ] = pow(id[numid1], (double)sum);
      state.set( 4 + 2 * virtMinus1 + numid1 );
      //C.T.
      weights[ 4 + 4 * virtMinus1 + numid1 ] = pow(id[numid1], (double)sum);
      state.set( 4 + 4 * virtMinus1 + numid1 );
      //G.T.
      weights[ 4 + 5 * virtMinus1 + numid1 ] = pow(id[numid1], (double)sum);
      state.set( 4 + 5 * virtMinus1 + numid1 );
    }
  }
  else  {
    throw RbException() << "PoMo string state not correct. We found " << id1;
  }

  for (size_t i =0; i < weights.size(); ++i) {
    if (weights[i] < 1e-8) {
      weights[i] = 1e-8;
    }
  //    std::cout << weights[i] << " ";
  }


  return;
}


void PoMoState4::setState(const std::string &symbol)
{
    /* Example with only ten states:
     A C G T A10C90 A20C80 A30C70...A90C10 A10G90 A20G80...A10T90...C10G90...C10T90...G10T90
     The preferred format is that of counts: e.g.:
     0,1,4,0
     meaning 0 A, 1 C, 4 G, 0 T were sampled at that position.
     */
    size_t index = 0;
    size_t virtMinus1 = virtualPopulationSize_-1;
    //Checking if we have the preferred format, i.e. counts.
    if ( symbol.find(",") != std::string::npos ) {
        std::stringstream ss(symbol);
        std::vector<int> vect;
        int j;
        while (ss >> j)
        {
            vect.push_back(j);
            if (ss.peek() == ',' || ss.peek() == ' ')
                ss.ignore();
        }
        if (vect.size() != 4)
            throw RbException() << "PoMo string state not correctly formatted. We found " << symbol << ", but the preferred format is that of counts, e.g. 0,1,4,0 meaning 0 A, 1 C, 4 G, 0 T were sampled at that position." ; 
        //We have the counts, now we create the state.
        //code borrowed in part from IQ-Tree
        int sum = 0;
        int count = 0;
        int id1 = -1;
        int id2 = -1;
        int num = 0;
        // Sum over elements and count non-zero elements.
        for (std::vector<int>::iterator i = vect.begin(); i != vect.end(); ++i) {
            // `i` is an iterator object that points to some
            // element of `vect`.
            if (*i != 0) {
                // `i - vect.begin()` ranges from 0 to 3 and
                // determines the nucleotide or allele type.
                if (id1 == -1)
                {
                  id1 = i - vect.begin();
                  num = *i;
                }
                else id2 = i - vect.begin();
                count++;
                sum += *i;
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
          int basicIndex = 4;
          if (id1==0) {
            if (id2 == 1) {
              basicIndex = 4;
            }
            else if (id2 == 2) {
              basicIndex = 4 + (virtMinus1);
            }
            else if (id2 == 3) {
              basicIndex = 4 + 2*(virtMinus1);
            }
            else  {
              throw RbException() << "PoMo string state not correct. We found " << symbol  ; 
            }
          }
          else if (id1==1) {
            if (id2 == 2) {
              basicIndex = 4 + 3*(virtMinus1);
            }
            else if (id2 == 3) {
              basicIndex = 4 + 4*(virtMinus1);
            }
            else  {
              throw RbException() << "PoMo string state not correct. We found " << symbol  ; 
            }
          }
          else if (id1==2) {
            if (id2 == 3) {
              basicIndex = 4 + 5*(virtMinus1);
            }
            else  {
              throw RbException() << "PoMo string state not correct. We found " << symbol  ; 
            }
          }
          else  {
            throw RbException() << "PoMo string state not correct. We found " << symbol  ; 
          }

          std::cout << "basicIndex: "<<basicIndex <<std::endl;


          size_t situation = (id1+1)*2+(id2+1);

          state.reset();
          // index corresponds to the closest place where the pomo state is.
          // In case the virtual population size is inferior to the counts in the state, or the reverse,
          // we have to do some maths.
          // We have to get the closest numbers to the observed counts.
          // Basically, the lowest count has to be >=1, and as close as possible
          // to the observed count.
          if (sum != virtualPopulationSize_) {
            double obs_proportion = (double) (sum - num) / (double)sum ;
            size_t corrected_num = 0;
            std::cout<< "obs_proportion*virtualPopulationSize_: " << obs_proportion*virtualPopulationSize_ <<std::endl;
            // A minima we decide there must be at least
            // one sample from the least frequent allele
            if  (obs_proportion*virtualPopulationSize_ <= 1.0) {
              corrected_num = virtualPopulationSize_ - 1;
            }
            else {
              corrected_num = virtualPopulationSize_ - (size_t)round(obs_proportion*virtualPopulationSize_);
            }
            std::cout << "corrected_num: " << corrected_num << " vs num: " << num << std::endl;
            index = basicIndex + virtualPopulationSize_ - corrected_num - 1;
          }
          else {
            index = basicIndex + virtualPopulationSize_ - num - 1;
          }
          std::cout << "index: "<<index <<std::endl;

          //index = 3 + id1*virtualPopulationSize_ + id2 ;
          // Let's try PoMo state averaging.
          // Basically all cells in the weight matrix that contain only id1, only id2, or a combination of both should have a non-zero weight.
          double n = (double)sum;
          double p = (double)num/(double)sum;

          std::vector<double> prob (virtualPopulationSize_ );
          for (size_t j =0; j <= virtMinus1; ++j) {
            prob[j] = RbMath::choose(sum, num) *  pow( ( (double)j/double(virtualPopulationSize_)) , num) * pow( ( (double)(virtualPopulationSize_ - j)/double(virtualPopulationSize_) ) , (double)(sum-num)) ;
            if (prob[j] < 1e-8) {
              prob[j] = 0.0;
            }
            //RbStatistics::Binomial::pdf(n, (double)j/double(virtualPopulationSize_), (double)(virtualPopulationSize_));
          }

          if (situation == 4) { //A.C.
            for (size_t j =0; j < virtMinus1; ++j) {
                weights[j+4] = prob[j+1];
                state.set( j+4 );
            }
          }
          else if (situation == 5) { //A.G.
            for (size_t j =0; j < virtMinus1; ++j) {
                weights[virtMinus1+j+4] = prob[j+1];
                state.set( virtMinus1+j+4 );
            }
          }
          else if (situation == 6) { //A.T.
            for (size_t j =0; j < virtMinus1; ++j) {
                weights[4 + 2 * virtMinus1 + j] = prob[j+1];
                state.set( 4 + 2 * virtMinus1 + j );
            }
          }
          else if (situation == 7) { //C.G.
            for (size_t j =0; j < virtMinus1; ++j) {
                weights[4 + 3 * virtMinus1 + j] = prob[j+1];
                state.set( 4 + 3 * virtMinus1 + j );
            }
          }
          else if (situation == 8) { //C.T.
            for (size_t j =0; j < virtMinus1; ++j) {
                weights[4 + 4 * virtMinus1 + j] = prob[j+1];
                state.set( 4 + 4 * virtMinus1 + j );
            }
          }
          else if (situation == 10) { //G.T.
            for (size_t j =0; j < virtMinus1; ++j) {
                weights[4 + 5 * virtMinus1 + j] = prob[j+1];
                std::cout << " at " << 4 + 5 * virtMinus1 + j << " we have : " << prob[j+1] <<std::endl;
                state.set( 4 + 5 * virtMinus1 + j );
            }
          }
          else {
            throw RbException() << "PoMo string state not correct. We found " << symbol  ; 
          }
        }
        else if ( count == 0 )
        {
            setGapState( true ); //We say we have a gap
        }
        else
        {
            throw RbException() << "PoMo string state not correct. We found " << symbol << ", but only 2 states can be non-0." ; 
        }
      }
    index_single_state = index;
    num_observed_states = 1;

    std::cout << "stringvalue: " << stringValue_ <<std::endl;
    for (size_t i =0; i < weights.size(); ++i) {
      std::cout << weights[i] << " ";
    }
    std::cout << "\n\n "<<std::endl;
}


std::string PoMoState4::getDataType( void ) const
{

    return "PoMo";
}


std::string PoMoState4::getStateLabels( void ) const
{

    static std::string labels = "A C G T ";
    std::string acgt( "ACGT" );
    std::vector< size_t > frequencies;
    int stepSize = 100 / virtualPopulationSize_;
    for (size_t i = 1; i < virtualPopulationSize_; ++i)
    {
        frequencies.push_back(i*stepSize);
    }
    for ( size_t k = 0; k < acgt.size(); ++k )
    {
        char ch = acgt[k];

        for ( size_t j = k+1; j < acgt.size(); ++j )
        {
            char ch2 = acgt[j];
            for (size_t i = 0; i < virtualPopulationSize_-1; ++i)
            {
                labels += ch + boost::lexical_cast<std::string>(frequencies[i]) + ch2 + boost::lexical_cast<std::string>(frequencies[virtualPopulationSize_ - 2 - i]) + " ";
            }
        }
    }

    return labels;
}

std::string PoMoState4::getStringValue(void) const
{

    if ( isMissingState() )
    {
        return "?";
    }

    if ( isGapState() )
    {
        return "-";
    }

    return stringValue_;
}

void PoMoState4::setVirtualPopulationSize(size_t populationSize)
{
    if (populationSize >= 100)
    {
        throw RbException( "The virtual population size should be < 100 and should be a divisor of 100." );
    }
    if (100 % populationSize != 0)
    {
        throw RbException( "The virtual population size should be a divisor of 100." );
    }
    virtualPopulationSize_ = populationSize;
}

void PoMoState4::setChromosome(std::string chromosome){
    chromosome_ = chromosome;
}

void PoMoState4::setPosition(size_t position){
    position_ = position;
}

const std::string PoMoState4::getChromosome( void )
{
    return chromosome_;
}

const size_t PoMoState4::getPosition( void )
{
    return position_;
}


void PoMoState4::addState(const std::string &symbol)
{
    ++num_observed_states;

    std::string labels = getStateLabels();
    size_t pos = labels.find(symbol);

    state.set( pos );
    index_single_state = pos;
}


RbBitSet PoMoState4::getState(void) const
{
    // we need to clear the bits first
//    RbBitSet bs = RbBitSet( state.size());
  //  bs.set(index_single_state);

    return state;
//    return bs;
}


bool PoMoState4::isGapState( void ) const
{
    return is_gap;
}


bool PoMoState4::isMissingState( void ) const
{
    return is_missing;
}


void PoMoState4::setGapState( bool tf )
{
    is_gap = tf;
}


void PoMoState4::setMissingState( bool tf )
{
    is_missing = tf;
}


void PoMoState4::setToFirstState(void)
{
    num_observed_states = 1;
    index_single_state = 0;
    state.reset();
    state.set( 0 );
}


void PoMoState4::setStateByIndex(size_t index)
{

    num_observed_states = 1;
    index_single_state = index;
    state.reset();
    state.set( index );
}

const std::vector<double>& PoMoState4::getWeights( void ) const
{
    return weights;
}

bool PoMoState4::isWeighted( void ) const
{
    return weighted;
}

void PoMoState4::setWeighted( bool tf )
{
    weighted = tf;
}
