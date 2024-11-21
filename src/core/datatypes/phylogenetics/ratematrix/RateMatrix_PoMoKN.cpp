#include "RateMatrix_PoMoKN.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "RbException.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/*
RateMatrix_PoMoKN::RateMatrix_PoMoKN(size_t num_states) : 
AbstractRateMatrix( num_states ), 
K( in_k ),
N( in_n ),
mu( in_nmr, 0.01 ),
phi( in_k, 1.0 )
{
    update();
}
*/

/** Construct rate matrix with n states, an exchangeability matrix, a simplex of equilibrium frequencies, and a virtual population size */
RateMatrix_PoMoKN::RateMatrix_PoMoKN(std::int64_t num_states, std::int64_t in_k, std::int64_t in_n, std::int64_t in_nmr)  : 
AbstractRateMatrix( num_states ), 
K( in_k ),
N( in_n ),
mu( in_nmr, 0.01 ),
phi( in_k, 1.0 )
{
  update();
}

/** Copy constructor */
RateMatrix_PoMoKN::RateMatrix_PoMoKN(const RateMatrix_PoMoKN& m) : 
AbstractRateMatrix( m ), 
K( m.K ),
N( m.N ),
mu( m.mu ),
phi( m.phi ) 
{

}


RateMatrix_PoMoKN& RateMatrix_PoMoKN::operator=(const RateMatrix_PoMoKN &r)
{

  if (this != &r)
  {
    AbstractRateMatrix::operator=( r );
    K                   = r.K;
    N                   = r.N;
    mu                  = r.mu;
    phi                 = r.phi;
  }

  return *this;
}




/** Destructor */
RateMatrix_PoMoKN::~RateMatrix_PoMoKN(void)
{

}



double RateMatrix_PoMoKN::averageRate(void) const
{
    return 1.0;
}





void RateMatrix_PoMoKN::buildRateMatrix(void)
{

 /*  
  INFORMATION ABOUT THE PoMoKN RATE MATRIX

  It includes both fixed and polymorphic sites, but only biallelic states are considered.

  The pomo rate matrices defined here first list the fixed states {Na0}, {Na1} ...,
  these occupying positions 0:(K-1), and then polymorphic states.

  K alleles comprise (K*K-K)/2 pairwise combinations of alleles.
  This is the number of edges in the pomo state-space. Each edge comprises N-1 polymorphic states, 
  each of which represents a state of allelic frequencies in the population (summing to N).
  Example for a random allele pair aiaj: {(N-1)ai,1aj}, {(N-2)ai,2aj}, ..., {1ai,(N-1)aj}

  The polymorphic edges are listed in the following order a0a1, a0a2, a0a3, ..., a(K-2)aK-1
  We say the a0a1 polymorphic states sit at edge 0. Thus, {(N-1)a0,1a1} sits at position K and 
  {1a0,(N-1)a1} sits N-2 positions further (i.e., K+N-2).
  More generally, the polymorphic states of the allele pair sitting at edge E occupy the positions
  [K+E*N-E]:[K+(E+1)*(N-1)-1].
  */

   MatrixReal& m = *the_rate_matrix;

  //populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++){
    for (int j=0; j< num_states; j++){
      m[i][j] = 0.0;        
    }
  }


  // set the diagonal vector
  std::vector< double > diagonal(num_states,0.0); 

  //first edge
  int E = 0;

  //these for loops go through the (K*K-K) edges of the pomo state-space
  //their represent all the possible pairwise combinations of alleles: a0a1, a1a2, ..., aK-2aK-1
  for (int i=0; i<K; i++){
    for (int j=i+1; j<K; j++){

      //mutations
      m[i][K+E*N-E]          = N*mu[2*E];    //{Nai} -> {(N-1)ai,1aj}
      m[j][K+(E+1)*(N-1)-1]  = N*mu[2*E+1];  //{Naj} -> {1ai,(N-1)aj}
      diagonal[i] += m[i][K+E*N-E];
      diagonal[j] += m[j][K+(E+1)*(N-1)-1] ;

      //fixations
      m[K+E*N-E]        [i]  = (N-1.0)*phi[i]/(phi[j] + (N-1.0)*phi[i]);  //{(N-1)ai,1aj} -> {Nai} 
      m[K+(E+1)*(N-1)-1][j]  = (N-1.0)*phi[j]/(phi[i] + (N-1.0)*phi[j]);  //{1ai,(N-1)aj} -> {Naj} 
      diagonal[K+E*N-E]         += m[K+E*N-E]        [i];
      diagonal[K+(E+1)*(N-1)-1] += m[K+(E+1)*(N-1)-1][j];

      //the pomo rate matrix is entirely defined by fixations and mutations if N=2
      if (N>2) {

        //frequency shifts from singletons
        m[K+E*N-E]        [K+E*N-E+1]       = (N-1.0)*phi[j]/(phi[j] + (N-1.0)*phi[i]);  //{(N-1)ai,1aj} -> {(N-2)ai,2aj}
        m[K+(E+1)*(N-1)-1][K+(E+1)*(N-1)-2] = (N-1.0)*phi[i]/(phi[i] + (N-1.0)*phi[j]);  //{1ai,(N-1)aj} -> {2ai,(N-2)aj}
        diagonal[K+E*N-E]         += m[K+E*N-E]        [K+E*N-E+1] ;
        diagonal[K+(E+1)*(N-1)-1] += m[K+(E+1)*(N-1)-1][K+(E+1)*(N-1)-2];

        //frequency shifts for all the other polymorphic states
        if (N>3) {

          //polymorphic states are populated
          for (int n=2; n<(N-1); n++){

            //populates the first half of the polymorphic edge aiaj
            m[K+E*N-E+n-1]    [K+E*N-E+n]         = n*(N-n)*phi[j]/(n*phi[j] + (N-n)*phi[i]); //{(N-n)ai,naj} -> {(N-n-1)ai,(n+1)aj}
            m[K+E*N-E+n-1]    [K+E*N-E+n-2]       = n*(N-n)*phi[i]/(n*phi[j] + (N-n)*phi[i]); //{(N-n)ai,naj} -> {(N-n+1)ai,(n-1)aj}
            diagonal[K+E*N-E+n-1] += m[K+E*N-E+n-1][K+E*N-E+n] + m[K+E*N-E+n-1][K+E*N-E+n-2];

          }

        }

      }

      //update edge
      E += 1; 

    }
  }

  // set the diagonal values and calculate the reciprocal of the average rate

  double rRate = 0.0;

  for (int j=0; j< num_states; j++){
    m[j][j] = -diagonal[j];
  }
  
    
  /*
  //I cannot make the stationary frequencies to work, so I will just skiping the matrix noralization
  //std::vector< double > stationaryFrequencies = getStationaryFrequencies();
  for (int j=0; j< num_states; j++){
    rRate += diagonal[j]*stationaryVector[j];   
    std::cout << "rate:" << stationaryVector[j] << "\n";
  }


  rRate = 1.0;///rRate;

  // rescaling the matrix
  //first edge
  E = 0;

  //these for loops go through the (K*K-K) edges of the pomo state-space
  //their represent all the possible pairwise combinations of alleles: a0a1, a1a2, ..., aK-2aK-1
  for (int i=0; i<K; i++){
    for (int j=i+1; j<K; j++){

      //mutations
      m[i][K+E*N-E]          *= rRate;  //{Nai} -> {(N-1)ai,1aj}
      m[j][K+(E+1)*(N-1)-1]  *= rRate;  //{Naj} -> {1ai,(N-1)aj}


      //fixations
      m[K+E*N-E]        [i]  *= rRate;  //{(N-1)ai,1aj} -> {Nai} 
      m[K+(E+1)*(N-1)-1][j]  *= rRate;  //{1ai,(N-1)aj} -> {Naj} 

      //the pomo rate matrix is entirely defined by fixations and mutations if N=2
      if (N>2) {

        //frequency shifts from singletons
        m[K+E*N-E]        [K+E*N-E+1]       *= rRate;  //{(N-1)ai,1aj} -> {(N-2)ai,2aj}
        m[K+(E+1)*(N-1)-1][K+(E+1)*(N-1)-2] *= rRate;  //{1ai,(N-1)aj} -> {2ai,(N-2)aj}

        //frequency shifts for all the other polymorphic states
        if (N>3) {

          //polymorphic states are populated
          for (int n=2; n<(N-1); n++){

            //populates the first half of the polymorphic edge aiaj
            m[K+E*N-E+n-1]    [K+E*N-E+n]         *= rRate; //{(N-n)ai,naj} -> {(N-n-1)ai,(n+1)aj}
            m[K+E*N-E+n-1]    [K+E*N-E+n-2]       *= rRate; //{(N-n)ai,naj} -> {(N-n+1)ai,(n-1)aj}

          }

        }

      }

      //update edge
      E += 1; 

    }
  }
  */


}


/** Calculate the transition probabilities */
void RateMatrix_PoMoKN::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{

  // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
  // Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
  double t = rate * (startAge - endAge);
  exponentiateMatrixByScalingAndSquaring(t, P );

  return;
}


RateMatrix_PoMoKN* RateMatrix_PoMoKN::clone( void ) const
{
    return new RateMatrix_PoMoKN( *this );
}


std::vector<double> RateMatrix_PoMoKN::getStationaryFrequencies( void ) const
{
  return stationaryVector;
}



void RateMatrix_PoMoKN::setK( std::int64_t & na )
{   
    K = na;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMoKN::setN( std::int64_t & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMoKN::setMu( const std::vector<double> &m )
{
    mu = m;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoKN::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}


void RateMatrix_PoMoKN::update( void )
{

    if ( needs_update )
    {

        buildRateMatrix();

        // rescale: not useful, same loglk.
        //rescaleToAverageRate( 1.0 );

        // clean flags
        needs_update = false;
    }
}
