#include "RateMatrix_PoMoBalanceKN.h"

#include <assert.h>
#include <cstddef>

#include "RbException.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with k alleles, a virtual population size, n states, a mutation matrix, fitnesses,
 * and balancing selection*/
RateMatrix_PoMoBalanceKN::RateMatrix_PoMoBalanceKN( long num_states, long in_k, long in_n, long in_nmr)  :
AbstractRateMatrix( num_states ),
K( in_k ),
N( in_n ),
mu( in_nmr , 0.01 ),
phi( in_k , 1.0 ),
beta( in_nmr*0.5, 1.0 ),
B( in_nmr*0.5, 1 )
{
  update();
}

/** Copy constructor */
RateMatrix_PoMoBalanceKN::RateMatrix_PoMoBalanceKN(const RateMatrix_PoMoBalanceKN& m) :
AbstractRateMatrix( m ),
K( m.K ),
N( m.N ),
mu( m.mu ),
phi( m.phi ),
beta( m.beta ),
B( m.B )
{

}


RateMatrix_PoMoBalanceKN& RateMatrix_PoMoBalanceKN::operator=(const RateMatrix_PoMoBalanceKN &r)
{

  if (this != &r)
  {
    AbstractRateMatrix::operator=( r );
    K                   = r.K;
    N                   = r.N;
    mu                  = r.mu;
    phi                 = r.phi;
    beta                = r.beta;
    B                   = r.B;
  }

  return *this;
}




/** Destructor */
RateMatrix_PoMoBalanceKN::~RateMatrix_PoMoBalanceKN(void)
{

}



RateMatrix_PoMoBalanceKN& RateMatrix_PoMoBalanceKN::assign(const Assignable &m)
{

    const RateMatrix_PoMoBalanceKN *rm = dynamic_cast<const RateMatrix_PoMoBalanceKN*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}


double RateMatrix_PoMoBalanceKN::averageRate(void) const
{
    return 1.0;
}





void RateMatrix_PoMoBalanceKN::buildRateMatrix(void)
{

 /*  
  INFORMATION ABOUT THE PoMoBalanceKN RATE MATRIX

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

  Compared to PoMoKN, we add the relative fitness phi and balancing selection defined with beta and B.
  */

   MatrixReal& m = *the_rate_matrix;

  //populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++){
    for (int j=0; j< num_states; j++){
      m[i][j] = 0.0;        
    }
  }
  //reciprocal of the population size
  double rN = 1.0 / N;
  //first edge
  int E = 0;

  //these for loops go through the (K*K-K) edges of the pomo state-space
  //their represent all the possible pairwise combinations of alleles: a0a1, a1a2, ..., aK-2aK-1
  for (int i=0; i<K; i++) {
      for (int j = i + 1; j < K; j++) {
          // mutations
          m[i][K + E * N - E] = mu[2 * E];                  //{Nai} -> {(N-1)ai,1aj}
          m[j][K + (E + 1) * (N - 1) - 1] = mu[2 * E + 1];  //{Naj} -> {1ai,(N-1)aj}

          // fixations
          m[K + E * N - E][i] = (N - 1.0) * phi[i] * rN;  //{(N-1)ai,1aj} -> {Nai}
          m[K + (E + 1) * (N - 1) - 1][j] = (N - 1.0) * phi[j] * rN;  //{1ai,(N-1)aj} -> {Naj}

          //the pomo rate matrix is entirely defined by fixations and mutations if N=2
          if (N > 2) {

              //frequency shifts from singletons
              m[K + E * N - E][K + E * N - E + 1] =
                      (N - 1.0) * phi[j] * pow(beta[E], 0.5 * (abs(N - 1 - B[E]) - abs(N - 2 - B[E]) + 1.0)) *
                      rN;             //{(N-1)ai,1aj} -> {(N-2)ai,2aj}
              m[K + (E + 1) * (N - 1) - 1][K + (E + 1) * (N - 1) - 2] =
                      (N - 1.0) * phi[i] * pow(beta[E], 0.5 * (abs(1 - B[E]) - abs(2 - B[E]) + 1.0)) *
                      rN;   //{1ai,(N-1)aj} -> {2ai,(N-2)aj}

              //frequency shifts for all the other polymorphic states
              if (N > 3) {

                  //polymorphic states are populated in two fronts, thus the need for the middle frequency (some of the values populated twice)
                  int S = N / 2 + 1;

                  for (int n = 2; n < S; n++) {

                      //populates the first half of the polymorphic edges
                      m[K + E * N - E + n - 1][K + E * N - E + n] = n * (N - n) * phi[j] * pow(beta[E],
                           0.5 *(abs(N - n - B[E]) - abs(N - n - 1 - B[E]) + 1.0)) * rN;
                      //{(N-n)ai,naj} -> {(N-n-1)ai,(n-1)aj}
                      m[K + E * N - E + n - 1][K + E * N - E + n - 2] = n * (N - n) * phi[i] * pow(beta[E],
                           0.5 * (abs(N - n - B[E]) - abs(N - n + 1 - B[E]) + 1.0)) * rN;
                      //{(N-n)ai,naj} -> {(N-n+1)ai,(n+1)aj}

                      //populates the second half of the polymorphic edges
                      m[K + (E + 1) * N - E - n - 1][K + (E + 1) * N - E - n] =
                          n * (N - n) * phi[j] * pow(beta[E], 0.5 * (abs(n - B[E]) - abs(n - 1 - B[E]) + 1.0)) * rN;
                      //{nai,(N-n)aj} -> {(n-1)ai,(N-n-1)aj}
                      m[K + (E + 1) * N - E - n - 1][K + (E + 1) * N - E - n - 2] =
                          n * (N - n) * phi[i] * pow(beta[E], 0.5 * (abs(n - B[E]) - abs(n + 1 - B[E]) + 1.0)) * rN;
                      //{nai,(N-n)aj} -> {(n+1)ai,(N-n+1)aj}
                  }
              }
          }
          //update edge
          E += 1;

      }
  }

    // set the diagonal values and calculate the reciprocal of the average rate
    for (int j = 0; j < num_states; j++) {
        for (int k = 0; k < num_states; k++) {
            if (j != k) {
                m[j][j] = m[j][j] - m[j][k];
            }
        }
    }

    
  /*
  //I cannot get the stationary frequencies by any means. So, I am just skiping the matrix normalization.
  std::vector< double > stationaryFrequencies = getStationaryFrequencies();
  double rRate = 0.0;
  for (int j=0; j< num_states; j++){
    rRate += diagonal[j]*stationaryVector[j];   
    std::cout << "rate:" << stationaryVector[j] << "\n";
  }

  rRate = 1.0/rRate;
  */


}


/** Calculate the transition probabilities */
void RateMatrix_PoMoBalanceKN::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{

  // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
  // Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
  double t = rate * (startAge - endAge);
  exponentiateMatrixByScalingAndSquaring(t, P );

  return;
}


RateMatrix_PoMoBalanceKN* RateMatrix_PoMoBalanceKN::clone( void ) const
{
    return new RateMatrix_PoMoBalanceKN( *this );
}


// We comment out this function because the stationary frequencies in this case seem intractable, will need to review
// this in the future
//std::vector<double> RateMatrix_PoMoBalanceKN::getStationaryFrequencies( void ) const
//{
//  return stationaryVector;
//}

std::vector<double> RateMatrix_PoMoBalanceKN::getStationaryFrequencies( void ) const
{
    return calculateStationaryFrequencies();
}


void RateMatrix_PoMoBalanceKN::setK( long & na )
{
    K = na;

    // set flags
    needs_update = true;

}

void RateMatrix_PoMoBalanceKN::setN( long & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMoBalanceKN::setMu( const std::vector<double> &m )
{
    mu = m;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoBalanceKN::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoBalanceKN::setBeta( const std::vector<double> &b )
{
    beta = b;

    // set flags
    needs_update = true;
}

void RateMatrix_PoMoBalanceKN::setB( const std::vector<long> &Bf )
{
    B = Bf;

    // set flags
    needs_update = true;
}


void RateMatrix_PoMoBalanceKN::update( void )
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
