#include "RateMatrix_PoMo2N.h"

#include <cassert>
#include <cstddef>

#include "RbException.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/*
RateMatrix_PoMo2N::RateMatrix_PoMo2N(size_t num_states) : 
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
RateMatrix_PoMo2N::RateMatrix_PoMo2N(long num_states, long in_n )  : 
AbstractRateMatrix( num_states ), 
N( in_n ),
mu( 2 , 0.01 ),
phi( 2 , 1.0 )
{
  update();
}

/** Copy constructor */
RateMatrix_PoMo2N::RateMatrix_PoMo2N(const RateMatrix_PoMo2N& m) : 
AbstractRateMatrix( m ), 
N( m.N ),
mu( m.mu ),
phi( m.phi )
{

}


RateMatrix_PoMo2N& RateMatrix_PoMo2N::operator=(const RateMatrix_PoMo2N &r)
{

  if (this != &r)
  {
    AbstractRateMatrix::operator=( r );
    N                   = r.N;
    mu                  = r.mu;
    phi                 = r.phi;
  }

  return *this;
}




/** Destructor */
RateMatrix_PoMo2N::~RateMatrix_PoMo2N(void)
{

}



RateMatrix_PoMo2N& RateMatrix_PoMo2N::assign(const Assignable &m)
{

    const RateMatrix_PoMo2N *rm = dynamic_cast<const RateMatrix_PoMo2N*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}


double RateMatrix_PoMo2N::averageRate(void) const
{
    return 1.0;
}





void RateMatrix_PoMo2N::buildRateMatrix(void)
{

 /*  
  INFORMATION ABOUT THE PoMo2N RATE MATRIX

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


  //mutations
  //mutations
  m[0][2]  = N*mu[0];  //{Na0} -> {(N-1)a0,1a1}
  m[1][N]  = N*mu[1];  //{Na1} -> {1a0,(N-1)a1}
  diagonal[0] = m[0][2];
  diagonal[1] = m[1][N];

  //fixations
  //fixations
  m[2][0]  = (N-1.0)*phi[0]/((N-1)*phi[0]+phi[1]);  //{(N-1)a0,1a1} -> {Na0} 
  m[N][1]  = (N-1.0)*phi[1]/((N-1)*phi[1]+phi[0]);  //{1a0,(N-1)a1} -> {Na1} 
  diagonal[2] += m[2][0];
  diagonal[N] += m[N][1];

  //the pomo rate matrix is entirely defined by fixations and mutations if N=2
  if (N>2) {

    //frequency shifts from singletons
    m[2][3]   = (N-1.0)*phi[1]/((N-1)*phi[0]+phi[1]);  //{(N-1)a0,1a1} -> {(N-2)a0,2a1}
    m[N][N-1] = (N-1.0)*phi[0]/((N-1)*phi[1]+phi[0]);  //{1a0,(N-1)a1} -> {2a0,(N-2)a1}
    diagonal[2] += m[2][3] ;
    diagonal[N] += m[N][N-1];

     //frequency shifts for all the other polymorphic states
    if (N>3) {

      //polymorphic states are populated
      for (int n=2; n<(N-1); n++){

        //ai aj
        m[n+1][2+n]  = n*(N-n)*phi[1]/(n*phi[1]+(N-n)*phi[0]); //{naj,(N-n)ai} -> {(n+1)aj,(N-n-1)ai}
        m[n+1][n]    = n*(N-n)*phi[0]/(n*phi[1]+(N-n)*phi[0]); //{naj,(N-n)ai} -> {(n-1)aj,(N-n+1)ai}
        diagonal[n+1] += m[n+1][2+n] + m[n+1][n];

      }

    }

  }


  // set the diagonal values and calculate the reciprocal of the average rate
  for (int j=0; j< num_states; j++){
    m[j][j] = -diagonal[j];
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
void RateMatrix_PoMo2N::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{

  // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
  // Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
  double t = rate * (startAge - endAge);
  exponentiateMatrixByScalingAndSquaring(t, P );

  return;
}


RateMatrix_PoMo2N* RateMatrix_PoMo2N::clone( void ) const
{
    return new RateMatrix_PoMo2N( *this );
}


std::vector<double> RateMatrix_PoMo2N::getStationaryFrequencies( void ) const
{
  return stationaryVector;
}




void RateMatrix_PoMo2N::setN( long & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMo2N::setMu( const std::vector<double> &m )
{
    mu = m;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMo2N::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}


void RateMatrix_PoMo2N::update( void )
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
