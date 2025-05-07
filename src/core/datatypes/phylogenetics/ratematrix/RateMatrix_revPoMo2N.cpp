#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMo2N.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_revPoMo2N::RateMatrix_revPoMo2N(long num_states, long in_n ) : TimeReversibleRateMatrix( num_states ),
    N( in_n ),
    pi( 2, 0.5 ),
    rho( 0.01 ),
    phi( 2, 1.0 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMo2N::RateMatrix_revPoMo2N(const RateMatrix_revPoMo2N& m) : TimeReversibleRateMatrix( m ),
    N( m.N ),
    pi( m.pi ),
    rho( m.rho ),
    phi( m.phi )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_revPoMo2N::~RateMatrix_revPoMo2N(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMo2N& RateMatrix_revPoMo2N::operator=(const RateMatrix_revPoMo2N &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        N                   = r.N;
        pi                  = r.pi;
        rho                 = r.rho;
        phi                 = r.phi;

        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_revPoMo2N::calculateCijk(void)
{
    
    if ( eigen_system->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = eigen_system->getEigenvectors();
        const MatrixReal& iev = eigen_system->getInverseEigenvectors();
        double* pc = &c_ijk[0];
        for (size_t i=0; i<num_states; i++)
        {
            for (size_t j=0; j<num_states; j++)
            {
                for (size_t k=0; k<num_states; k++)
                {
                    *(pc++) = ev[i][k] * iev[k][j];
                }
            }
        }
    }
    else
    {
        // complex case
        const MatrixComplex& cev  = eigen_system->getComplexEigenvectors();
        const MatrixComplex& ciev = eigen_system->getComplexInverseEigenvectors();
        std::complex<double>* pc = &cc_ijk[0];
        for (size_t i=0; i<num_states; i++)
        {
            for (size_t j=0; j<num_states; j++)
            {
                for (size_t k=0; k<num_states; k++)
                {
                    *(pc++) = cev[i][k] * ciev[k][j];
                }
            }
        }
    }
}


/** Calculate the transition probabilities */
void RateMatrix_revPoMo2N::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    double t = rate * (startAge - endAge);
    if ( eigen_system->isComplex() == false )
    {
        tiProbsEigens(t, P);
    }
    else
    {
        tiProbsComplexEigens(t, P);
    }
}


RateMatrix_revPoMo2N* RateMatrix_revPoMo2N::clone( void ) const
{
    return new RateMatrix_revPoMo2N( *this );
}


double RateMatrix_revPoMo2N::calculateReciprocalExpectedDivergence( void )
{

  // first the numerator

  double sum_n = 0.0;

  for (int n=1; n<(N+1); n++) {
    sum_n += pi[0]*pi[1]*rho*pow(phi[1],n-1)*pow(phi[0],N-n);
  }

  // then the denominator

  double sum_d = 0.0;
  double drift_coefficient;

  for (int n=1; n<N; n++) {
    drift_coefficient = 1.0/(n*(N-n));
    sum_d += pi[0]*pi[1]*rho*pow(phi[1],n-1)*pow(phi[0],N-n-1)*(n*phi[1]+(N-n)*phi[0])*drift_coefficient;
  }


  // finaly the rate (the reciprocal)
  double rRate;
  rRate = ( pi[0]*pow(phi[0],N-1) +
            pi[1]*pow(phi[1],N-1) + N*sum_d ) / ( 2.0*N*sum_n );

  return rRate;

}


/*populating the rate matrix*/
void RateMatrix_revPoMo2N::computeOffDiagonal( void )
{
    
  MatrixReal& m = *the_rate_matrix;

  /*  
  INFORMATION ABOUT THE PoMo2N RATE MATRIX

  It includes both fixed and polymorphic sites, but only biallelic states are considered.

  The pomo rate matrices defined here first list the fixed states {Na0}, {Na1} ...,
  these occupying positions 0:1 (only two fxed states with two alles), and then polymorphic states.

  2 alleles comprise 1 pairwise combinations of alleles.
  This is the number of edges in the PoMo2 state-space. This single edge comprises N-1 polymorphic states, 
  each of which represents a state of allelic frequencies in the population (summing to N).
  Example for a random allele pair aiaj: {(N-1)a0,1a1}, {(N-2)a0,2a1}, ..., {1a0,(N-1)a1}

  We say the a0a1 polymorphic states sit at edge 0. Thus, {(N-1)a0,1a1} sits at position 2 and 
  {1a0,(N-1)a1} sits N-2 positions further (i.e., N).
  */


  //populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++){
    for (int j=0; j< num_states; j++){
      m[i][j] = 0.0;        
    }
  }

  // get the expected divergence (or number of evens) per unit of time
  // normalize rate matrix such that one event happens per unit time.
  double rRate = calculateReciprocalExpectedDivergence();
  //std::cout << "rate:" << rRate << "\n";


  //mutations
  m[0][2]  = N*rho*pi[1]*rRate;  //{Na0} -> {(N-1)a0,1a1}
  m[1][N]  = N*rho*pi[0]*rRate;  //{Na1} -> {1a0,(N-1)a1}

  //fixations
  m[2][0]  = (N-1.0)*phi[0]*rRate/((N-1)*phi[0]+phi[1]);  //{(N-1)a0,1a1} -> {Na0} 
  m[N][1]  = (N-1.0)*phi[1]*rRate/((N-1)*phi[1]+phi[0]);  //{1a0,(N-1)a1} -> {Na1} 


  //the pomo rate matrix is entirely defined by fixations and mutations if N=2
  if (N>2) {

    //frequency shifts from singletons
    m[2][3]   = (N-1.0)*phi[1]*rRate/((N-1)*phi[0]+phi[1]);  //{(N-1)a0,1a1} -> {(N-2)a0,2a1}
    m[N][N-1] = (N-1.0)*phi[0]*rRate/((N-1)*phi[1]+phi[0]);  //{1a0,(N-1)a1} -> {2a0,(N-2)a1}

    //frequency shifts for all the other polymorphic states
    if (N>3) {

      for (int n=2; n<(N-1); n++){

        //ai aj
        m[n+1][2+n]  = n*(N-n)*phi[1]*rRate/(n*phi[1]+(N-n)*phi[0]); //{naj,(N-n)ai} -> {(n+1)aj,(N-n-1)ai}
        m[n+1][n]    = n*(N-n)*phi[0]*rRate/(n*phi[1]+(N-n)*phi[0]); //{naj,(N-n)ai} -> {(n-1)aj,(N-n+1)ai}

      }

    }
  }


  //int n_states = N+1;
  //std::vector<double> stationary_freqs; 
  //stationary_freqs = getStationaryFrequencies();

  //for (int i=0; i<n_states; i++) {
  //  std::cout << "sf" << i << ":" << stationary_freqs[i] << "\n";
  //}


  // set flags
  needs_update = true;

}

std::vector<double> RateMatrix_revPoMo2N::getStationaryFrequencies( void ) const
{

  // calculating the normalization constant

  double nc = pi[0]*pow(phi[0],N-1) +
              pi[1]*pow(phi[1],N-1) ;

  double drift_coefficient;

  for (int n=1; n<N; n++) {
    drift_coefficient = 1.0*N/(n*(N-n));
    nc += pi[0]*pi[1]*rho*pow(phi[1],n-1)*pow(phi[0],N-n-1)*(n*phi[1]+(N-n)*phi[0])*drift_coefficient;
  }

  
  // calculating the stationary vector

  double rnc = 1.0/nc;
  std::vector<double> stationary_freqs(N+1,0.0);

  stationary_freqs[0]  = pi[0]*pow(phi[0],N-1)*rnc;
  stationary_freqs[1]  = pi[1]*pow(phi[1],N-1)*rnc;

  for (int n=1; n<N; n++) {

    drift_coefficient = 1.0*N/(n*(N-n));
    stationary_freqs[1+n    ] = pi[0]*pi[1]*rho*pow(phi[1],n-1)*pow(phi[0],N-n-1)*(n*phi[1]+(N-n)*phi[0])*drift_coefficient*rnc;

  }


  return stationary_freqs;

}

/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMo2N::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = eigen_system->getRealEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<double> eigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        eigValExp[s] = exp(eigenValue[s] * t);
    }
    
    // calculate the transition probabilities
    const double* ptr = &c_ijk[0];
    double*         p = P.theMatrix;
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++, ++p)
        {
            double sum = 0.0;
            for (size_t s=0; s<num_states; s++)
            {
                sum += (*ptr++) * eigValExp[s];
            }
            
            //            P[i][j] = (sum < 0.0) ? 0.0 : sum;
            (*p) = (sum < 0.0) ? 0.0 : sum;
        }
    }
}


/** Calculate the transition probabilities for the complex case */
void RateMatrix_revPoMo2N::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = eigen_system->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = eigen_system->getImagEigenvalues();
    
    // precalculate the product of the eigenvalue and the branch length
    std::vector<std::complex<double> > ceigValExp(num_states);
    for (size_t s=0; s<num_states; s++)
    {
        std::complex<double> ev = std::complex<double>(eigenValueReal[s], eigenValueComp[s]);
        ceigValExp[s] = exp(ev * t);
    }
    
    // calculate the transition probabilities
    const std::complex<double>* ptr = &cc_ijk[0];
    for (size_t i=0; i<num_states; i++)
    {
        for (size_t j=0; j<num_states; j++)
        {
            std::complex<double> sum = std::complex<double>(0.0, 0.0);
            for (size_t s=0; s<num_states; s++)
                sum += (*ptr++) * ceigValExp[s];
            P[i][j] = (sum.real() < 0.0) ? 0.0 : sum.real();
        }
    }
}


void RateMatrix_revPoMo2N::setN( long & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_revPoMo2N::setPi( const Simplex &bf )
{
    pi = bf;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMo2N::setRho( double &ex )
{
    rho = ex;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMo2N::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}



/** Update the eigen system */
void RateMatrix_revPoMo2N::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMo2N::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale
        //rescaleToAverageRate( 3.0 );
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}



