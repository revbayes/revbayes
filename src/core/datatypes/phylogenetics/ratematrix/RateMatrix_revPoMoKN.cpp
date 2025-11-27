#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMoKN.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_revPoMoKN::RateMatrix_revPoMoKN(std::int64_t num_states, std::int64_t in_k, std::int64_t in_n, std::int64_t in_nex) : TimeReversibleRateMatrix( num_states ),
    K( in_k ),
    N( in_n ),
    pi( in_k, 1.0/in_k ),
    rho( in_nex, 0.01 ),
    phi( in_k, 1.0 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMoKN::RateMatrix_revPoMoKN(const RateMatrix_revPoMoKN& m) : TimeReversibleRateMatrix( m ),
    K( m.K ),
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
RateMatrix_revPoMoKN::~RateMatrix_revPoMoKN(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMoKN& RateMatrix_revPoMoKN::operator=(const RateMatrix_revPoMoKN &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        K                   = r.K;
        N                   = r.N;
        pi                  = r.pi;
        rho                 = r.rho;
        phi                 = r.phi;

        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_revPoMoKN::calculateCijk(void)
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
void RateMatrix_revPoMoKN::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_revPoMoKN* RateMatrix_revPoMoKN::clone( void ) const
{
    return new RateMatrix_revPoMoKN( *this );
}

double RateMatrix_revPoMoKN::calculateReciprocalExpectedDivergence( void )
{

  // calculating the denominator of the expected rate
  double sum_d = 0.0;

  // summing the stationary weights for the fixed states 
  for (int i=0; i<K; i++) {
    sum_d += pi[i]*pow(phi[i],N-1);
  }

  // now for the polymorphic states 
  //first edge
  int E = 0;

  double drift_coefficient;

  for (int i=0; i<K; i++){
    for (int j=i+1; j<K; j++){
      for (int n=1; n<N; n++) {

        // sum weight
        drift_coefficient = 1.0*N/(n*(N-n));
        sum_d += pi[i]*pi[j]*rho[E]*pow(phi[j],n-1)*pow(phi[i],N-n-1)*(n*phi[j]+(N-n)*phi[i])*drift_coefficient;

      }
      E += 1;
    }
  }



  // now the numerator of the expected rate
  double sum_n = 0.0;
  E = 0;

  for (int i=0; i<K; i++){
    for (int j=i+1; j<K; j++){
      for (int n=1; n<(N+1); n++) {

        sum_n += pi[i]*pi[j]*rho[E]*pow(phi[j],n-1)*pow(phi[i],N-n);

      }
      E += 1;

    }
  }

  // finaly the rate (the reciprocal)
  double rRate;
  rRate = sum_d / ( 2.0*N*sum_n ) ;
  return rRate;

}


/*populating the rate matrix*/
void RateMatrix_revPoMoKN::computeOffDiagonal( void )
{
    
  MatrixReal& m = *the_rate_matrix;

 /*  
  INFORMATION ABOUT THE revPoMoKN RATE MATRIX

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



  //populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++){
    for (int j=0; j< num_states; j++){
      m[i][j] = 0.0;        
    }
  }

  double rRate;
  rRate = calculateReciprocalExpectedDivergence();
  //std::cout << "rate:" << rRate << "\n";

  //first edge
  int E = 0;


  //these for loops go through the (K*K-K) edges of the pomo state-space
  //their represent all the possible pairwise combinations of alleles: a0a1, a1a2, ..., aK-2aK-1
  for (int i=0; i<K; i++){
    for (int j=i+1; j<K; j++){

      //mutations
      m[i][K+E*N-E]          = N*pi[j]*rho[E]*rRate;    //{Nai} -> {(N-1)ai,1aj}
      m[j][K+(E+1)*(N-1)-1]  = N*pi[i]*rho[E]*rRate;    //{Naj} -> {1ai,(N-1)aj}

      //fixations
      m[K+E*N-E]        [i]  = (N-1.0)*phi[i]*rRate/(phi[j] + (N-1.0)*phi[i]);  //{(N-1)ai,1aj} -> {Nai} 
      m[K+(E+1)*(N-1)-1][j]  = (N-1.0)*phi[j]*rRate/(phi[i] + (N-1.0)*phi[j]);  //{1ai,(N-1)aj} -> {Naj} 


      //the pomo rate matrix is entirely defined by fixations and mutations if N=2
      if (N>2) {

        //frequency shifts from singletons
        m[K+E*N-E]        [K+E*N-E+1]       = (N-1.0)*phi[j]*rRate/(phi[j] + (N-1.0)*phi[i]);  //{(N-1)ai,1aj} -> {(N-2)ai,2aj}
        m[K+(E+1)*(N-1)-1][K+(E+1)*(N-1)-2] = (N-1.0)*phi[i]*rRate/(phi[i] + (N-1.0)*phi[j]);  //{1ai,(N-1)aj} -> {2ai,(N-2)aj}

        //frequency shifts for all the other polymorphic states
        if (N>3) {

          //polymorphic states are populated
          for (int n=2; n<(N-1); n++){

            //populates the first half of the polymorphic edge aiaj
            m[K+E*N-E+n-1]    [K+E*N-E+n]         = n*(N-n)*phi[j]*rRate/(n*phi[j] + (N-n)*phi[i]); //{(N-n)ai,naj} -> {(N-n-1)ai,(n+1)aj}
            m[K+E*N-E+n-1]    [K+E*N-E+n-2]       = n*(N-n)*phi[i]*rRate/(n*phi[j] + (N-n)*phi[i]); //{(N-n)ai,naj} -> {(N-n+1)ai,(n-1)aj}



          }

        }

      }

      //update edge
      E += 1; 

    }
  }


  //int n_states = K+(K*K-K)*(N-1)/2;
  //std::vector<double> stationary_freqs; 
  //stationary_freqs = getStationaryFrequencies();

  //for (int i=0; i<n_states; i++) {
  //  std::cout << "sf" << i << ":" << stationary_freqs[i] << "\n";
  //}



  // set flags
  needs_update = true;

}



std::vector<double> RateMatrix_revPoMoKN::getStationaryFrequencies( void ) const
{

  // calculating the normalization constant and stationary weights
  int n_states = K+(K*K-K)*(N-1)/2;
  std::vector<double> stationary_freqs(n_states,0.0);

  // normalization constant
  double nc = 0.0;

  // the fixed states first
  for (int i=0; i<K; i++) {
    stationary_freqs[i] = pi[i]*pow(phi[i],N-1);
    nc += stationary_freqs[i];
  }

  // now for the polymorphic states 
  //first edge
  int index = K;
  int E = 0;
  double drift_coefficient;

  for (int i=0; i<K; i++){
    for (int j=i+1; j<K; j++){
      for (int n=1; n<N; n++) {

        // sum weight
        drift_coefficient = 1.0*N/(n*(N-n));
        stationary_freqs[index] = pi[i]*pi[j]*rho[E]*pow(phi[j],n-1)*pow(phi[i],N-n-1)*(n*phi[j]+(N-n)*phi[i])*drift_coefficient;
        nc += stationary_freqs[index];
        index += 1;

      }

      //update edge
      E += 1; 

    }
  }


  // normalizing the stationary vector
  double rnc = 1.0/nc;

  for (int i=0; i<n_states; i++) {
    stationary_freqs[i] = stationary_freqs[i]*rnc;
  }

  return stationary_freqs;

}


/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMoKN::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_revPoMoKN::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_revPoMoKN::setK( std::int64_t & na )
{   
    K = na;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_revPoMoKN::setN( std::int64_t & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_revPoMoKN::setPi( const Simplex &bf )
{
    pi = bf;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoKN::setRho( const std::vector<double> &ex )
{
    rho = ex;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoKN::setPhi( const std::vector<double> &f )
{
    phi = f;
    
    // set flags
    needs_update = true;
}



/** Update the eigen system */
void RateMatrix_revPoMoKN::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMoKN::update( void )
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



