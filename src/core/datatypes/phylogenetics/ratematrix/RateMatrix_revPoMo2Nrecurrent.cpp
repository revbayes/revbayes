#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMo2Nrecurrent.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_revPoMo2Nrecurrent::RateMatrix_revPoMo2Nrecurrent(long num_states, long in_n ) : TimeReversibleRateMatrix( num_states ),
    N( in_n ),
    pi( 2, 0.5 ),
    rho( 0.01 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMo2Nrecurrent::RateMatrix_revPoMo2Nrecurrent(const RateMatrix_revPoMo2Nrecurrent& m) : TimeReversibleRateMatrix( m ),
    N( m.N ),
    pi( m.pi ),
    rho( m.rho )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_revPoMo2Nrecurrent::~RateMatrix_revPoMo2Nrecurrent(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMo2Nrecurrent& RateMatrix_revPoMo2Nrecurrent::operator=(const RateMatrix_revPoMo2Nrecurrent &r)
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

        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_revPoMo2Nrecurrent::calculateCijk(void)
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
void RateMatrix_revPoMo2Nrecurrent::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_revPoMo2Nrecurrent* RateMatrix_revPoMo2Nrecurrent::clone( void ) const
{
    return new RateMatrix_revPoMo2Nrecurrent( *this );
}




/*populating the rate matrix*/
void RateMatrix_revPoMo2Nrecurrent::computeOffDiagonal( void )
{
    
  MatrixReal& m = *the_rate_matrix;

 /*  
  INFORMATION ABOUT THE revPoMo2Nrecurrent RATE MATRIX

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



  // populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++)
  {
      for (int j=0; j< num_states; j++)
      {
          m[i][j] = 0.0;
      }
  }
  

  //This needs to be CAREFULLY CHECKED!
  // first edge
  //int E = 0;
  //i=0, j=1


  // these for loops go through the (K*K-K) edges of the pomo state-space
  // their represent all the possible pairwise combinations of alleles: a0a1, a1a2, ..., aK-2aK-1
  // mutations

  double mu_01 = rho*pi[1];
  double mu_10 = rho*pi[0];


  m[0][2]    = N*mu_01;    //{Nai} -> {(N-1)ai,1aj}
  m[1][N]    = N*mu_10;    //{Naj} -> {1ai,(N-1)aj}
  
  // diagonal elements
  m[0][0] = -N*mu_01;
  m[1][1] = -N*mu_10;

  //fixations
  double fixation_rate = (N-1.0)/N; 

  m[2][0]  = fixation_rate + mu_10;  //{(N-1)ai,1aj} -> {Nai}
  m[N][1]  = fixation_rate + mu_01;  //{1ai,(N-1)aj} -> {Naj}
  
  // diagonal elements 
  m[2][2] = -2.0*fixation_rate -         mu_10 - (N-1.0)*mu_01;
  m[N][N] = -2.0*fixation_rate - (N-1.0)*mu_10 -         mu_01;

  // the pomo rate matrix is entirely defined by fixations and mutations if N=2
  
  double drift_rate;
  
  if (N>2)
  {
      
      //polymorphic states are populated
      for (int n=1; n<(N-1); n++)
      {

          drift_rate = 1.0*n*(N-n)/N;
                      
          //populates the first half of the polymorphic edge aiaj
          m[n+1][n+2]   = drift_rate + (N-n)*mu_01; //{(N-n)ai,naj} -> {(N-n-1)ai,(n+1)aj}
          m[N-n+1][N-n] = drift_rate + (N-n)*mu_10; //{(N-n)ai,naj} -> {(N-n+1)ai,(n-1)aj}

          // diagonal elements
          m[n+1]  [n+1]   = -2.0*drift_rate - (N-n)*mu_01 -     n*mu_10;
          //m[N-n+1][N-n+1] = -2.0*drift_rate -     n*mu_10 - (N-n)*mu_01;
                                      
      }
                                
  }
          
  // set flags
  needs_update = true;

}






/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMo2Nrecurrent::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_revPoMo2Nrecurrent::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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




void RateMatrix_revPoMo2Nrecurrent::setN( long & ni )
{
    N = ni;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_revPoMo2Nrecurrent::setPi( const Simplex &bf )
{
    pi = bf;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMo2Nrecurrent::setRho( double &ex )
{
    rho = ex;
    
    // set flags
    needs_update = true;
}




/** Update the eigen system */
void RateMatrix_revPoMo2Nrecurrent::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMo2Nrecurrent::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        //setDiagonal();
        
        // rescale
        //rescaleToAverageRate( 3.0 );
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}



