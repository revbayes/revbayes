#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMoTwo4N.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>
#include <boost/math/special_functions/digamma.hpp>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_revPoMoTwo4N::RateMatrix_revPoMoTwo4N( void ) : TimeReversibleRateMatrix( 10 ),
    N( 2 ),
    pi(4, 0.25),
    rho( 6, 0.01 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMoTwo4N::RateMatrix_revPoMoTwo4N(const RateMatrix_revPoMoTwo4N& m) : TimeReversibleRateMatrix( m ),
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
RateMatrix_revPoMoTwo4N::~RateMatrix_revPoMoTwo4N(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMoTwo4N& RateMatrix_revPoMoTwo4N::operator=(const RateMatrix_revPoMoTwo4N &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        N    				= r.N;
        pi                  = r.pi;
        rho                 = r.rho;

        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_revPoMoTwo4N::calculateCijk(void)
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
void RateMatrix_revPoMoTwo4N::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_revPoMoTwo4N* RateMatrix_revPoMoTwo4N::clone( void ) const
{
    return new RateMatrix_revPoMoTwo4N( *this );
}




void RateMatrix_revPoMoTwo4N::computeOffDiagonal( void )
{
    
  MatrixReal& m = *the_rate_matrix;

  /*  
  INFORMATION ABOUT PoMoTwo AND PoMoThree

  The idea of the virtual PoMos (Two and Three) resumes to mimicking a population dynamic that unfolds 
  on the effective population N, using a virtual population of smaller size M.
  M is equal to 2 in PoMoTwo and 3 in PoMoThree. PoMo thee additionally includes selection.

  M defines a lighter and more efficient state-space. 
  By matching the expected diversity (i.e., the proportion of fixed and polymorphic sites) in both 
  the effective and the virtual populations, one can obtain scaling laws for the mutation rates and 
  selection coefficients (see Borges et al. 2019 Genetics). 

  These are intuitive for the M=2 case, in which the mutation rates are scaled by the 
  harmonic number of N-1:
  mu'_ij = mu_ij*H_{N-1}  or equivalently rho'_ij = rho_ij*H_{N-1}  for the reverisible PoMo 
  mu' and rho' correspond to the mutation rate and exchageability in a virtual population of 2 individuals

  Like the standard PoMos the virtual PoMos include both fixed and polymorphic sites: 
            virtual_individuals_M   n_states   fixed_states   polymorphic_states
  PoMoTwo   2                       10         {2ai}          {1ai,1aj}
  PoMoThree 3                       16         {3ai}          {1ai,2aj} and {2ai,1aj}

  The pomo rate matrices defined here first list the fixed states {Na0}, {Na1} ...,
  these occupying positions 0:(K-1), and then polymorphic states, where K is the number of alleles

  K alleles comprise (K*K-K)/2 pairwise combinations of alleles.
  This is the number of edges in the pomo state-space. Each edge comprises M-1 polymorphic states.
  The polymorphic edges are listed in the following order a0a1, a0a2, a0a3, ..., a(K-2)aK-1
  */

    
    
  // populating the rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i<10 ; i++){
    for (int j=0; j<10; j++ ){
        m[i][j] = 0.0;
    }
  }


  // calculating the harmonic number of N-1
  // used to scale the mutation rates (or exchangeabilities)
  double harmonic_number = boost::math::digamma(N) - boost::math::digamma(1.0);
  double r = N*harmonic_number/2.0;	


  // get the expected divergence (or number of evens) per unit of time
  // normalize rate matrix such that one event happens per unit time.
  // a common quantity to the numerator and denominator
  double value = pi[0]*pi[1]*rho[0] +
                 pi[0]*pi[2]*rho[1] +
                 pi[0]*pi[3]*rho[2] +
                 pi[1]*pi[2]*rho[3] +
                 pi[1]*pi[3]*rho[4] +
                 pi[2]*pi[3]*rho[5] ;

  // receiprocal of the rate
  double rRate = ( 1.0 + 4.0*r*value ) / ( 8.0*r*value );

  
  // Mutations
  m[0][4] = 2*rho[0]*r*pi[1]*rRate;    //mutation AC
  m[0][5] = 2*rho[1]*r*pi[2]*rRate;    //mutation AG
  m[0][6] = 2*rho[2]*r*pi[3]*rRate;    //mutation AT

  m[1][4] = 2*rho[0]*r*pi[0]*rRate;    //mutation CA
  m[1][7] = 2*rho[3]*r*pi[2]*rRate;    //mutation CG
  m[1][8] = 2*rho[4]*r*pi[3]*rRate;    //mutation CT

  m[2][5] = 2*rho[1]*r*pi[0]*rRate;    //mutation GA
  m[2][7] = 2*rho[3]*r*pi[1]*rRate;    //mutation GC
  m[2][9] = 2*rho[5]*r*pi[3]*rRate;    //mutation GT

  m[3][6] = 2*rho[2]*r*pi[0]*rRate;    //mutation TA
  m[3][8] = 2*rho[4]*r*pi[1]*rRate;    //mutation TC
  m[3][9] = 2*rho[5]*r*pi[2]*rRate;    //mutation TG


  // Fixations
  // PoMoTwo only accouts for genetic drift
  // selection is not indentifyable with two virtual individuals
  m[4][0]   = 0.5*rRate;             //A fixed
  m[4][1]   = 0.5*rRate;             //C fixed

  m[5][0]   = 0.5*rRate;             //A fixed
  m[5][2]   = 0.5*rRate;             //G fixed

  m[6][0]   = 0.5*rRate;             //A fixed
  m[6][3]   = 0.5*rRate;             //T fixed

  m[7][1]   = 0.5*rRate;             //C fixed
  m[7][2]   = 0.5*rRate;             //G fixed

  m[8][1]   = 0.5*rRate;             //C fixed
  m[8][3]   = 0.5*rRate;             //T fixed
  
  m[9][2]   = 0.5*rRate;             //G fixed
  m[9][3]   = 0.5*rRate;             //T fixed


  // set flags
  needs_update = true;

}

std::vector<double> RateMatrix_revPoMoTwo4N::getStationaryFrequencies( void ) const
{

  // calculating the harmonic number of N-1
  // used to scale the mutation rates (or exchangeabilities)
  double harmonic_number = boost::math::digamma(N) - boost::math::digamma(1.0);
  double r = N*harmonic_number/2.0;

  // calculating the normalization constant

  double nc = 1.0 +
              4.0*pi[0]*pi[1]*rho[0]*r + 
              4.0*pi[0]*pi[2]*rho[1]*r + 
              4.0*pi[0]*pi[3]*rho[2]*r + 
              4.0*pi[1]*pi[2]*rho[3]*r + 
              4.0*pi[1]*pi[3]*rho[4]*r + 
              4.0*pi[2]*pi[3]*rho[5]*r ;


  // calculating the stationary vector

  double rnc = 1.0/nc;

  std::vector<double> stationary_freqs(16,0.0);

  stationary_freqs[0]  = pi[0]*rnc;
  stationary_freqs[1]  = pi[1]*rnc;
  stationary_freqs[2]  = pi[2]*rnc;
  stationary_freqs[3]  = pi[3]*rnc;
  stationary_freqs[4]  = pi[0]*pi[1]*rho[0]*r*4.0*rnc;
  stationary_freqs[5]  = pi[0]*pi[2]*rho[1]*r*4.0*rnc;
  stationary_freqs[6]  = pi[0]*pi[3]*rho[2]*r*4.0*rnc;
  stationary_freqs[7]  = pi[1]*pi[2]*rho[3]*r*4.0*rnc;
  stationary_freqs[8]  = pi[1]*pi[3]*rho[4]*r*4.0*rnc;
  stationary_freqs[9]  = pi[2]*pi[3]*rho[5]*r*4.0*rnc;

  return stationary_freqs;

}


/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMoTwo4N::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_revPoMoTwo4N::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_revPoMoTwo4N::setN(double ps)
{
    N = ps;
    // set flags
    needs_update = true;
}


void RateMatrix_revPoMoTwo4N::setPi( const std::vector<double> &f )
{
    pi = f;
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoTwo4N::setRho( const std::vector<double> &r )
{
    rho = r;
    // set flags
    needs_update = true;
}

/** Update the eigen system */
void RateMatrix_revPoMoTwo4N::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMoTwo4N::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale
        //rescaleToAverageRate(1.0);
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}



