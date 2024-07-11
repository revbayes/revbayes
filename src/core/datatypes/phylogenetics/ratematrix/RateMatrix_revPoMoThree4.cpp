#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMoThree4.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>
#include <boost/math/special_functions/digamma.hpp>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_revPoMoThree4::RateMatrix_revPoMoThree4( void ) : TimeReversibleRateMatrix( 16 ),
    pi(4, 0.25),
    rho( 6, 0.01 ),
    phi( 4, 1.0 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMoThree4::RateMatrix_revPoMoThree4(const RateMatrix_revPoMoThree4& m) : TimeReversibleRateMatrix( m ),
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
RateMatrix_revPoMoThree4::~RateMatrix_revPoMoThree4(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMoThree4& RateMatrix_revPoMoThree4::operator=(const RateMatrix_revPoMoThree4 &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        pi                  = r.pi;
        rho                 = r.rho;
        phi                 = r.phi;

        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_revPoMoThree4::calculateCijk(void)
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
void RateMatrix_revPoMoThree4::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_revPoMoThree4* RateMatrix_revPoMoThree4::clone( void ) const
{
    return new RateMatrix_revPoMoThree4( *this );
}


void RateMatrix_revPoMoThree4::computeOffDiagonal( void )
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
  for (int i=0; i<16 ; i++){
    for (int j=0; j<16; j++ ){
        m[i][j] = 0.0;
    }
  }

  
	
  //mutations
  //AC
  m[0][4]   = 3*rho[0]*pi[1];    //{3A} -> {2A,1C}
  m[1][5]   = 3*rho[0]*pi[0];    //{3C} -> {1A,2C}

  //AG
  m[0][6]   = 3*rho[1]*pi[2];    //{3A} -> {2A,1G}
  m[2][7]   = 3*rho[1]*pi[0];    //{3G} -> {1A,2G}

  //AT
  m[0][8]   = 3*rho[2]*pi[3];    //{3A} -> {2A,1T}
  m[3][9]   = 3*rho[2]*pi[0];    //{3T} -> {1A,2T}

  //CG
  m[1][10]  = 3*rho[3]*pi[2];    //{3C} -> {2C,1G}
  m[2][11]  = 3*rho[3]*pi[1];    //{3G} -> {1C,2G}

  //CT
  m[1][12]  = 3*rho[4]*pi[3];    //{3C} -> {2C,1T}
  m[3][13]  = 3*rho[4]*pi[1];    //{3T} -> {1C,2T}

  //GT
  m[2][14]  = 3*rho[5]*pi[3];    //{3G} -> {2G,1T}
  m[3][15]  = 3*rho[5]*pi[2];    //{3T} -> {1G,2T}


  //fixations
  //AC
  m[4][0]   = phi[0]/(2.0*phi[0]+1.0*phi[1]);  //{2A,1C} -> {3A} 
  m[5][1]   = phi[1]/(2.0*phi[1]+1.0*phi[0]);  //{1A,2C} -> {3C} 

  //AG
  m[6][0]   = phi[0]/(2.0*phi[0]+1.0*phi[2]);  //{2A,1G} -> {3A} 
  m[7][2]   = phi[2]/(2.0*phi[2]+1.0*phi[0]);  //{1A,2G} -> {3G} 

  //AT
  m[8][0]   = phi[0]/(2.0*phi[0]+1.0*phi[3]);  //{2A,1T} -> {3A} 
  m[9][3]   = phi[3]/(2.0*phi[3]+1.0*phi[0]);  //{1A,2T} -> {3T} 

  //CG
  m[10][1]  = phi[1]/(2.0*phi[1]+1.0*phi[2]);  //{2C,1G} -> {3C} 
  m[11][2]  = phi[2]/(2.0*phi[2]+1.0*phi[1]);  //{1C,2G} -> {3G}

  //CT
  m[12][1]  = phi[1]/(2.0*phi[1]+1.0*phi[3]);  //{2C,1T} -> {3C} 
  m[13][3]  = phi[3]/(2.0*phi[3]+1.0*phi[1]);  //{1C,2T} -> {3T} 

  //GT
  m[14][2]  = phi[2]/(2.0*phi[2]+1.0*phi[3]);  //{2G,1T} -> {3G} 
  m[15][3]  = phi[3]/(2.0*phi[3]+1.0*phi[2]);  //{1G,2T} -> {3T} 

  // frequency shifts
  m[4][5]    = phi[1]/(2.0*phi[0]+1.0*phi[1]);  //{2A,1C} -> {1A,2C}
  m[5][4]    = phi[0]/(2.0*phi[1]+1.0*phi[0]);  //{1A,2C} -> {2A,1C}

  //AG
  m[6][7]    = phi[2]/(2.0*phi[0]+1.0*phi[2]);  //{2A,1G} -> {1A,2G}
  m[7][6]    = phi[0]/(2.0*phi[2]+1.0*phi[0]);  //{1A,2G} -> {2A,1G}

  //AT
  m[8][9]    = phi[3]/(2.0*phi[0]+1.0*phi[3]);  //{2A,1T} -> {1A,2T}
  m[9][8]    = phi[0]/(2.0*phi[3]+1.0*phi[0]);  //{1A,2T} -> {2A,1T}

  //CG
  m[10][11]  = phi[2]/(2.0*phi[1]+1.0*phi[2]);  //{2C,1G} -> {1C,2G}
  m[11][10]  = phi[1]/(2.0*phi[2]+1.0*phi[1]);  //{1C,2G} -> {2C,1G}

  //CT
  m[12][13]  = phi[3]/(2.0*phi[1]+1.0*phi[3]);  //{2C,1T} -> {1C,2T}
  m[13][12]  = phi[1]/(2.0*phi[3]+1.0*phi[1]);  //{1C,2T} -> {2C,1T}

  //GT
  m[14][15]  = phi[3]/(2.0*phi[2]+1.0*phi[3]);  //{2G,1T} -> {1G,2T}
  m[15] [14] = phi[2]/(2.0*phi[3]+1.0*phi[2]);  //{1G,2T} -> {2G,1T}

  // set flags
  needs_update = true;

}


/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMoThree4::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_revPoMoThree4::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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





void RateMatrix_revPoMoThree4::setPi( const std::vector<double> &f )
{
    pi = f;
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoThree4::setRho( const std::vector<double> &r )
{
    rho = r;
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoThree4::setPhi( const std::vector<double> &fc )
{
    phi = fc;
    // set flags
    needs_update = true;
}

/** Update the eigen system */
void RateMatrix_revPoMoThree4::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMoThree4::update( void )
{
    
    if ( needs_update )
    {
        // compute the off-diagonal values
        computeOffDiagonal();
        
        // set the diagonal values
        setDiagonal();
        
        // rescale
        //rescaleToAverageRate(e_rate);
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}



