#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMoNeutralM4N.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>
#include <boost/math/special_functions/digamma.hpp>

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_revPoMoNeutralM4N::RateMatrix_revPoMoNeutralM4N( std::int64_t num_states, std::int64_t in_m ) : TimeReversibleRateMatrix( num_states ),
    N( 2 ),
    M( in_m ),
    pi(4, 0.25),
    rho( 6, 0.01 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMoNeutralM4N::RateMatrix_revPoMoNeutralM4N(const RateMatrix_revPoMoNeutralM4N& m) : TimeReversibleRateMatrix( m ),
    N( m.N ),
    M( m.M ),
    pi( m.pi ),
    rho( m.rho )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_revPoMoNeutralM4N::~RateMatrix_revPoMoNeutralM4N(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMoNeutralM4N& RateMatrix_revPoMoNeutralM4N::operator=(const RateMatrix_revPoMoNeutralM4N &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        N    				        = r.N;
        M                   = r.M;
        pi                  = r.pi;
        rho                 = r.rho;

        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_revPoMoNeutralM4N::calculateCijk(void)
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
void RateMatrix_revPoMoNeutralM4N::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_revPoMoNeutralM4N* RateMatrix_revPoMoNeutralM4N::clone( void ) const
{
    return new RateMatrix_revPoMoNeutralM4N( *this );
}




void RateMatrix_revPoMoNeutralM4N::computeOffDiagonal( void )
{
    
  MatrixReal& rm = *the_rate_matrix;

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
  for (int i=0; i<num_states ; i++){
    for (int j=0; j<num_states; j++ ){
        rm[i][j] = 0.0;
    }
  }


  // calculating the ration of harmonic numbers
  // this scaling helps scaling the virtual pomo matrix so times is in terms of 
  // moran generations in the effective population
  // plus guarantees that both have the same stationary distribution 
  double ratio_harmonic_numbers = (boost::math::digamma(N) - boost::math::digamma(1.0))/(boost::math::digamma(M) - boost::math::digamma(1.0));
  double ratio_harmonic_numbers2 = ratio_harmonic_numbers*ratio_harmonic_numbers;	

  
  // Mutations
  // from A
  rm[0][4]     = M*rho[0]*pi[1]*ratio_harmonic_numbers2;
  rm[0][M+3]   = M*rho[1]*pi[2]*ratio_harmonic_numbers2;
  rm[0][2*M+2] = M*rho[2]*pi[3]*ratio_harmonic_numbers2;
  
  // from C
  rm[1][M+2]   = M*rho[0]*pi[0]*ratio_harmonic_numbers2;
  rm[1][3*M+1] = M*rho[3]*pi[2]*ratio_harmonic_numbers2;
  rm[1][4*M]   = M*rho[4]*pi[3]*ratio_harmonic_numbers2;
  
  // from G
  rm[2][2*M+1] = M*rho[1]*pi[0]*ratio_harmonic_numbers2;
  rm[2][4*M-1] = M*rho[3]*pi[1]*ratio_harmonic_numbers2;
  rm[2][5*M-1] = M*rho[5]*pi[3]*ratio_harmonic_numbers2;
  
  // from T
  rm[3][3*M]   = M*rho[2]*pi[0]*ratio_harmonic_numbers2;
  rm[3][5*M-2] = M*rho[4]*pi[1]*ratio_harmonic_numbers2;
  rm[3][6*M-3] = M*rho[5]*pi[2]*ratio_harmonic_numbers2;


  // Fixations
  // frequency shifts from singletons
  double p = (M-1)*ratio_harmonic_numbers/N;

  // AC
  rm[4][0]     = p;
  rm[M+2][1]   = p;
  
  // AC
  rm[M+3][0]   = p;
  rm[2*M+1][2] = p;
  
  // AT
  rm[2*M+2][0] = p;
  rm[3*M][3]   = p;
  
  // CG
  rm[3*M+1][1]   = p;
  rm[4*M-1][2]   = p;
  
  // CT
  rm[4*M][1]    = p;
  rm[5*M-2][3]  = p;
  
  // GT
  rm[5*M-1][2]    = p;
  rm[6*M-3][3]    = p;


if(M>2){
    
    //AC
    rm[4][5]       = p;
    rm[M+2][M+1]   = p;
    
    //AG
    rm[M+3][M+4]     = p;
    rm[2*M+1][2*M] = p;
    
    //AT
    rm[2*M+2][2*M+3] = p;
    rm[3*M][3*M-1]   = p;
    
    //CG
    rm[3*M+1][3*M+2] = p;
    rm[4*M-1][4*M-2]   = p;
    
    //CT
    rm[4*M][4*M+1] = p;
    rm[5*M-2][5*M-3] = p;
    
    //GT
    rm[5*M-1][5*M]   = p;
    rm[6*M-3][6*M-4] = p;
    
    if (M>3){

      // remaining p;olymorp;hic states
      for (int m=2; m<(M-1); ++m ){
        
        p = m*(M-m)*ratio_harmonic_numbers/N;
        
        // AC
        rm[3+m][4+m] = p;
        rm[3+m][2+m] = p;
        
        // AG
        rm[M+2+m][M+3+m] = p;
        rm[M+2+m][M+1+m] = p;
        
        // AT
        rm[2*M+1+m][2*M+2+m] = p;
        rm[2*M+1+m][2*M+m] = p;
        
        // CG
        rm[3*M+m][3*M+1+m] = p;
        rm[3*M+m][3*M+m -1] = p;
        
        // CT
        rm[4*M+m-1][4*M+m] = p;
        rm[4*M+m-1][4*M+m-2] = p;
        
        // GT
        rm[5*M+m-2][5*M+m-1] = p;
        rm[5*M+m-2][5*M+m-3] = p;
        
      }
    }
  }


  // set flags
  needs_update = true;

}

std::vector<double> RateMatrix_revPoMoNeutralM4N::getStationaryFrequencies( void ) const
{

  // calculating the harmonic number of N-1
  // used to scale the mutation rates (or exchangeabilities)
  double scaling = N*(boost::math::digamma(N) - boost::math::digamma(1.0))/(M*(boost::math::digamma(M) - boost::math::digamma(1.0)));

  // calculating the normalization constant

  double nc = 1.0 +
              ( pi[0]*pi[1]*rho[0] + 
                pi[0]*pi[2]*rho[1] + 
                pi[0]*pi[3]*rho[2] + 
                pi[1]*pi[2]*rho[3] + 
                pi[1]*pi[3]*rho[4] + 
                pi[2]*pi[3]*rho[5] )*2*N*(boost::math::digamma(N) - boost::math::digamma(1.0));


  // calculating the stationary vector

  double rnc = 1.0/nc;

  std::vector<double> stationary_freqs(num_states,0.0);

  stationary_freqs[0]  = pi[0]*rnc;
  stationary_freqs[1]  = pi[1]*rnc;
  stationary_freqs[2]  = pi[2]*rnc;
  stationary_freqs[3]  = pi[3]*rnc;

  double drift_coefficient;

  for (int m=1; m<M; m++) {

    drift_coefficient = 1.0*M*M/(m*(M-m));

    stationary_freqs[3+m    ]  = pi[0]*pi[1]*rho[0]*scaling*drift_coefficient*rnc;
    stationary_freqs[M+m+2  ]  = pi[0]*pi[2]*rho[1]*scaling*drift_coefficient*rnc;
    stationary_freqs[2*M+m+1]  = pi[0]*pi[3]*rho[2]*scaling*drift_coefficient*rnc;
    stationary_freqs[3*M+m  ]  = pi[1]*pi[2]*rho[3]*scaling*drift_coefficient*rnc;
    stationary_freqs[4*M+m-1]  = pi[1]*pi[3]*rho[4]*scaling*drift_coefficient*rnc;
    stationary_freqs[5*M+m-2]  = pi[2]*pi[3]*rho[5]*scaling*drift_coefficient*rnc;

  }

  //std::vector<double> stationary_freqs = getStationaryFrequencies();

  //for (int i=0; i<num_states; i++) {
  //  std::cout << " " << i << "  =  " << stationary_freqs[i] << "\n";
  //}

  return stationary_freqs;

}


/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMoNeutralM4N::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_revPoMoNeutralM4N::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_revPoMoNeutralM4N::setN(std::int64_t ps)
{
    N = ps;
    // set flags
    needs_update = true;
}


void RateMatrix_revPoMoNeutralM4N::setM(std::int64_t vps)
{
    M = vps;
    // set flags
    needs_update = true;
}



void RateMatrix_revPoMoNeutralM4N::setPi( const std::vector<double> &f )
{
    pi = f;
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoNeutralM4N::setRho( const std::vector<double> &r )
{
    rho = r;
    // set flags
    needs_update = true;
}

/** Update the eigen system */
void RateMatrix_revPoMoNeutralM4N::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMoNeutralM4N::update( void )
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



