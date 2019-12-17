#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_PoMoBalance.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_PoMoBalance::RateMatrix_PoMoBalance( size_t ss, double in_n ) : TimeReversibleRateMatrix( ss ),
    N( in_n ),
    pi( 4, 0.25),
    rho( 6, 0.001 ),
    sigma(  4, 0.0 ),
    beta(  1.0 )
    
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_PoMoBalance::RateMatrix_PoMoBalance(const RateMatrix_PoMoBalance& m) : TimeReversibleRateMatrix( m ),
    N( m.N ),
    pi( m.pi ),
    rho( m.rho ),
    sigma( m.sigma ),
    beta( m.beta )

{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_PoMoBalance::~RateMatrix_PoMoBalance(void)
{
    
    delete eigen_system;
}


RateMatrix_PoMoBalance& RateMatrix_PoMoBalance::operator=(const RateMatrix_PoMoBalance &r)
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
        sigma               = r.sigma;
        beta                = r.beta;

        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


RateMatrix_PoMoBalance& RateMatrix_PoMoBalance::assign(const Assignable &m)
{
    
    const RateMatrix_PoMoBalance *rm = dynamic_cast<const RateMatrix_PoMoBalance*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
    
}



/** Do precalculations on eigenvectors */
void RateMatrix_PoMoBalance::calculateCijk(void)
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
void RateMatrix_PoMoBalance::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_PoMoBalance* RateMatrix_PoMoBalance::clone( void ) const
{
    return new RateMatrix_PoMoBalance( *this );
}


/*
To build thie PoMoBalance rate matrix we need 4 vectors that organize 8 parameters describing genetic drift,
mutation, allele flow and selection between population 1 and 2:

effective population size:  nu     = ( N_1         , N_2 )
mutation:                   mu     = ( mu_{Aa}     , mu_{aA} )
allele flow:                lambda = ( lambda_{12} , lambda_{21} )
allelic selection:          sigma  = ( sigma_a     , sigma_A )
*/

void RateMatrix_PoMoBalance::computeOffDiagonal( void )
{
    
    
    MatrixReal& m = *the_rate_matrix;
    // frequency of allele A in population 1 (n1) and 2 (n2)
    // create the combine state-space between the two populations
    // {n1A,(N1-n1)a}{n2A,(N2-n2)a}

    for (int i=0; i<num_states; i++){

        for (int j=0; j<num_states; j++){

            m[i][j] = 0.0;
        
        }

    }


    // Mutations

    m[0][N+2]   = rho[0]*pi[1];
    m[0][2*N+1] = rho[1]*pi[2];
    m[0][3*N]   = rho[2]*pi[3];
    m[1][4]     = rho[0]*pi[0];
    m[1][4*N-1] = rho[3]*pi[2];
    m[1][5*N-2] = rho[4]*pi[3];
    m[2][N+3]   = rho[1]*pi[0];
    m[2][3*N+1] = rho[3]*pi[1];
    m[2][6*N-3] = rho[5]*pi[3];
    m[3][2*N+2] = rho[2]*pi[0];
    m[3][4*N]   = rho[4]*pi[1];
    m[3][5*N-1] = rho[5]*pi[2];

    // Fixations

    m[4][1]     = (1.0+sigma[1])*(N-1.0)/(N*beta);
    m[N+2][0]   = (1.0+sigma[0])*(N-1.0)/(N*beta);
    m[N+3][2]   = (1.0+sigma[2])*(N-1.0)/(N*beta);
    m[2*N+1][0] = (1.0+sigma[0])*(N-1.0)/(N*beta);
    m[2*N+2][3] = (1.0+sigma[3])*(N-1.0)/(N*beta);
    m[3*N+1][2] = (1.0+sigma[2])*(N-1.0)/(N*beta);
    m[4*N][3]   = (1.0+sigma[3])*(N-1.0)/(N*beta);
    m[5*N-1][3] = (1.0+sigma[3])*(N-1.0)/(N*beta);
    m[3*N][0]   = (1.0+sigma[0])*(N-1.0)/(N*beta);
    m[4*N-1][1] = (1.0+sigma[1])*(N-1.0)/(N*beta);
    m[5*N-2][1] = (1.0+sigma[1])*(N-1.0)/(N*beta);
    m[6*N-3][2] = (1.0+sigma[2])*(N-1.0)/(N*beta);

    // reamining polymorphic sites

    if (N > 2){

        for (int n=1; n<(N-1); n++) {

        m[4+n-1][4+n-1+1]             = n*(N-n)*(1.0+sigma[0])*midpoint_beta(n,1)/N;
        m[N+3-1+1-n][N+3-1+1-n+1]     = n*(N-n)*(1.0+sigma[1])*midpoint_beta(N-n,0)/N;

        m[N+n+2+n-1][N+n+2+n-1]       = n*(N-n)*(1.0+sigma[0])*midpoint_beta(n,1)/N;
        m[2*N+2-1+1-n][2*N+2-1-n]     = n*(N-n)*(1.0+sigma[2])*midpoint_beta(N-n,0)/N;

        m[2*N+2+n-1][2*N+2+n]         = n*(N-n)*(1.0+sigma[0])*midpoint_beta(n,1)/N;
        m[3*N+1-1+1-n][3*N+1-1+1-n-1] = n*(N-n)*(1.0+sigma[3])*midpoint_beta(N-n,0)/N;

        m[3*N+1+n-1][3*N+1+n]         = n*(N-n)*(1.0+sigma[1])*midpoint_beta(n,1)/N;
        m[4*N-1+1-n][4*N-1+1-n-1]     = n*(N-n)*(1.0+sigma[2])*midpoint_beta(N-n,0)/N;

        m[4*N+n-1][4*N+n]             = n*(N-n)*(1.0+sigma[1])*midpoint_beta(n,1)/N;
        m[5*N-1-1+1-n][5*N-1-1+1-n-1] = n*(N-n)*(1.0+sigma[3])*midpoint_beta(N-n,0)/N;

        m[5*N-1+n-1][5*N-1+n]         = n*(N-n)*(1.0+sigma[2])*midpoint_beta(n,1)/N;
        m[6*N-3-1+1-n][6*N-3-1+1-n-1] = n*(N-n)*(1.0+sigma[3])*midpoint_beta(N-n,0)/N;

      }

    }


    needs_update = true;
}




double RateMatrix_PoMoBalance::midpoint_beta( double n, int indicator ) const
{

  // allele frequency n is in the midpoint
  // will only happen for N even
  if ( n == N/2) {
    
    if (indicator == 0) {
        return 1/beta;
    } else {
        return 1/beta;
    }

  // allele frequency n is lower than the midpoint 
  } else if ( n < N/2 ) {

    if (indicator == 0) {
        return 1/beta;
    } else {
        return beta;
    }

  // allele frequency n is higher than the midpoint 
  } else {

    if (indicator == 0) {
        return beta;
    } else {
        return 1/beta;
    }

  }


}


/** Calculate the transition probabilities for the real case */
void RateMatrix_PoMoBalance::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_PoMoBalance::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_PoMoBalance::setN( double n )
{
    
    N = n;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMoBalance::setPi(const std::vector<double> &p )
{
    pi = p;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_PoMoBalance::setRho( const std::vector<double> &r )
{
    rho = r;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoBalance::setSigma( const std::vector<double> &s )
{
    sigma = s;
    
    // set flags
    needs_update = true;
}

void RateMatrix_PoMoBalance::setBeta( double b )
{
    beta = b;
    
    // set flags
    needs_update = true;
}


/** Update the eigen system */
void RateMatrix_PoMoBalance::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_PoMoBalance::update( void )
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



