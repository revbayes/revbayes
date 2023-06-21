#include "RateMatrix_PoMo2N.h"

#include <assert.h>
#include <cstddef>

#include "RbException.h"
#include "Assignable.h"
#include "Cloneable.h"
#include "EigenSystem.h"
#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "RbMathCombinatorialFunctions.h"
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
RateMatrix_PoMo2N::RateMatrix_PoMo2N(long num_states, long in_n, bool mu_corr, bool d_corr )  :
AbstractRateMatrix( num_states ), 
N( in_n ),
mu( 2 , 0.01 ),
phi( 2 , 1.0 ),
use_mutation_correction( mu_corr ),
use_drift_correction( d_corr ),
N_eff( 1000.0 )
{
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    harmonic_number_M = RbMath::fastHarmonicNumber(N-1);
    
  
    update();
}

/** Copy constructor */
RateMatrix_PoMo2N::RateMatrix_PoMo2N(const RateMatrix_PoMo2N& m) : 
AbstractRateMatrix( m ), 
N( m.N ),
mu( m.mu ),
phi( m.phi ),
use_mutation_correction( m.use_mutation_correction ),
use_drift_correction( m.use_drift_correction ),
N_eff( m.N_eff),
harmonic_number_M( m.harmonic_number_M ),
harmonic_number_N( m.harmonic_number_N )
{

    
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
    
}


/** Destructor */
RateMatrix_PoMo2N::~RateMatrix_PoMo2N(void)
{
    
    delete theEigenSystem;
}


RateMatrix_PoMo2N& RateMatrix_PoMo2N::operator=(const RateMatrix_PoMo2N &r)
{

    if (this != &r)
    {
        AbstractRateMatrix::operator=( r );
        N                   = r.N;
        mu                  = r.mu;
        phi                 = r.phi;
        
        use_mutation_correction = r.use_mutation_correction;
        use_drift_correction    = r.use_drift_correction;
        
        N_eff               = r.N_eff;
        harmonic_number_M   = r.harmonic_number_M;
        harmonic_number_N   = r.harmonic_number_N;
        
        
        delete theEigenSystem;
        
        theEigenSystem       = new EigenSystem( *r.theEigenSystem );
        c_ijk                = r.c_ijk;
        cc_ijk               = r.cc_ijk;
        
        theEigenSystem->setRateMatrixPtr(the_rate_matrix);
        
  }

  return *this;
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
    throw RbException("Missing implementation of average rate in PoMo2N.");
    
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

    // populate rate matrix with 0.0
    // **total waste of time with sparse matrices like pomos**
    for (int i=0; i< num_states; i++)
    {
        for (int j=0; j< num_states; j++)
        {
            m[i][j] = 0.0;
        }
    }

    // set the diagonal vector
    std::vector< double > diagonal(num_states,0.0);
    
    harmonic_number_N = RbMath::fastHarmonicNumber(N_eff-1);
    
    double rescaling_factor_to_effective_population = (N_eff) / (N)   *    (harmonic_number_N/harmonic_number_M);
    if ( use_drift_correction == false )
    {
        rescaling_factor_to_effective_population = 1.0;
    }
    
    // mutations
    if ( use_mutation_correction == false )
    {
        m[0][2]  = N*mu[0];  //{Na0} -> {(N-1)a0,1a1}
        m[1][N]  = N*mu[1];  //{Na1} -> {1a0,(N-1)a1}
    }
    else
    {
//        m[0][2]  = N_eff*mu[0];  //{Na0} -> {(N-1)a0,1a1}
//        m[1][N]  = N_eff*mu[1];  //{Na1} -> {1a0,(N-1)a1}
//        m[0][2]  = N*mu[0]*harmonic_number_M/harmonic_number_N;  //{Na0} -> {(N-1)a0,1a1}
//        m[1][N]  = N*mu[1]*harmonic_number_M/harmonic_number_N;  //{Na1} -> {1a0,(N-1)a1}
        m[0][2]  = N*mu[0]*N/N_eff;  //{Na0} -> {(N-1)a0,1a1}
        m[1][N]  = N*mu[1]*N/N_eff;  //{Na1} -> {1a0,(N-1)a1}
    }
    diagonal[0] = m[0][2];
    diagonal[1] = m[1][N];

       

    // the pomo rate matrix is entirely defined by fixations and mutations if N=2
    if (N>2)
    {

        //frequency shifts from singletons

//        m[2][3]   = (N-1.0)*phi[1]/((N-1)*phi[0]+phi[1]);  //{(N-1)a0,1a1} -> {(N-2)a0,2a1}
//        m[N][N-1] = (N-1.0)*phi[0]/((N-1)*phi[1]+phi[0]);  //{1a0,(N-1)a1} -> {2a0,(N-2)a1}
//        m[2][3]   = (N-1.0)*phi[1]/((N-1)*phi[0]+phi[1]) / rescaling_factor_to_effective_population;  //{(N-1)a0,1a1} -> {(N-2)a0,2a1}
//        m[N][N-1] = (N-1.0)*phi[0]/((N-1)*phi[1]+phi[0]) / rescaling_factor_to_effective_population;  //{1a0,(N-1)a1} -> {2a0,(N-2)a1}
        m[2][3]   = ((N-1.0)/double(N)) / rescaling_factor_to_effective_population;  //{(N-1)a0,1a1} -> {(N-2)a0,2a1}
        m[N][N-1] = ((N-1.0)/double(N)) / rescaling_factor_to_effective_population;  //{1a0,(N-1)a1} -> {2a0,(N-2)a1}
        diagonal[2] += m[2][3] ;
        diagonal[N] += m[N][N-1];

        //frequency shifts for all the other polymorphic states
        if (N>3)
        {

            //polymorphic states are populated
            for (int n=2; n<(N-1); n++)
            {

                //ai aj
//                m[n+1][2+n]  = n*(N-n)*phi[1]/(n*phi[1]+(N-n)*phi[0]); //{naj,(N-n)ai} -> {(n+1)aj,(N-n-1)ai}
//                m[n+1][n]    = n*(N-n)*phi[0]/(n*phi[1]+(N-n)*phi[0]); //{naj,(N-n)ai} -> {(n-1)aj,(N-n+1)ai}
//                m[n+1][2+n]  = (n*(N-n)*phi[1]/(n*phi[1]+(N-n)*phi[0])) / rescaling_factor_to_effective_population; //{naj,(N-n)ai} -> {(n+1)aj,(N-n-1)ai}
//                m[n+1][n]    = (n*(N-n)*phi[0]/(n*phi[1]+(N-n)*phi[0])) / rescaling_factor_to_effective_population; //{naj,(N-n)ai} -> {(n-1)aj,(N-n+1)ai}
                m[n+1][2+n]  = (n*(N-n)/double(N)) / rescaling_factor_to_effective_population; //{naj,(N-n)ai} -> {(n+1)aj,(N-n-1)ai}
                m[n+1][n]    = (n*(N-n)/double(N)) / rescaling_factor_to_effective_population; //{naj,(N-n)ai} -> {(n-1)aj,(N-n+1)ai}
                diagonal[n+1] += m[n+1][2+n] + m[n+1][n];

            }

        }

    }


    // set the diagonal values and calculate the reciprocal of the average rate
    for (int j=0; j< num_states; j++)
    {
        m[j][j] = -diagonal[j];
    }
  
    
}



/** Do precalculations on eigenvectors */
void RateMatrix_PoMo2N::calculateCijk(void)
{
    
    if ( theEigenSystem->isComplex() == false )
    {
        // real case
        const MatrixReal& ev  = theEigenSystem->getEigenvectors();
        const MatrixReal& iev = theEigenSystem->getInverseEigenvectors();
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
        const MatrixComplex& cev  = theEigenSystem->getComplexEigenvectors();
        const MatrixComplex& ciev = theEigenSystem->getComplexInverseEigenvectors();
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
void RateMatrix_PoMo2N::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{

//  // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
//  // Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
//  double t = rate * (startAge - endAge);
//  exponentiateMatrixByScalingAndSquaring(t, P );

    
    double t = rate * (startAge - endAge);
    if ( theEigenSystem->isComplex() == false )
    {
        tiProbsEigens(t, P);
    }
    else
    {
        tiProbsComplexEigens(t, P);
    }
}


RateMatrix_PoMo2N* RateMatrix_PoMo2N::clone( void ) const
{
    return new RateMatrix_PoMo2N( *this );
}


std::vector<double> RateMatrix_PoMo2N::getStationaryFrequencies( void ) const
{
    return calculateStationaryFrequencies();

//    return stationaryVector;
}


void RateMatrix_PoMo2N::setNeff( double ni )
{
    N_eff = ni;
    
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


/** Calculate the transition probabilities for the real case */
void RateMatrix_PoMo2N::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValue = theEigenSystem->getRealEigenvalues();
    
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
        double rowsum = 0.0;
        for (size_t j=0; j<num_states; j++, ++p)
        {
            double sum = 0.0;
            for (size_t s=0; s<num_states; s++)
            {
                sum += (*ptr++) * eigValExp[s];
            }
            
            sum = (sum < 0.0) ? 0.0 : sum;
            rowsum += sum;
            (*p) = sum;
        }

        // Normalize transition probabilities for row to sum to 1.0
        double* p2 = p - num_states;
        for (size_t j=0; j<num_states; j++, ++p2)
            *p2 /= rowsum;
    }
}



/** Calculate the transition probabilities for the complex case */
void RateMatrix_PoMo2N::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
{
    
    // get a reference to the eigenvalues
    const std::vector<double>& eigenValueReal = theEigenSystem->getRealEigenvalues();
    const std::vector<double>& eigenValueComp = theEigenSystem->getImagEigenvalues();
    
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
        double rowsum = 0.0;
        for (size_t j=0; j<num_states; j++)
        {
            std::complex<double> sum = std::complex<double>(0.0, 0.0);
            for (size_t s=0; s<num_states; s++)
                sum += (*ptr++) * ceigValExp[s];

            double real_sum = (sum.real() < 0.0) ? 0.0 : sum.real();
            P[i][j] = real_sum;
            rowsum += real_sum;
        }
        // Normalize transition probabilities for row to sum to 1.0
        for (size_t j=0; j<num_states; j++)
            P[i][j] /= rowsum;
    }
}


/** Update the eigen system */
void RateMatrix_PoMo2N::updateEigenSystem(void)
{
    
    theEigenSystem->update();
    calculateCijk();
    
}



void RateMatrix_PoMo2N::update( void )
{

    if ( needs_update )
    {

        buildRateMatrix();
        
        // now update the eigensystem
        updateEigenSystem();

        // clean flags
        needs_update = false;
    }
}
