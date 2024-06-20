#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_revPoMoBalanceKN.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_revPoMoBalanceKN::RateMatrix_revPoMoBalanceKN( long num_states, long in_k, long in_n, long in_nex) : TimeReversibleRateMatrix( num_states ),
    K( in_k ),
    N( in_n ),
    pi( in_k, 1.0/in_k),
    rho( in_nex, 0.01 ),
    phi(  in_k, 1.0 ),
    beta( in_nex, 1.0 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_revPoMoBalanceKN::RateMatrix_revPoMoBalanceKN(const RateMatrix_revPoMoBalanceKN& m) : TimeReversibleRateMatrix( m ),
    K( m.K ),
    N( m.N ),
    pi( m.pi ),
    rho( m.rho ),
    phi( m.phi ),
    beta( m.beta )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_revPoMoBalanceKN::~RateMatrix_revPoMoBalanceKN(void)
{
    
    delete eigen_system;
}


RateMatrix_revPoMoBalanceKN& RateMatrix_revPoMoBalanceKN::operator=(const RateMatrix_revPoMoBalanceKN &r)
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
        beta                = r.beta;

        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}



/** Do precalculations on eigenvectors */
void RateMatrix_revPoMoBalanceKN::calculateCijk(void)
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
void RateMatrix_revPoMoBalanceKN::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_revPoMoBalanceKN* RateMatrix_revPoMoBalanceKN::clone( void ) const
{
    return new RateMatrix_revPoMoBalanceKN( *this );
}



void RateMatrix_revPoMoBalanceKN::computeOffDiagonal( void )
{
    
    
  MatrixReal& m = *the_rate_matrix;

  //populating the rate matrix

  //populate rate matrix with 0.0
  // **total waste of time with sparse matrices like pomos**
  for (int i=0; i< num_states; i++){
    for (int j=0; j< num_states; j++){
      m[i][j] = 0.0;        
    }
  }


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

    // Original by Rui Borges (without time in numbers of generations)
    //reciprocal of the population size
    double rN = 1.0 / N;
    //first edge
    int E = 0;

    // Preferred frequencies are in the middle for the reversible case
    std::vector<int> B(beta.size(),round(0.5*N));

    for (int i=0; i<K; i++) {
        for (int j = i + 1; j < K; j++) {
            // reversible mutation rates
            m[i][K + E * N - E] = pi[j]*rho[E];              //{Nai} -> {(N-1)ai,1aj}
            m[j][K + (E + 1) * (N - 1) - 1] = pi[i]*rho[E];  //{Naj} -> {1ai,(N-1)aj}

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

  needs_update = true;
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_revPoMoBalanceKN::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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

std::vector<double> RateMatrix_revPoMoBalanceKN::getStationaryFrequencies( void ) const
{
//    std::cout<<"Frequencies "<<pi<<std::endl;
//    std::cout<<"Selection coefficients "<<phi<<std::endl;
//    std::cout<<"Population size "<<N<<std::endl;

    // calculating the normalization constant and stationary weights
    int n_states = K+(K*K-K)*(N-1)/2;
    std::vector<double> stationary_freqs(n_states,0.0);

    // Preferred frequencies are in the middle for the reversible case
    std::vector<int> B(beta.size(),round(0.5*N));

    // normalization constant
    double nc = 0.0;

    //calculating the normalization constant psi_0 and psi_N for no selection
    for (int i=0; i<K; i++) {
        stationary_freqs[i] = pi[i]*pow(phi[i],N-1);
        nc += stationary_freqs[i];
    }

    // now for the polymorphic states
    //first edge
    int index = K;
    int E = 0;
    double drift_coefficient;
    // psi_n
    for (int i=0; i<K; i++){
        for (int j=i+1; j<K; j++){
            for (int n=1; n<N; n++) {

                // sum weight
                drift_coefficient = 1.0*N/(n*(N-n));
                stationary_freqs[index] = drift_coefficient*(pi[i]*pi[j]*rho[E]*pow(phi[j],n-1)*pow(phi[i],N-n-1)*pow(beta[E], B[E]-abs(n-B[E])-1));
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

/** Calculate the transition probabilities for the complex case */
void RateMatrix_revPoMoBalanceKN::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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

void RateMatrix_revPoMoBalanceKN::setK( long & na )
{
    K = na;

    // set flags
    needs_update = true;

}

void RateMatrix_revPoMoBalanceKN::setN( long & ni )
{
    
    N = ni;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_revPoMoBalanceKN::setPi(const std::vector<double> &p )
{
    pi = p;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_revPoMoBalanceKN::setRho( const std::vector<double> &r )
{
    rho = r;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoBalanceKN::setPhi( const std::vector<double> &s )
{
    phi = s;
    
    // set flags
    needs_update = true;
}

void RateMatrix_revPoMoBalanceKN::setBeta( const std::vector<double> &b )
{
    beta = b;
    
    // set flags
    needs_update = true;
}

/** Update the eigen system */
void RateMatrix_revPoMoBalanceKN::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_revPoMoBalanceKN::update( void )
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



