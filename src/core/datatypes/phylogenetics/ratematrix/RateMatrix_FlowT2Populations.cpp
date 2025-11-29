#include "CodonState.h"
#include "EigenSystem.h"
#include "MatrixComplex.h"
#include "MatrixReal.h"
#include "RateMatrix_FlowT2Populations.h"
#include "RbException.h"
#include "RbMathMatrix.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <string>
#include <iomanip>

using namespace RevBayesCore;




/** Construct rate matrix with n states */
RateMatrix_FlowT2Populations::RateMatrix_FlowT2Populations( size_t ss, const std::vector<double>& in_nu ) : TimeReversibleRateMatrix( ss ),
    nu( in_nu ),
    mu( 2, 0.001),
    lambda( 2, 0.0 ),
    sigma(  2, 0.0 )
{
    
    eigen_system       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(num_states * num_states * num_states);
    cc_ijk.resize(num_states * num_states * num_states);
    
    update();
}


/** Copy constructor */
RateMatrix_FlowT2Populations::RateMatrix_FlowT2Populations(const RateMatrix_FlowT2Populations& m) : TimeReversibleRateMatrix( m ),
    nu( m.nu ),
    mu( m.mu ),
    lambda( m.lambda ),
    sigma( m.sigma )
{
    
    eigen_system        = new EigenSystem( *m.eigen_system );
    c_ijk               = m.c_ijk;
    cc_ijk              = m.cc_ijk;
    
    eigen_system->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_FlowT2Populations::~RateMatrix_FlowT2Populations(void)
{
    
    delete eigen_system;
}


RateMatrix_FlowT2Populations& RateMatrix_FlowT2Populations::operator=(const RateMatrix_FlowT2Populations &r)
{
    
    if (this != &r)
    {
        TimeReversibleRateMatrix::operator=( r );
        
        delete eigen_system;
        
        eigen_system        = new EigenSystem( *r.eigen_system );
        c_ijk               = r.c_ijk;
        cc_ijk              = r.cc_ijk;
        nu                  = r.nu;
        mu                  = r.mu;
        lambda              = r.lambda;
        sigma               = r.sigma;

        
        eigen_system->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/** Do precalculations on eigenvectors */
void RateMatrix_FlowT2Populations::calculateCijk(void)
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
void RateMatrix_FlowT2Populations::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
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


RateMatrix_FlowT2Populations* RateMatrix_FlowT2Populations::clone( void ) const
{
    return new RateMatrix_FlowT2Populations( *this );
}


/*
To build thie FlowT2Populations rate matrix we need 4 vectors that organize 8 parameters describing genetic drift,
mutation, allele flow and selection between population 1 and 2:

effective population size:  nu     = ( N_1         , N_2 )
mutation:                   mu     = ( mu_{Aa}     , mu_{aA} )
allele flow:                lambda = ( lambda_{12} , lambda_{21} )
allelic selection:          sigma  = ( sigma_A     , sigma_a )
*/

void RateMatrix_FlowT2Populations::computeOffDiagonal( void )
{
    
    
    MatrixReal& m = *the_rate_matrix;
    std::cout << "matrix was created";

    // frequency of allele A in population 1 (n1) and 2 (n2)
    // create the combine state-space between the two populations
    // {n1A,(N1-n1)a}{n2A,(N2-n2)a}

    std::vector<int> frequency_1 = std::vector<int>(num_states,0);

    std::vector<int> frequency_2 = std::vector<int>(num_states,0);

    int index = 0;

    for (int i=0; i<nu[0]+1; ++i){

        for (int j=0; j<nu[1]+1; ++j){

            frequency_1[index] = i;

            frequency_2[index] = j;

            index=index+1;

        }
    }

    for (int i=0; i<num_states; i++ ) {

        for (int j=0; j<num_states; j++ ) {

            // some quantities of interest:
            // initial frequency of allele A in both populations
            // allele frequency change in both populations

            int initial_frequency_1 = frequency_1[i];

            int initial_frequency_2 = frequency_2[i];

            int change_frequency_1 = frequency_1[j] - frequency_1[i];

            int change_frequency_2 = frequency_2[j] - frequency_2[i];

            if ( (abs(change_frequency_1) + abs(change_frequency_2)) != 1 ) {

                // populating all state changes that require more the one allele frequency shift with 0.0
                // assigns diagonal elementes to 0.0
                m[i][j] = 0.0;

            } else {

                if ( initial_frequency_1 == nu[0] && change_frequency_1 == -1 ) {

                    // a emerges in population 1
                    // mutation A to a OR allele flow from population 2
                    // note that the second term is null if a is absent in population 2 (initial_frequency_2 = nu[1])
                    m[i][j] = mu[0]/nu[0] + (nu[1]-initial_frequency_2)*lambda[1]/nu[1];

                } else if ( initial_frequency_2 == nu[1] && change_frequency_2 == -1 ) {

                    // a emerges in population 2
                    // mutation A to a OR allele flow from population 1
                    // note that the second term is null if a is absent in population 1 (initial_frequency_1 = nu[0])
                    m[i][j] = mu[0]/nu[1] + (nu[0]-initial_frequency_1)*lambda[0]/nu[0];

                } else if ( initial_frequency_1 == 0 && change_frequency_1 == 1 ) {

                    // A emerges in population 1
                    // mutation a to A OR allele flow from population 2
                    // note that the second term is null if A is absent in population 2 (initial_frequency_2 = 0)
                    m[i][j] = mu[1]/nu[0] + initial_frequency_2*lambda[1]/nu[1];

                } else if ( initial_frequency_2 == 0 && change_frequency_2 == 1  ) {

                    // A emerges in population 2
                    // mutation a to A OR allele flow from population 1
                    // note that the second term is null if A is absent in population 1 (initial_frequency_1 = 0)
                    m[i][j] = mu[1]/nu[1] + initial_frequency_1*lambda[0]/nu[0];

                } else if ( initial_frequency_1 > 0 && initial_frequency_1 < nu[0] && change_frequency_1 == 1 ) {

                    // A increases in frequency in population 1
                    // drift and selection OR allele flow from population 2
                    // note that the second term is null if A is absent in population 2 (initial_frequency_2 = 0)
                    m[i][j] = initial_frequency_1*(nu[0]-initial_frequency_1)*sigma[1]/(nu[0]*nu[0]) + (nu[0]-initial_frequency_1)*initial_frequency_2*lambda[1]/(nu[0]*nu[1]);

                } else if ( initial_frequency_1 > 0 && initial_frequency_1 < nu[0] && change_frequency_1 == -1 ) {

                    // a increases in frequency in population 1
                    // drift and selection OR allele flow from population 2
                    // note that the second term is null if a is absent in population 2 (initial_frequency_2 = nu[1])
                    m[i][j] = initial_frequency_1*(nu[0]-initial_frequency_1)*sigma[0]/(nu[0]*nu[0]) + initial_frequency_1*(nu[1]-initial_frequency_2)*lambda[1]/(nu[0]*nu[1]);

                } else if ( initial_frequency_2 > 0 && initial_frequency_2 < nu[1] && change_frequency_2 == 1 ) {

                    // A increases in frequency in population 2
                    // drift and selection OR allele flow from population 1
                    // note that the second term is null if A is absent in population 1 (initial_frequency_1 = 0)
                    m[i][j] = initial_frequency_2*(nu[1]-initial_frequency_2)*sigma[1]/(nu[1]*nu[1]) + (nu[1]-initial_frequency_2)*initial_frequency_1*lambda[0]/(nu[1]*nu[0]);


                } else if ( initial_frequency_2 > 0 && initial_frequency_2 < nu[1] && change_frequency_2 == -1 ) {

                    // a increases in frequency in population 2
                    // drift and selection OR allele flow from population 1
                    // note that the second term is null if a is absent in population 1 (initial_frequency_1 = nu[0])
                    m[i][j] = initial_frequency_2*(nu[1]-initial_frequency_2)*sigma[0]/(nu[1]*nu[1]) + initial_frequency_2*(nu[0]-initial_frequency_1)*lambda[0]/(nu[1]*nu[0]);

                } else {

                     RbException("Something went wrong!");
                
                }

            }

        }
    }
    
    // set flags
    needs_update = true;
}



/** Calculate the transition probabilities for the real case */
void RateMatrix_FlowT2Populations::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_FlowT2Populations::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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


void RateMatrix_FlowT2Populations::setNu(const std::vector<double> &n )
{
    
    nu = n;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_FlowT2Populations::setMu(const std::vector<double> &m )
{
    
    mu = m;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_FlowT2Populations::setLambda( const std::vector<double> &l )
{
    lambda = l;
    
    // set flags
    needs_update = true;
}

void RateMatrix_FlowT2Populations::setSigma( const std::vector<double> &s )
{
    sigma = s;
    
    // set flags
    needs_update = true;
}



/** Update the eigen system */
void RateMatrix_FlowT2Populations::updateEigenSystem(void)
{
    
    eigen_system->update();
    calculateCijk();
    
}


void RateMatrix_FlowT2Populations::update( void )
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



