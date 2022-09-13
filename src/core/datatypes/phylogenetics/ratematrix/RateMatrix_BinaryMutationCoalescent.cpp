#include "RateMatrix_BinaryMutationCoalescent.h"
#include "EigenSystem.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states, virtual population size, mutation rates, selection coefficients */
RateMatrix_BinaryMutationCoalescent::RateMatrix_BinaryMutationCoalescent(size_t n) : AbstractRateMatrix( n * (n+3) / 2.0 ),
    num_ind( n ),
    matrix_size( (n*(n+3)/2.0) ),
    mu( 1.0 ),
    Ne( 1.0 )
{
    
    
    theEigenSystem       = new EigenSystem(the_rate_matrix);
    c_ijk.resize(matrix_size * matrix_size * matrix_size);
    cc_ijk.resize(matrix_size * matrix_size * matrix_size);
    
    update();
}


/** Copy constructor */
RateMatrix_BinaryMutationCoalescent::RateMatrix_BinaryMutationCoalescent(const RateMatrix_BinaryMutationCoalescent& m) : AbstractRateMatrix( m ),
    num_ind( m.num_ind ),
    matrix_size( m.matrix_size ),
    mu( m.mu ),
    Ne( m.Ne )
{
    
    theEigenSystem       = new EigenSystem( *m.theEigenSystem );
    c_ijk                = m.c_ijk;
    cc_ijk               = m.cc_ijk;
    
    theEigenSystem->setRateMatrixPtr(the_rate_matrix);
}


/** Destructor */
RateMatrix_BinaryMutationCoalescent::~RateMatrix_BinaryMutationCoalescent(void)
{
    
    delete theEigenSystem;
}


RateMatrix_BinaryMutationCoalescent& RateMatrix_BinaryMutationCoalescent::operator=(const RateMatrix_BinaryMutationCoalescent &r)
{
    
    if (this != &r)
    {
        AbstractRateMatrix::operator=( r );
        
        num_ind         = r.num_ind;
        matrix_size     = r.matrix_size;
        mu              = r.mu;
        Ne              = r.Ne;
        
        delete theEigenSystem;
        
        theEigenSystem       = new EigenSystem( *r.theEigenSystem );
        c_ijk                = r.c_ijk;
        cc_ijk               = r.cc_ijk;
        
        theEigenSystem->setRateMatrixPtr(the_rate_matrix);
    }
    
    return *this;
}


/**
 * Assign the value of m to this instance. This function is our mechanism to call the assignment operator.
 *
 *
 */
RateMatrix_BinaryMutationCoalescent& RateMatrix_BinaryMutationCoalescent::assign(const Assignable &m)
{
    
    const RateMatrix_BinaryMutationCoalescent *rm = dynamic_cast<const RateMatrix_BinaryMutationCoalescent*>(&m);
    if ( rm != NULL )
    {
        return operator=(*rm);
    }
    else
    {
        throw RbException("Could not assign rate matrix.");
    }
}

double RateMatrix_BinaryMutationCoalescent::averageRate(void) const
{
    return 1.0;
}


void RateMatrix_BinaryMutationCoalescent::buildRateMatrix(void)
{

    // calculate the transition probabilities
    for (size_t i=0; i< matrix_size; i++)
    {
        //The first 4 states are the monomorphic states; we can't directly change from one into another one
        for (size_t j=0; j< matrix_size; j++)
        {
            (*the_rate_matrix)[i][j] = 0.0;
        }
        
    }
    
    size_t row_index = 0;
    for ( size_t num_total_ind=num_ind; num_total_ind>=1; --num_total_ind )
    {
        for (size_t num_derived_ind=0; num_derived_ind<=num_total_ind; ++num_derived_ind)
        {
            
            if ( num_total_ind > num_derived_ind )
            {
                // we can have a mutation from the 0 to the 1 state
                (*the_rate_matrix)[row_index][row_index+1] = (num_total_ind-num_derived_ind)*mu;
            }
            
            if ( num_derived_ind > 0 )
            {
                // we can have a mutation from the 1 to the 0 state
                (*the_rate_matrix)[row_index][row_index-1] = num_derived_ind*mu;
            }
            
            if ( (num_total_ind-num_derived_ind) > 1 )
            {
                // we can have a coalescent event in the ancestral state
                (*the_rate_matrix)[row_index][row_index+num_total_ind+1] = (num_total_ind-num_derived_ind)*(num_total_ind-num_derived_ind-1)/(2.0*2.0*Ne);
            }
            
            if ( num_derived_ind > 1 )
            {
                // we can have a coalescent event in the ancestral state
                (*the_rate_matrix)[row_index][row_index+num_total_ind] = num_derived_ind*(num_derived_ind-1)/(2.0*2.0*Ne);
            }
            
            ++row_index;
        }
    }
    
    // set the diagonal values
    setDiagonal();
    
    
    
    // rescale
    //rescaleToAverageRate( 1.0 );
    
    
}



/** Do precalculations on eigenvectors */
void RateMatrix_BinaryMutationCoalescent::calculateCijk(void)
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
void RateMatrix_BinaryMutationCoalescent::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    // Now the instantaneous rate matrix has been filled up entirely.
    // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
    double t = rate * (startAge - endAge);
//    computeExponentialMatrixByRepeatedSquaring(t, P);
    
    if ( theEigenSystem->isComplex() == false )
    {
        tiProbsEigens(t, P);
    }
    else
    {
        tiProbsComplexEigens(t, P);
    }
    
    return;
}

void RateMatrix_BinaryMutationCoalescent::computeExponentialMatrixByRepeatedSquaring(double t,  TransitionProbabilityMatrix& P ) const
{
    //We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
    //Ideally one should dynamically decide how many squarings are necessary.
    //For the moment, we arbitrarily do 10 such squarings, as it seems to perform well in practice (N. Lartillot, personal communication).
    //first, multiply the matrix by the right scalar
    //2^10 = 1024
    double tOver2s = t/(1024);
    for ( size_t i = 0; i < matrix_size; i++ )
    {
        for ( size_t j = 0; j < matrix_size; j++ )
        {
            P[i][j] = (*the_rate_matrix)[i][j] * tOver2s;
        }
    }
    
    //Add the identity matrix:
    for ( size_t i = 0; i < matrix_size; i++ )
    {
        P[i][i] += 1;
    }
    //Now we can do the multiplications
    TransitionProbabilityMatrix P2 (matrix_size);
    squareMatrix (P, P2); //P2 at power 2
    squareMatrix (P2, P); //P at power 4
    squareMatrix (P, P2); //P2 at power 8
    squareMatrix (P2, P); //P at power 16
    squareMatrix (P, P2); //P2 at power 32
    squareMatrix (P2, P); //P at power 64
    squareMatrix (P, P2); //P2 at power 128
    squareMatrix (P2, P); //P at power 256
    squareMatrix (P, P2); //P2 at power 512
    squareMatrix (P2, P); //P at power 1024
    
    return;
}

inline void RateMatrix_BinaryMutationCoalescent::squareMatrix( TransitionProbabilityMatrix& P,  TransitionProbabilityMatrix& P2) const
{
    //Could probably use boost::ublas here, for the moment we do it ourselves.
    for ( size_t i = 0; i < matrix_size; i++ )
    {
        for ( size_t j = 0; j < matrix_size; j++ )
        {
            P2.getElement ( i, j ) = 0;
            for ( size_t k = 0; k < matrix_size; k++ )
            {
                P2.getElement ( i, j ) += P.getElement ( i, k ) * P.getElement ( k, j );
                
            }
        }
    }
}



RateMatrix_BinaryMutationCoalescent* RateMatrix_BinaryMutationCoalescent::clone( void ) const
{
    return new RateMatrix_BinaryMutationCoalescent( *this );
}

std::vector<double> RateMatrix_BinaryMutationCoalescent::getStationaryFrequencies( void ) const
{
    
    return stationary_freqs;
}


/** Calculate the transition probabilities for the real case */
void RateMatrix_BinaryMutationCoalescent::tiProbsEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_BinaryMutationCoalescent::tiProbsComplexEigens(double t, TransitionProbabilityMatrix& P) const
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
void RateMatrix_BinaryMutationCoalescent::updateEigenSystem(void)
{
    
    theEigenSystem->update();
    calculateCijk();
    
}


void RateMatrix_BinaryMutationCoalescent::update( void )
{
    
    if ( needs_update || true )
    {
        buildRateMatrix();
        
        // now update the eigensystem
        updateEigenSystem();
        
        // clean flags
        needs_update = false;
    }
}


void RateMatrix_BinaryMutationCoalescent::setMutationRate(double m)
{
    
    mu = m;
    needs_update = true;
    
}


void RateMatrix_BinaryMutationCoalescent::setEffectivePopulationSize(double n)
{

    Ne = n;
    needs_update = true;
    
}
