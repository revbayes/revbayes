#include "RateMatrix_BinaryMutationCoalescent.h"
#include "MatrixReal.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "TransitionProbabilityMatrix.h"

#include <cmath>
#include <iomanip>

using namespace RevBayesCore;

/** Construct rate matrix with n states, virtual population size, mutation rates, selection coefficients */
RateMatrix_BinaryMutationCoalescent::RateMatrix_BinaryMutationCoalescent(size_t n) : AbstractRateMatrix( n * (n+3) / 2.0 ),
    N( n ),
    matrix_size( (n*(n+3)/2.0) ),
    mu( 1.0 ),
    Ne( 1.0 )
{
    
    update();
}


/** Destructor */
RateMatrix_BinaryMutationCoalescent::~RateMatrix_BinaryMutationCoalescent(void)
{
    
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
    for ( size_t i=N; i>=1; --i )
    {
        for (size_t j=0; j<=i; ++j)
        {
            
            if ( i > j )
            {
                // we can have a mutation
                (*the_rate_matrix)[row_index][row_index+1] = (i-j)*mu;
            }
            
            if ( (i-j) > 1 )
            {
                // we can have a coalescent event in the ancestral state
                (*the_rate_matrix)[row_index][row_index+i+1] = (i-j)*(i-j-1)/(2.0*2.0*Ne);
            }
            
            if ( j > 1 )
            {
                // we can have a coalescent event in the ancestral state
                (*the_rate_matrix)[row_index][row_index+i] = j*(j-1)/(2.0*2.0*Ne);
            }
            
            ++row_index;
        }
    }
    
    // set the diagonal values
    setDiagonal();
    
    
    
    // rescale
    //rescaleToAverageRate( 1.0 );
    
    
}


/** Calculate the transition probabilities */
void RateMatrix_BinaryMutationCoalescent::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const
{
    
    //Now the instantaneous rate matrix has been filled up entirely.
    //We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
    double t = rate * (startAge - endAge);
    computeExponentialMatrixByRepeatedSquaring(t, P);
    
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


void RateMatrix_BinaryMutationCoalescent::update( void )
{
    
    if ( needs_update )
    {
        buildRateMatrix();
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
