#include "RateMatrix_ChromosomesPloidy.h"

#include <cmath>

#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix where n = max_chromo */
RateMatrix_ChromosomesPloidy::RateMatrix_ChromosomesPloidy(size_t n) : AbstractRateMatrix( 2 * (n + 1), false, AbstractRateMatrix::EIGEN ),
    matrix_size( 2 * (n + 1) ),
    max_chromo( n )
{
    setGamma_d(1.0);
    setRho_d(1.0);
    setDelta_d(1.0);
    setEta_d(1.0);
    setGamma_dl(1.0);
    setDelta_dl(1.0);
    setGamma_p(1.0);
    setRho_p(1.0);
    setDelta_p(1.0);
    setEta_p(1.0);
    setGamma_pl(1.0);
    setDelta_pl(1.0);
    setBeta(1.0);
    
    update();
}


/** Destructor */
RateMatrix_ChromosomesPloidy::~RateMatrix_ChromosomesPloidy(void) 
{
    
}

double RateMatrix_ChromosomesPloidy::averageRate(void) const 
{
    return 1.0;
}


void RateMatrix_ChromosomesPloidy::buildRateMatrix(void) 
{
    
    for (size_t i = 0; i < matrix_size; i++)
    {
        for (size_t j = 0; j < matrix_size; j++)
        {
			(*the_rate_matrix)[i][j] = 0.0;
			if ( j != 0 && i != 0 )
            {
                // diploid chromosome evolution
                if ( i <= max_chromo && j <= max_chromo )
                {
                    if (j == i + 1)
                    {
                        (*the_rate_matrix)[i][j] += gamma_d * exp( gamma_dl * (i - 1) );
                    }
                    if (j == i - 1)
                    {
                        (*the_rate_matrix)[i][j] += delta_d * exp( delta_dl * (i - 1) );
                    }
                }
                // diploid to polyploid
                else if ( i <= max_chromo && j > (max_chromo + 1) )
                {
                    if (j == ((2 * i) + (max_chromo + 1)))
                    {
                        (*the_rate_matrix)[i][j] += rho_d;
                    }
                    if ( (i % 2 == 0) && (j == (size_t)((1.5 * i) + (max_chromo + 1))) )
                    {
                        (*the_rate_matrix)[i][j] += eta_d;
                    }
                    if ( (i % 2 != 0) && ( (j == (size_t)((1.5 * i) - (0.5 + max_chromo + 1))) || (j == (size_t)((1.5 * i) + (0.5 + max_chromo + 1)) ) ) )
                    {
                        (*the_rate_matrix)[i][j] += eta_d;
                    }
                }
                // polyploid chromosome evolution
                else if ( i > (max_chromo + 1) && j > (max_chromo + 1) )
                {
                    if ( j == i + 1 )
                    {
                        (*the_rate_matrix)[i][j] += gamma_p * exp( gamma_pl * (i - 2 - max_chromo) );
                    }
                    if ( j == i - 1 )
                    {
                        (*the_rate_matrix)[i][j] += delta_p * exp( delta_pl * (i - 2 - max_chromo) );
                    }
                    if ( (j - (max_chromo + 1)) == (2 * (i - (max_chromo + 1))) )
                    {
                        (*the_rate_matrix)[i][j] += rho_p;
                    }
                    if ( ((i - (max_chromo + 1)) % 2 == 0) && ( (j - (max_chromo + 1)) == (size_t)(1.5 * (i - (max_chromo + 1))) ) )
                    {
                        (*the_rate_matrix)[i][j] += eta_p;
                    }
                    if ( ((i - (max_chromo + 1)) % 2 != 0) && ( ((j - (max_chromo + 1)) == (size_t)(1.5 * (i - (max_chromo + 1)) - 0.5)) || ( (j - (max_chromo + 1)) == (size_t)(1.5 * (i - (max_chromo + 1)) + 0.5) ) ) )
                    {
                        (*the_rate_matrix)[i][j] += eta_p;
                    }
                }
                // rediploidization
                else if (i > (max_chromo + 1) && j <= max_chromo)
                {
                    if ( (i - (max_chromo + 1)) == j )
                    {
                        (*the_rate_matrix)[i][j] += beta;
                    }
                }
			}
        }
    }	
    // set the diagonal values
    setDiagonal();
    
    // rescale rates
    //rescaleToAverageRate( 1.0 );
}




/** Calculate the transition probabilities */
void RateMatrix_ChromosomesPloidy::calculateTransitionProbabilities(double startAge, double endAge, double rate, TransitionProbabilityMatrix& P) const {
    
    // We use repeated squaring to quickly obtain exponentials, as in Poujol and Lartillot, Bioinformatics 2014.
	// Mayrose et al. 2010 also used this method for chromosome evolution (named the squaring and scaling method in Moler and Van Loan 2003).
    double t = rate * (startAge - endAge);
    exponentiateMatrixByScalingAndSquaring(t, P);
    
}

RateMatrix_ChromosomesPloidy* RateMatrix_ChromosomesPloidy::clone( void ) const
{
    return new RateMatrix_ChromosomesPloidy( *this );
}


std::vector<double> RateMatrix_ChromosomesPloidy::getStationaryFrequencies( void ) const
{
    
    return stationary_freqs;
}

void RateMatrix_ChromosomesPloidy::setGamma_d( double g )
{
    
    gamma_d = g;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setRho_d( double r )
{
    
    rho_d = r;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setDelta_d( double d ) 
{
    
    delta_d = d;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setEta_d( double e ) 
{
    
    eta_d = e;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setGamma_dl( double l ) 
{
    
    gamma_dl = l;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setDelta_dl( double d ) 
{
    
    delta_dl = d;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setGamma_p( double g )
{
    
    gamma_p = g;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setRho_p( double r )
{
    
    rho_p = r;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setDelta_p( double d ) 
{
    
    delta_p = d;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setEta_p( double e ) 
{
    
    eta_p = e;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setGamma_pl( double l ) 
{
    
    gamma_pl = l;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setDelta_pl( double d ) 
{
    
    delta_pl = d;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::setBeta( double b ) 
{
    
    beta = b;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_ChromosomesPloidy::update( void ) {
    
    if ( needs_update )
    {
        buildRateMatrix();
        // clean flags
        needs_update = false;
    }
}


