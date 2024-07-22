#include "RateMatrix_Chromosomes.h"

#include <cmath>

#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_Chromosomes::RateMatrix_Chromosomes(size_t n, AbstractRateMatrix::METHOD m) : AbstractRateMatrix( n+1, false, m ),
    matrixSize( n+1 )
{
    setGamma(1.0);
    setRho(1.0);
    setDelta(1.0);
    setEta(1.0);
    setGamma_l(1.0);
    setDelta_l(1.0);
    
    update();
}


/** Destructor */
RateMatrix_Chromosomes::~RateMatrix_Chromosomes(void) 
{
    
}

double RateMatrix_Chromosomes::averageRate(void) const 
{
    return 1.0;
}


void RateMatrix_Chromosomes::buildRateMatrix(void) 
{
    
    for (size_t i=0; i< matrixSize; i++)
    {
        for (size_t j=0; j< matrixSize; j++)
        {
			(*the_rate_matrix)[i][j] = 0.0;
//            if (j != 0 && i != 0)
			if (i != 0)
            {
				if (j == i+1)
                {
					(*the_rate_matrix)[i][j] += gamma * exp( gamma_l * (i-1) );
				}
                if (j == i-1)
                {
					(*the_rate_matrix)[i][j] += delta * exp( delta_l * (i-1) );
				}
                if (j == (2*i))
                {
					(*the_rate_matrix)[i][j] += rho;
                }
                if ( (i % 2 == 0) && (j == (size_t)(1.5*i)) )
                {
                    (*the_rate_matrix)[i][j] += eta;
                }
                if ( (i % 2 != 0) && ( (j == (size_t)(1.5*i - 0.5)) || (j == (size_t)(1.5*i + 0.5) ) ) )
                {
                    (*the_rate_matrix)[i][j] += eta;
                }
			}
        }
    }	
    // set the diagonal values
    setDiagonal();
    
    // rescale rates
    //rescaleToAverageRate( 1.0 );
}



RateMatrix_Chromosomes* RateMatrix_Chromosomes::clone( void ) const
{
    return new RateMatrix_Chromosomes( *this );
}


std::vector<double> RateMatrix_Chromosomes::getStationaryFrequencies( void ) const
{
    
    return stationary_freqs;
}

void RateMatrix_Chromosomes::setGamma( double g )
{
    
    gamma = g;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_Chromosomes::setRho( double r )
{
    
    rho = r;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_Chromosomes::setDelta( double d ) 
{
    
    delta = d;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_Chromosomes::setEta( double e ) 
{
    
    eta = e;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_Chromosomes::setGamma_l( double l ) 
{
    
    gamma_l = l;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_Chromosomes::setDelta_l( double d ) 
{
    
    delta_l = d;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_Chromosomes::updateInternalRateMatrix( void )
{
    
    if ( needs_update )
    {
        buildRateMatrix();
    }
}


