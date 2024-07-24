#include "RateMatrix_GeneGainLoss.h"

#include <cmath>

#include "MatrixReal.h"
#include "TransitionProbabilityMatrix.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

using namespace RevBayesCore;

/** Construct rate matrix with n states */
RateMatrix_GeneGainLoss::RateMatrix_GeneGainLoss(size_t n, AbstractRateMatrix::METHOD m) : AbstractRateMatrix( n+1, false, m )
{
    birth = 2.0;
    death = 1.0;
    
    update();
}


/** Destructor */
RateMatrix_GeneGainLoss::~RateMatrix_GeneGainLoss(void)
{
    
}

double RateMatrix_GeneGainLoss::averageRate(void) const
{
    return 1.0;
}


void RateMatrix_GeneGainLoss::buildRateMatrix(void)
{
    
    for (size_t i=0; i< num_states; i++)
    {
        for (size_t j=0; j< num_states; j++)
        {
            (*the_rate_matrix)[i][j] = 0.0;
//            if (j != 0 && i != 0)
            if (i != 0)
            {
                if (j == i+1)
                {
                    (*the_rate_matrix)[i][j] += birth * i;
                }
                if (j == i-1)
                {
                    (*the_rate_matrix)[i][j] += death * i;
                }
//                if (j == (2*i))
//                {
//                    (*the_rate_matrix)[i][j] += rho;
//                }
//                if ( (i % 2 == 0) && (j == (size_t)(1.5*i)) )
//                {
//                    (*the_rate_matrix)[i][j] += eta;
//                }
//                if ( (i % 2 != 0) && ( (j == (size_t)(1.5*i - 0.5)) || (j == (size_t)(1.5*i + 0.5) ) ) )
//                {
//                    (*the_rate_matrix)[i][j] += eta;
//                }
            }
        }
    }
    // set the diagonal values
    setDiagonal();
    
    // rescale rates
    //rescaleToAverageRate( 1.0 );
}



RateMatrix_GeneGainLoss* RateMatrix_GeneGainLoss::clone( void ) const
{
    return new RateMatrix_GeneGainLoss( *this );
}


void RateMatrix_GeneGainLoss::setBirth( double b )
{
    
    birth = b;
    
    // set flags
    needs_update = true;
    
}

void RateMatrix_GeneGainLoss::setDeath( double d )
{
    
    death = d;
    
    // set flags
    needs_update = true;
    
}


void RateMatrix_GeneGainLoss::updateInternalRateMatrix( void )
{
    
    if ( needs_update )
    {
        buildRateMatrix();
    }
}


