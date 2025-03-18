#include "DiscretizeBetaFunction.h"
#include "DistributionBeta.h"
#include "RbMathFunctions.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

/**
 * Wrapper for dealing with DiscretizeBetaFunction used the TypedDagNode classes of type doubles
 * @param a double value
 * @param b double value
 * @param nc number of categories to use in the approximation
 * @param med a bool of whether to use the median values to represent each category. If false then the mean values are used
 */

RevBayesCore::DiscretizeBetaFunction::DiscretizeBetaFunction(const TypedDagNode<double> *a, const TypedDagNode<double> *b, const TypedDagNode<std::int64_t> *nc, bool med) : TypedFunction< RbVector<double> >( new RbVector<double>(nc->getValue(), 1.0) ),
alpha( a ),
beta( b ),
numCats(nc),
median(med)
{
    
    addParameter( alpha );
    addParameter( beta );
    addParameter( numCats );
    
}



RevBayesCore::DiscretizeBetaFunction* RevBayesCore::DiscretizeBetaFunction::clone( void ) const
{
    
    return new DiscretizeBetaFunction(*this);
}


void RevBayesCore::DiscretizeBetaFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == alpha)
    {
        alpha = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == beta)
    {
        beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}

void RevBayesCore::DiscretizeBetaFunction::update( void ) {
    
    double a = alpha->getValue();
    double b = beta->getValue();
    int nCats = numCats->getValue();
    double interval = 1.0 / (2.0 * nCats);
    double factor = (a /(a + b)) * nCats;
    
    if (median) {
        /* the median value for each category is used to represent all of the values
         in that category */
        for (int i=0; i<nCats; i++)
            (*value)[i] = RbStatistics::Beta::quantile(a, b, (i/double (nCats)) + interval);
        double t = 0.0;
        for (int i=0; i<nCats; i++)
            t += (*value)[i];
        for (int i=0; i<nCats; i++)
            (*value)[i] *= factor / t;
    }
    else
    {
        /* the mean value for each category is used to represent all of the values
         in that category */
        for (int i=0; i<nCats-1; i++) {
            /* calculate the points in the gamma distribution */
            (*value)[i] = RbStatistics::Beta::quantile(a, b, (i+1)/double (nCats));
            /* calculate the cumulative values */
            (*value)[i] = RbMath::incompleteBeta(a+1, b, (*value)[i]);
        }
        (*value)[nCats-1] = 1.0;
        /* calculate the relative values and rescale */
        for (int i=nCats-1; i>0; i--){
            (*value)[i] -= (*value)[i-1];
            (*value)[i] *= factor;
        }
        (*value)[0] *= factor;
    }
    
    //    *value = dist->quantile( p->getValue() );
}
