#include "DiscretizeGammaFunction.h"

#include "DistributionChisq.h"
#include "RbMathFunctions.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

/**
 * Wrapper for dealing with DiscretizeGammaFunction used the TypedDagNode classes of type doubles
 * @param s double value for the shape parameter
 * @param r double value for the rate parameter
 * @param nc number of categories to use in the approximation
 * @param med a bool of whether to use the median values to represent each category. If false then the mean values are used
 */


RevBayesCore::DiscretizeGammaFunction::DiscretizeGammaFunction(const TypedDagNode<double> *s, const TypedDagNode<double> *r, const TypedDagNode<std::int64_t> *nc, bool med) : TypedFunction< RbVector<double> >( new RbVector<double>(nc->getValue(), 1.0) ),
    shape( s ),
    rate( r ),
    numCats(nc),
    median(med)
{
    
    addParameter( shape );
    addParameter( rate );
    addParameter( numCats );
        
}



RevBayesCore::DiscretizeGammaFunction* RevBayesCore::DiscretizeGammaFunction::clone( void ) const
{
    
    return new DiscretizeGammaFunction(*this);
}


void RevBayesCore::DiscretizeGammaFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {
    
    if (oldP == shape) 
    {
        shape = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == rate)
    {
        rate = static_cast<const TypedDagNode<double>* >( newP );
    }
    
    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    
}

void RevBayesCore::DiscretizeGammaFunction::update( void )
{
    
    double a = shape->getValue();
    double b = rate->getValue();
    int nCats = (int)numCats->getValue();
    
    double factor = a / b * nCats;

    if (median) {
        /* the median value for each category is used to represent all of the values
        in that category */
        double interval = 1.0 / (2.0 * nCats);
        for (int i=0; i<nCats; i++) 
            (*value)[i] = RbStatistics::ChiSquare::quantile((i * 2.0 + 1.0) * interval, 2.0 * a) / (2.0 * b);
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
        /* calculate the points in the gamma distribution */
        for (int i=0; i<nCats-1; i++) 
            (*value)[i] = RbStatistics::ChiSquare::quantile((i + 1.0) / nCats, 2.0 * a) / (2.0 * b);
        /* calculate the cumulative values */
        for (int i=0; i<nCats-1; i++) 
            (*value)[i] = RbMath::incompleteGamma((*value)[i] * b, a + 1.0);
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
