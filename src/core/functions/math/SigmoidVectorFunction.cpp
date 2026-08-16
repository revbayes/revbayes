#include "SigmoidVectorFunction.h"
#include "DagNodeTypeUtilities.h"

#include <cmath>

#include "Cloneable.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * SigmoidVectorFunction of a TypedDagNode containing a value of type double
 *
 * @param a the value of type double
 *
 */
SigmoidVectorFunction::SigmoidVectorFunction(const TypedDagNode<RbVector<double> > *a, const TypedDagNode<double> *min_, const TypedDagNode<double> *max_, const TypedDagNode<double> *middle_, const TypedDagNode<double> *slope_) : TypedFunction<RbVector<double> >( new RbVector<double>(a->getValue().size(), 0.0) ),
    x( a ),
	min(min_),
	max(max_),
	middle(middle_),
	slope(slope_)
{
    addParameter( a );
    addParameter( min_ );
    addParameter( max_ );
    addParameter( middle_ );
    addParameter( slope_ );
}


SigmoidVectorFunction* SigmoidVectorFunction::clone( void ) const
{
    return new SigmoidVectorFunction(*this);
}


void SigmoidVectorFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == x)
    {
        replaceNodeReference(x, newP);
    }
    if (oldP == min)
    {
        replaceNodeReference(min, newP);
    }
    if (oldP == max)
    {
        replaceNodeReference(max, newP);
    }
    if (oldP == middle)
    {
        replaceNodeReference(middle, newP);
    }
    if (oldP == slope)
    {
        replaceNodeReference(slope, newP);
    }

}

void SigmoidVectorFunction::update( void )
{
    // recompute and set the internal value
	double y0 = min->getValue();
	double y1 = max->getValue();
	double m  = middle->getValue();
	double r  = slope->getValue();

    size_t n = x->getValue().size();
    for (size_t i = 0; i < n; ++i)
    {
    	const double &this_x = x->getValue()[i];
        (*value)[i] = y0 + (y1 - y0) / (1.0 + std::exp(r * (m - this_x)));
    }

}


