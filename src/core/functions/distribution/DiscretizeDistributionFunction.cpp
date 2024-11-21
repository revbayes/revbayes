#include "DiscretizeDistributionFunction.h"

#include <vector>
#include <cstdint>

#include "Cloneable.h"
#include "ContinuousDistribution.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }


using namespace RevBayesCore;

/*
 * Default Constructor for the DiscretizeDistributionFunction
 *
 * @param d a type ContinuousDistribution for the distribution to be discretized
 * @param nc A type std::int64_t value for the number of categories
 *
 */


DiscretizeDistributionFunction::DiscretizeDistributionFunction(ContinuousDistribution *d, const TypedDagNode<std::int64_t> *nc) : TypedFunction< RbVector<double> >( new RbVector<double>(nc->getValue(), 1.0) ),
    dist( d ),
    num_cats(nc)
{
    
    addParameter( num_cats );

    const std::vector<const DagNode*>& params = dist->getParameters();
    for (std::vector<const DagNode* >::const_iterator it = params.begin(); it != params.end(); ++it)
    {
        addParameter( *it );
    }

}


DiscretizeDistributionFunction::DiscretizeDistributionFunction(const DiscretizeDistributionFunction &df) : TypedFunction< RbVector<double> >( df ),
    dist( df.dist->clone() ),
    num_cats(df.num_cats)
{
    
}



DiscretizeDistributionFunction::~DiscretizeDistributionFunction(void)
{
    
    delete dist;
    
}


DiscretizeDistributionFunction& DiscretizeDistributionFunction::operator=(const DiscretizeDistributionFunction &df)
{
    
    if ( this != &df )
    {
        TypedFunction< RbVector<double> >::operator=(df);
        
        delete dist;
        
        num_cats    = df.num_cats;
        dist        = df.dist->clone();

    }
    
    return *this;
}




DiscretizeDistributionFunction* DiscretizeDistributionFunction::clone( void ) const
{
    
    return new DiscretizeDistributionFunction(*this);
}


void DiscretizeDistributionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == num_cats)
    {
        num_cats = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
    else
    {
        dist->swapParameter(oldP, newP);
    }
    
}

void RevBayesCore::DiscretizeDistributionFunction::update( void )
{
    
    int n_cats = (int)num_cats->getValue();
    
    for (int i=0; i<n_cats; ++i)
    {
        double p = (i+0.5)/n_cats;
        (*value)[i] = dist->quantile( p );
    }
    
}
