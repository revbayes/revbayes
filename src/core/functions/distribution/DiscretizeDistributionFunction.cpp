#include "DiscretizeDistributionFunction.h"

#include <vector>

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
 * @param nc A type long value for the number of categories
 *
 */


DiscretizeDistributionFunction::DiscretizeDistributionFunction(ContinuousDistribution *d, const TypedDagNode<long> *nc) : TypedFunction< RbVector<double> >( new RbVector<double>(nc->getValue(), 1.0) ),
    dist( d ),
    numCats(nc)
{
    
    addParameter( numCats );

    const std::vector<const DagNode*>& params = dist->getParameters();
    for (std::vector<const DagNode* >::const_iterator it = params.begin(); it != params.end(); ++it)
    {
        addParameter( *it );
    }

}


DiscretizeDistributionFunction::DiscretizeDistributionFunction(const DiscretizeDistributionFunction &df) : TypedFunction< RbVector<double> >( df ),
    dist( df.dist->clone() ),
    numCats(df.numCats)
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
        
        numCats = df.numCats;
        dist = df.dist->clone();

    }
    
    return *this;
}




DiscretizeDistributionFunction* DiscretizeDistributionFunction::clone( void ) const
{
    
    return new DiscretizeDistributionFunction(*this);
}


void DiscretizeDistributionFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == numCats)
    {
        numCats = static_cast<const TypedDagNode<long>* >( newP );
    }
    else
    {
        dist->swapParameter(oldP, newP);
    }
    
}

void RevBayesCore::DiscretizeDistributionFunction::update( void )
{
    
    int nCats = (int)numCats->getValue();
    
    for (int i=0; i<nCats; ++i)
    {
        double p = (i+0.5)/nCats;
        (*value)[i] = dist->quantile( p );
    }
    
}
