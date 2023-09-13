#include "EarlyBurstRatesFunction.h"

#include <cstddef>
#include <cmath>

#include "Cloner.h"
#include "RbVectorImpl.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

EarlyBurstRatesFunction::EarlyBurstRatesFunction(const TypedDagNode<Tree> *t, const TypedDagNode<double> *s, const TypedDagNode<double> *l) : TypedFunction<RbVector<double> >( new RbVector<double>(t->getValue().getNumberOfNodes()-1,0.0) ),
    tau( t ),
    sigma( s ),
    lambda( l )
{
    // add the lambda parameter as a parent
    addParameter( tau );
    addParameter( sigma );
    addParameter( lambda );
    
    update();
}


EarlyBurstRatesFunction::~EarlyBurstRatesFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



EarlyBurstRatesFunction* EarlyBurstRatesFunction::clone( void ) const
{
    return new EarlyBurstRatesFunction( *this );
}


void EarlyBurstRatesFunction::update( void )
{
    const Tree &tree = tau->getValue();
    double root_age = tree.getRoot().getAge();
    double s = sigma->getValue();
    double l = lambda->getValue();

    // tip nodes have pre-set ages
    // NOTE: should be able to rescale with v* tree.getAge(i) under new setup...
    for (size_t i = 0; i < tree.getNumberOfNodes(); ++i)
    {
        const TopologyNode &node = tree.getNode( i );
        if ( node.isRoot() == false )
        {
            double my_age = node.getAge();
            double parent_age = node.getParent().getAge();
            if ( l == 0 )
            {
//                (*value)[i]  = s * ( exp(l*parent_age) * parent_age - exp(l*my_age) * my_age ) / (parent_age-my_age);
                (*value)[i]  = s;
            }
            else
            {
//                (*value)[i]  = s*s * ( exp(l*parent_age) * parent_age - exp(l*my_age) * my_age )  / l / (parent_age-my_age);
//                (*value)[i]  = s * ( exp(l*parent_age) * parent_age - exp(l*my_age) * my_age ) / (parent_age-my_age);
                (*value)[i]  = s/l * ( exp( (parent_age-root_age)*l ) - exp( (my_age-root_age)*l ) ) / (parent_age-my_age);
            }
            
        }
        
    }
    
}



void EarlyBurstRatesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    if (oldP == sigma)
    {
        sigma = static_cast<const TypedDagNode<double>* >( newP );
    }
    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
}


