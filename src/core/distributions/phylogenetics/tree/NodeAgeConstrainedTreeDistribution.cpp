#include "NodeAgeConstrainedTreeDistribution.h"

#include <cstddef>
#include <cmath>

#include "RbConstants.h"
#include "TopologyNode.h"
#include "RbBoolean.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbOrderedSet; }
namespace RevBayesCore { template <class variableType> class StochasticNode; }

using namespace RevBayesCore;


/**
 * Constructor.
 *
 * The constructor connects the parameters of the birth-death process (DAG structure)
 * and initializes the probability density by computing the combinatorial constant of the tree structure.
 *
 * \param[in]    o         Origin or time of the process.
 * \param[in]    cdt       The condition of the process (time/survival/nTaxa)
 * \param[in]    nTaxa     Number of taxa (used for initialization during simulation).
 * \param[in]    tn        Taxon names used during initialization.
 * \param[in]    c         Clade constraints.
 */
NodeAgeConstrainedTreeDistribution::NodeAgeConstrainedTreeDistribution(TypedDistribution<Tree> *base_dist, const std::vector< RbVector<Clade> > &c) : TypedDistribution<Tree>( new Tree() ),
    base_distribution( base_dist ),
    constraints( c )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    // we just set our pointer to the same pointer to the value of the base distribution
    value = &base_distribution->getValue();
    
}


NodeAgeConstrainedTreeDistribution::NodeAgeConstrainedTreeDistribution(const NodeAgeConstrainedTreeDistribution &d) : TypedDistribution<Tree>( d ),
    base_distribution( d.base_distribution->clone() ),
    constraints( d.constraints )
{
    // the copy constructor of the TypedDistribution creates a new copy of the value
    // however, here we want to hold exactly the same value as the base-distribution
    // thus, we delete the newly created value
    delete value;
    
    
    // we simply use the same pointer
    value = &base_distribution->getValue();
    
    
    // add the parameters of the base distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
}



NodeAgeConstrainedTreeDistribution::~NodeAgeConstrainedTreeDistribution()
{
    
    delete base_distribution;
    
    // DO NOT DELETE THE VALUE
    
}


NodeAgeConstrainedTreeDistribution* NodeAgeConstrainedTreeDistribution::clone( void ) const
{
    
    return new NodeAgeConstrainedTreeDistribution( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double NodeAgeConstrainedTreeDistribution::computeLnProbability( void )
{
    
    // first check if the current tree matches the clade constraints
    if ( matchesConstraints() == false )
    {
        return RbConstants::Double::neginf;
    }
    
    double ln_prob = base_distribution->computeLnProbability();
    
    return ln_prob;
}


/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void NodeAgeConstrainedTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    // delegate to the base distribution
    base_distribution->getAffected(affected, affecter);
}


/**
 * We check here if all the constraints are satisfied.
 * These are hard constraints, that is, the clades must be monophyletic.
 *
 * \return     True if the constraints are matched, false otherwise.
 */
bool NodeAgeConstrainedTreeDistribution::matchesConstraints( void )
{
    
    for (size_t j=0; j<constraints.size(); ++j)
    {

        if ( constraints[j].size() == 0 )
        {
            continue;
        }
        
        double constrain_age = value->getTmrca( constraints[j][0] );

        for (size_t i = 1; i < constraints[j].size() ; ++i)
        {
            double node_age      = value->getTmrca( constraints[j][i] );
        
            if ( fabs(constrain_age-node_age) > 1E-10 )
            {
                return false;
            }

        }

    }
    
    return true;
}



/**
 * Redraw the current value. We delegate this to the simulate method.
 */
void NodeAgeConstrainedTreeDistribution::redrawValue( void )
{
    
    base_distribution->redrawValue();
    // we just replace the current pointer with the pointer
    // to the new value of the base distribution
    value = &base_distribution->getValue();
    
    
}


/**
 * Set the DAG node.
 */
void NodeAgeConstrainedTreeDistribution::setStochasticNode( StochasticNode<Tree> *n )
{
    
    // delegate to base class first
    TypedDistribution<Tree>::setStochasticNode( n );
    
    if ( base_distribution != NULL )
    {
        base_distribution->setStochasticNode( n );
    }
    
}


/**
 * Set the current value.
 */
void NodeAgeConstrainedTreeDistribution::setValue(Tree *v, bool f )
{
    
    //  we set our value to the same value as the base distribution
    // but first we need to make sure that our base class doesn't delete the value
    value = NULL;
        
    // and the we can set it for both ourselves and the base distribution
    TypedDistribution<Tree>::setValue(v, f);
    base_distribution->setValue(v, f);
        
}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void NodeAgeConstrainedTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    base_distribution->swapParameter(oldP,newP);
    
}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void NodeAgeConstrainedTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    base_distribution->touch(affecter, touchAll);
    double a = base_distribution->getValue().getRoot().getAge();
    value->getRoot().setAge( a );
}

void NodeAgeConstrainedTreeDistribution::keepSpecialization(const DagNode *affecter)
{
    base_distribution->keep(affecter);
    double a = base_distribution->getValue().getRoot().getAge();
    value->getRoot().setAge( a );
}

void NodeAgeConstrainedTreeDistribution::restoreSpecialization(const DagNode *restorer)
{
    base_distribution->restore(restorer);
    double a = base_distribution->getValue().getRoot().getAge();
    value->getRoot().setAge( a );
    
}
