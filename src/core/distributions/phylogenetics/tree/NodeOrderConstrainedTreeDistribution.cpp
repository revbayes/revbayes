#include <cstddef>
#include <iosfwd>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "RelativeNodeAgeConstraints.h"
#include "NodeOrderConstrainedTreeDistribution.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "TreeUtilities.h"
#include "Cloneable.h"
#include "Tree.h"
#include "TypedDistribution.h"

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
NodeOrderConstrainedTreeDistribution::NodeOrderConstrainedTreeDistribution(TypedDistribution<Tree> *base_dist, const RelativeNodeAgeConstraints &c) : TypedDistribution<Tree>( new Tree() ),
    base_distribution( base_dist ),
    constraints( c ),
    constrained_nodes(),
    node_ages(),
    owns_tree( false )
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

    // if we own the tree, then we initialize our value with a true clone of the value of the base distribution
    if ( owns_tree == true )
    {
        value = base_distribution->getValue().clone();
    }
    else // otherwise we just set our pointer to the same pointer to the value of the base distribution
    {
        value = &base_distribution->getValue();
    }

    updateSetOfConstrainedNodes();
}


NodeOrderConstrainedTreeDistribution::NodeOrderConstrainedTreeDistribution(const NodeOrderConstrainedTreeDistribution &d) : TypedDistribution<Tree>( d ),
    base_distribution( d.base_distribution->clone() ),
    constraints( d.constraints ),
    constrained_nodes( d.constrained_nodes ),
    node_ages( d.node_ages ),
    owns_tree( d.owns_tree )
{
    // the copy constructor of the TypedDistribution creates a new copy of the value
    // however, here we want to hold exactly the same value as the base-distribution
    // thus, we delete the newly created value
    delete value;


    // and then set it to the value of the base distribution
    if ( owns_tree == true )
    {
        // if we own the tree, then we set it to a true copy
        value = base_distribution->getValue().clone();
    }
    else
    {
        // otherwise we simply use the same pointer
        value = &base_distribution->getValue();
    }


    // add the parameters of the base distribution
    const std::vector<const DagNode*>& pars = base_distribution->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }

}



NodeOrderConstrainedTreeDistribution::~NodeOrderConstrainedTreeDistribution()
{

    delete base_distribution;

    // DO NOT DELETE THE VALUE
    // the base distribution is the actual owner of the value!!!
    // we simply avoid the deletion of the value by setting its pointer to NULL
    // our base class, the TypedDistribution thinks that it owns the value and thus deletes it
    if ( owns_tree == false )
    {
        value = NULL;
    }

}


NodeOrderConstrainedTreeDistribution* NodeOrderConstrainedTreeDistribution::clone( void ) const
{

    return new NodeOrderConstrainedTreeDistribution( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double NodeOrderConstrainedTreeDistribution::computeLnProbability( void )
{

    // first check if the current tree matches the clade constraints
    if ( !matchesConstraints() )
    {
        return RbConstants::Double::neginf;
    }

    // since we and the base distribution own the same value,
    // we do not need to set the value of the base distribution
    if ( owns_tree == true )
    {
        base_distribution->setValue( value->clone() );
    }

    double lnProb = base_distribution->computeLnProbability();

    return lnProb;
}


/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void NodeOrderConstrainedTreeDistribution::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
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
bool NodeOrderConstrainedTreeDistribution::matchesConstraints( void )
{


    updateMapOfNodeAges();

    std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > > constra = constraints.getConstraints();
    /*
    for (size_t i = 0; i < constra.size() ; ++i)
    {
        constrained_nodes.insert(constra[i].first);
        constrained_nodes.insert(constra[i].second);
    }

    for (size_t i = 0; i < constra.size() ; ++i) {
        constrained_nodes.insert(constra[i].first);
        constrained_nodes.insert(constra[i].second);
    }*/

    for (size_t i = 0; i < constra.size() ; ++i)
    {
      //std::cout << "FIRST: "<<    node_ages.at(constra[i].first)<<" and SECOND: " << node_ages.at(constra[i].second) << std::endl;

      //std::cout << value->getPlainNewickRepresentation() << std::endl;

        if ( node_ages.at(constra[i].first) <  node_ages.at(constra[i].second) )
        {

          //std::cout << "NodeOrderConstrainedTreeDistribution::matchesConstraints: FALSE" <<std::endl;

            return false;
        }
    }
    //std::cout << "NodeOrderConstrainedTreeDistribution::matchesConstraints: TRUE; "<< constra.size() <<std::endl;

    return true;

}


void NodeOrderConstrainedTreeDistribution::updateSetOfConstrainedNodes()
{
    //std::vector <std::pair < std::pair<std::string, std::string>, std::pair<std::string, std::string> > >
    constra = constraints.getConstraints();
    for (size_t i = 0; i < constra.size() ; ++i)
    {
        constrained_nodes.insert(constra[i].first);
        constrained_nodes.insert(constra[i].second);
    }
    return;
}


//Here we compute node ages from the current tree.
void NodeOrderConstrainedTreeDistribution::updateMapOfNodeAges()
{

    node_ages.clear();
    for (std::set< std::pair < std::string, std::string > >::iterator elem=constrained_nodes.begin(); elem != constrained_nodes.end(); ++elem)
    {
        node_ages[(*elem)] = TreeUtilities::getAgeOfMRCA(*value, elem->first, elem->second);
    }

    
}



/**
 * Redraw the current value. We delegate this to the simulate method.
 */
void NodeOrderConstrainedTreeDistribution::redrawValue( void )
{

    base_distribution->redrawValue();
    // if we own the tree, then we need to free the memory before we create a new random variable
    if ( owns_tree == true )
    {
        delete value;
        value = base_distribution->getValue().clone();
    }
    else
    {
        // if we don't own the tree, then we just replace the current pointer with the pointer
        // to the new value of the base distribution
        value = &base_distribution->getValue();
    }


}


/**
 * Set the DAG node.
 */
void NodeOrderConstrainedTreeDistribution::setStochasticNode( StochasticNode<Tree> *n )
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
void NodeOrderConstrainedTreeDistribution::setValue(Tree *v, bool f )
{

    if ( owns_tree == true )
    {
        TypedDistribution<Tree>::setValue(v, f);

        // if we own the tree then we simply initialize the base distribution with a clone
        base_distribution->setValue(v->clone(), f);
    }
    else
    {
        // otherwise we set our value to the same value as the base distribution
        // but first we need to make sure that our base class doesn't delete the value
        value = NULL;

        // and the we can set it for both ourselves and the base distribution
        TypedDistribution<Tree>::setValue(v, f);
        base_distribution->setValue(v, f);
    }

    updateSetOfConstrainedNodes();


    //    if ( rootAge != NULL )
    //    {
    //        const StochasticNode<double> *stoch_root_age = dynamic_cast<const StochasticNode<double>* >(rootAge);
    //        if ( stoch_root_age != NULL )
    //        {
    //            const_cast<StochasticNode<double> *>(stoch_root_age)->setValue( new double( value->getRoot().getAge() ), f);
    //        }
    //        else
    //        {
    //            //            double factor = rootAge->getValue() / value->getRoot().getAge();
    //            //            TreeUtilities::rescaleTree( value, &value->getRoot(), factor);
    //
    //            value->getRoot().setAge( rootAge->getValue() );
    //        }
    //
    //    }

}


/**
 * Swap the parameters held by this distribution.
 *
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void NodeOrderConstrainedTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{

    base_distribution->swapParameter(oldP,newP);

}

/**
 * Touch the current value and reset some internal flags.
 * If the root age variable has been restored, then we need to change the root age of the tree too.
 */
void NodeOrderConstrainedTreeDistribution::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    base_distribution->touch(affecter, touchAll);
    double a = base_distribution->getValue().getRoot().getAge();
    value->getRoot().setAge( a );
}

void NodeOrderConstrainedTreeDistribution::keepSpecialization(const DagNode *affecter)
{
    base_distribution->keep(affecter);
    double a = base_distribution->getValue().getRoot().getAge();
    value->getRoot().setAge( a );
}

void NodeOrderConstrainedTreeDistribution::restoreSpecialization(const DagNode *restorer)
{
    base_distribution->restore(restorer);
    double a = base_distribution->getValue().getRoot().getAge();
    value->getRoot().setAge( a );

}
