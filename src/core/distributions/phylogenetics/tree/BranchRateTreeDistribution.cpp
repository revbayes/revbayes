#include "BranchRateTreeDistribution.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TreeChangeEventMessage.h"

using namespace RevBayesCore;

BranchRateTreeDistribution::BranchRateTreeDistribution(const TypedDagNode<Tree>* tt, TypedDistribution<double>* brp) : TypedDistribution<Tree>( new Tree() ),
    branch_rate_prior( brp ),
    time_tree( tt )
{
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    
    simulateTree();
    
}


BranchRateTreeDistribution::BranchRateTreeDistribution(const BranchRateTreeDistribution &d) : TypedDistribution<Tree>( d ),
    branch_rate_prior( d.branch_rate_prior->clone() ),
    time_tree( d.time_tree )
{
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    
}



BranchRateTreeDistribution::~BranchRateTreeDistribution()
{
    
    delete branch_rate_prior;
    // the tree will be deleted automatically by the base class
    
}


BranchRateTreeDistribution& BranchRateTreeDistribution::operator=(const BranchRateTreeDistribution &d)
{
    
    if ( this != &d )
    {
        TypedDistribution<Tree>::operator=( d );
        
        // remove the old branch-length-prior parameters
        const std::vector<const DagNode*>& old_pars = branch_rate_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = old_pars.begin(); it != old_pars.end(); ++it)
        {
            this->removeParameter( *it );
        }
        delete branch_rate_prior;
        
        branch_rate_prior           = d.branch_rate_prior->clone();
        time_tree                   = d.time_tree;

        // add the parameters of the distribution
        const std::vector<const DagNode*>& pars = branch_rate_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    
    }

    return *this;
}


BranchRateTreeDistribution* BranchRateTreeDistribution::clone( void ) const
{
    
    return new BranchRateTreeDistribution( *this );
}


double BranchRateTreeDistribution::computeLnProbability( void )
{
    
    double ln_prob = 0.0;
    
    // make the time tree unrooted
    // see code in UltrametricTreeDistribution
    
    // compare if it matches our topology
    
    // compute the branch rates as r = bl / t
    
    // compute the probability of each rate
    
    return ln_prob;
}


void BranchRateTreeDistribution::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
    
    if ( m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY )
    {
        
//        dirty_topology = true;
    }
    
    
}

void BranchRateTreeDistribution::redrawValue( void )
{
    simulateTree();
}


void BranchRateTreeDistribution::setValue(RevBayesCore::Tree *v, bool force)
{
    
    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );
    
    // Sebastian: check if anything special needs to be done
}


void BranchRateTreeDistribution::simulateTree( void )
{
    
    // make our own copy of this time tree as an unrooted tree
    
    // draw new branch rates
    
    
}


/** Swap a parameter of the distribution */
void BranchRateTreeDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    if ( branch_rate_prior != NULL )
    {
        branch_rate_prior->swapParameter(oldP,newP);
    }
}

