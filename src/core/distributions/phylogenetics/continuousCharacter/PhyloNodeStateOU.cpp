#include "PhyloNodeStateOU.h"

#include <cstddef>
#include <cmath>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "Cloner.h"
#include "RbVectorImpl.h"
#include "StochasticNode.h"
#include "Tree.h"
#include "TypedDagNode.h"



using namespace RevBayesCore;



// constructor(s)
PhyloNodeStateOU::PhyloNodeStateOU(const TypedDagNode< Tree > *tr, const TypedDagNode< double >* r, const TypedDagNode< double >* s, const TypedDagNode< double >* a, const TypedDagNode< double >* th): TypedDistribution< RbVector<double> >(new RbVector<double>(tr->getValue().getNumberOfNodes()-1,0.0)),
    tau( tr ),
    root_state( r ),
    sigma( s ),
    alpha( a ),
    theta( th )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( tau );
    addParameter( root_state );
    addParameter( sigma );
    addParameter( alpha );
    addParameter( theta );
    
    simulate();
}



PhyloNodeStateOU* PhyloNodeStateOU::clone(void) const
{
    return new PhyloNodeStateOU( *this );
}


double PhyloNodeStateOU::computeLnProbability(void)
{
    size_t n_nodes = tau->getValue().getNumberOfNodes();

    if ( this->value->size() != n_nodes )
    {
        throw RbException("The dimension of the rates vector and the tree don't match.");
    }
    
    double ln_prob = 0.0;
    if ( root_state->getValue() != (*this->value)[n_nodes-1] )
    {
        ln_prob = RbConstants::Double::neginf;
    }
    else
    {
        ln_prob = recursiveLnProb(tau->getValue().getRoot());
    }
    
    return ln_prob;
}


void PhyloNodeStateOU::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == root_state && dag_node != NULL )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}


void PhyloNodeStateOU::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_state && dag_node != NULL)
    {
        dag_node->keepAffected();
    }
    
}


void PhyloNodeStateOU::restoreSpecialization(const DagNode *affecter)
{
    if ( affecter == root_state )
    {
        (*this->value)[ tau->getValue().getRoot().getIndex() ] = root_state->getValue();
        
        if ( dag_node != NULL )
        {
            dag_node->restoreAffected();
        }
    }
    
}


double PhyloNodeStateOU::recursiveLnProb( const TopologyNode& node )
{
    
    double ln_prob = 0.0;
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        size_t parent_index = node.getParent().getIndex();
        double parent_value = (*this->value)[parent_index];
        double node_value = (*this->value)[index];
        
        double a = alpha->getValue();
        double t = node.getBranchLength();
        double e = exp(-a*t);
        double e2 = exp(-2.0*a*t);
        double m = e * parent_value + (1.0-e)*theta->getValue();
        double s = sigma->getValue() * sqrt( (1-e2) / 2 / a );

        ln_prob += RbStatistics::Normal::lnPdf(node_value, s, m);
        
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    
    for (size_t i = 0; i < num_children; ++i)
    {
        ln_prob += recursiveLnProb(node.getChild(i));
    }
    
    return ln_prob;
    
}

void PhyloNodeStateOU::redrawValue(void)
{
    simulate();
}


void PhyloNodeStateOU::simulate()
{
    
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    (*this->value) = RbVector<double>(n_nodes, 0.0);
    (*this->value)[n_nodes-1] = root_state->getValue();
    recursiveSimulate(tau->getValue().getRoot());
}


void PhyloNodeStateOU::recursiveSimulate(const TopologyNode& node)
{
    
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        
        size_t parent_index = node.getParent().getIndex();
        double parent_value = (*this->value)[parent_index];
        
        double a = alpha->getValue();
        double t = node.getBranchLength();
        double e = exp(-a*t);
        double e2 = exp(-2.0*a*t);
        double m = e * parent_value + (1.0-e)*theta->getValue();
        double s = sigma->getValue() * sqrt( (1-e2) / 2 / a );
        
        (*this->value)[index] = RbStatistics::Normal::rv( m, s, *GLOBAL_RNG);
        
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    for (size_t i = 0; i < num_children; ++i)
    {
        recursiveSimulate(node.getChild(i));
    }
    
}

/** Swap a parameter of the distribution */
void PhyloNodeStateOU::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if ( oldP == tau )
    {
        tau = static_cast< const TypedDagNode<Tree> * >( newP );
    }
    
    if ( oldP == root_state )
    {
        root_state = static_cast< const TypedDagNode<double> * >( newP );
    }
    
    if ( oldP == sigma )
    {
        sigma = static_cast< const TypedDagNode<double> * >( newP );
    }
    
    if ( oldP == alpha )
    {
        alpha = static_cast< const TypedDagNode< double > * >( newP );
    }
    
    if ( oldP == theta )
    {
        theta = static_cast< const TypedDagNode< double > * >( newP );
    }
    
}


void PhyloNodeStateOU::touchSpecialization(const DagNode *affecter, bool touchAll)
{
    
    if ( affecter == root_state )
    {
        (*this->value)[ tau->getValue().getRoot().getIndex() ] = root_state->getValue();

        
        if ( dag_node != NULL )
        {
            dag_node->touchAffected();
        }
        
    }
    
}
