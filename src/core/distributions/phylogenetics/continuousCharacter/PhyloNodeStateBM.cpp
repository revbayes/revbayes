#include "PhyloNodeStateBM.h"

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
PhyloNodeStateBM::PhyloNodeStateBM(const TypedDagNode< Tree > *t, const TypedDagNode< double >* r, const TypedDagNode< double >* s, const TypedDagNode< double >* d): TypedDistribution< RbVector<double> >(new RbVector<double>(t->getValue().getNumberOfNodes()-1,0.0)),
    tau( t ),
    root_state( r ),
    sigma( s ),
    drift( d )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    addParameter( tau );
    addParameter( root_state );
    addParameter( sigma );
    addParameter( drift );
    
    simulate();
}



PhyloNodeStateBM* PhyloNodeStateBM::clone(void) const
{
    return new PhyloNodeStateBM( *this );
}


double PhyloNodeStateBM::computeLnProbability(void)
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


void PhyloNodeStateBM::getAffected(RbOrderedSet<DagNode *> &affected, const DagNode *affecter)
{
    
    if ( affecter == root_state && dag_node != NULL )
    {
        dag_node->initiateGetAffectedNodes( affected );
    }
    
}


void PhyloNodeStateBM::keepSpecialization(const DagNode *affecter)
{
    
    if ( affecter == root_state && dag_node != NULL)
    {
        dag_node->keepAffected();
    }
    
}


void PhyloNodeStateBM::restoreSpecialization(const DagNode *affecter)
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


double PhyloNodeStateBM::recursiveLnProb( const TopologyNode& node )
{
    
    double ln_prob = 0.0;
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        size_t parent_index = node.getParent().getIndex();
        double parent_value = (*this->value)[parent_index];
        double node_value = (*this->value)[index];
        double stand_dev = sigma->getValue() * sqrt(node.getBranchLength());
        ln_prob += RbStatistics::Normal::lnPdf(node_value, stand_dev, parent_value);
        
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    
    for (size_t i = 0; i < num_children; ++i)
    {
        ln_prob += recursiveLnProb(node.getChild(i));
    }
    
    return ln_prob;
    
}

void PhyloNodeStateBM::redrawValue(void)
{
    simulate();
}


void PhyloNodeStateBM::simulate()
{
    
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    (*this->value) = RbVector<double>(n_nodes, 0.0);
    (*this->value)[n_nodes-1] = root_state->getValue();
    recursiveSimulate(tau->getValue().getRoot());
}


void PhyloNodeStateBM::recursiveSimulate(const TopologyNode& node)
{
    
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        
        size_t parent_index = node.getParent().getIndex();
        double parent_value = (*this->value)[parent_index];
        double stand_dev = sigma->getValue() * sqrt(node.getBranchLength());
        
        // simulate the new Val
        RandomNumberGenerator* rng = GLOBAL_RNG;
        (*this->value)[index] = RbStatistics::Normal::rv( parent_value, stand_dev, *rng);
        
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    for (size_t i = 0; i < num_children; ++i)
    {
        recursiveSimulate(node.getChild(i));
    }
    
}

/** Swap a parameter of the distribution */
void PhyloNodeStateBM::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
    
    if ( oldP == drift )
    {
        drift = static_cast< const TypedDagNode< double > * >( newP );
    }
    
}


void PhyloNodeStateBM::touchSpecialization(const DagNode *affecter, bool touchAll)
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
