#include "PhyloBranchRatesBM.h"

#include <cstddef>
#include <cmath>

#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "Cloner.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"



using namespace RevBayesCore;



// constructor(s)
PhyloBranchRatesBM::PhyloBranchRatesBM(const TypedDagNode< Tree > *t, const TypedDagNode< double >* r, const TypedDagNode< double >* s, const TypedDagNode< double >* d): TypedDistribution< RbVector<double> >(new RbVector<double>(t->getValue().getNumberOfNodes()-1,0.0)),
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



PhyloBranchRatesBM* PhyloBranchRatesBM::clone(void) const
{
    return new PhyloBranchRatesBM( *this );
}


double PhyloBranchRatesBM::computeLnProbability(void)
{
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    std::vector<double> node_values = std::vector<double>(n_nodes, 0.0);
    if ( this->value->size() != (n_nodes-1) )
    {
        throw RbException("The dimension of the rates vector and the tree don't match.");
    }
    node_values[n_nodes-1] = log(root_state->getValue());
    double ln_prob = recursiveLnProb(tau->getValue().getRoot(), node_values);
    
    return ln_prob;
}


double PhyloBranchRatesBM::recursiveLnProb( const TopologyNode& node, std::vector<double> &parent_values )
{
    
    double ln_prob = 0.0;
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        size_t parent_index = node.getParent().getIndex();
        double parent_value = parent_values[parent_index];
        double branch_rate = (*this->value)[ index ];
        double ln_node_value = 2*branch_rate - exp(parent_value);
        if ( ln_node_value < 0.0 )
        {
            return RbConstants::Double::neginf;
        }
        double node_value = log(ln_node_value);
        double stand_dev = sigma->getValue() * sqrt(node.getBranchLength());
        double mean = parent_value + drift->getValue() * node.getBranchLength();
        ln_prob += RbStatistics::Normal::lnPdf(node_value, stand_dev, mean) - std::log(ln_node_value);
        
        parent_values[index] = node_value;
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    
    for (size_t i = 0; i < num_children; ++i)
    {
        ln_prob += recursiveLnProb(node.getChild(i), parent_values);
    }
    
    return ln_prob;
    
}

void PhyloBranchRatesBM::redrawValue(void)
{
    simulate();
}


void PhyloBranchRatesBM::simulate()
{
    
    size_t n_nodes = tau->getValue().getNumberOfNodes();
    std::vector<double> node_values = std::vector<double>(n_nodes, 0.0);
    node_values[n_nodes-1] = root_state->getValue();
    recursiveSimulate(tau->getValue().getRoot(), node_values);
}


void PhyloBranchRatesBM::recursiveSimulate(const TopologyNode& node, std::vector<double> &parent_values)
{
    
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
        // x ~ normal(x_up, sigma^2 * branchLength)
        
        size_t parent_index = node.getParent().getIndex();
        double stand_dev = sigma->getValue() * sqrt(node.getBranchLength());
        double mean = log(parent_values[parent_index]) + drift->getValue() * node.getBranchLength();
        
        // simulate the new Val
        RandomNumberGenerator* rng = GLOBAL_RNG;
        parent_values[index] = exp(RbStatistics::Normal::rv( mean, stand_dev, *rng));
        (*this->value)[index] = (parent_values[parent_index] + parent_values[index]) / 2.0;
        
    }
    
    // propagate forward
    size_t num_children = node.getNumberOfChildren();
    for (size_t i = 0; i < num_children; ++i)
    {
        recursiveSimulate(node.getChild(i), parent_values);
    }
    
}

/** Swap a parameter of the distribution */
void PhyloBranchRatesBM::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
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
