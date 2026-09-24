#include "PhyloBranchRatesBM.h"

#include <cstddef>
#include <cmath>

#include "DistributionMultivariateNormal.h"
#include "DistributionNormal.h"
#include "DistributionLognormal.h"
#include "Cloner.h"
#include "MatrixReal.h"
#include "RandomNumberFactory.h"
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
    const TopologyNode& root_node = tau->getValue().getRoot();

    std::vector<double> node_values = std::vector<double>(n_nodes, 0.0);
    if ( this->value->size() != (n_nodes-1) )
    {
        throw RbException("The dimension of the rates vector and the tree don't match.");
    }
    node_values[n_nodes-1] = root_state->getValue();
    double ln_prob = recursiveLnProb(root_node, node_values);
    
    
//    ln_prob += (n_nodes-1) * RbConstants::LN2;
    
    return ln_prob;
}


double PhyloBranchRatesBM::recursiveLnProb( const TopologyNode& node, std::vector<double> &parent_values )
{
    
    double ln_prob = 0.0;
    size_t index = node.getIndex();
    
    if ( node.isRoot() == false )
    {
        
//        // x ~ normal(x_up, sigma^2 * branchLength)
//        size_t parent_index = node.getParent().getIndex();
//        double parent_value = parent_values[parent_index];
//        double ln_parent_value = log( parent_value );
//        double branch_rate = (*this->value)[ index ];
//        // rate = (x+x_parent) / 2
//        double node_value = 2*branch_rate - parent_value;
//        if ( node_value < 0.0 )
//        {
//            return RbConstants::Double::neginf;
//        }
//        double ln_node_value = log(node_value);
//        double stand_dev = sigma->getValue() * sqrt(node.getBranchLength());
//        double ln_mean = ln_parent_value; // + drift->getValue() * node.getBranchLength();
//        ln_prob += RbStatistics::Normal::lnPdf(ln_mean, stand_dev, ln_node_value) - ln_node_value;
//        
//        parent_values[index] = node_value;
        
        size_t parent_index = node.getParent().getIndex();
        double parent_value = 0.0;
        double parent_branch_length = 0.0;
        if ( node.getParent().isRoot() == true )
        {
            parent_value = root_state->getValue();
        }
        else
        {
            parent_value = (*this->value)[parent_index];
            parent_branch_length = node.getParent().getBranchLength();
        }
        double ln_parent_value = log( parent_value );
        double node_value = (*this->value)[index];
        double ln_node_value = log(node_value);
        double ln_mean = ln_parent_value; // + drift->getValue() * node.getBranchLength();
        double stand_dev = sigma->getValue() * sqrt( (node.getBranchLength() + parent_branch_length)/2.0);
//        log rate at branch j ~ Normal (log rate at upstream branch, sigma^2 (length of branch j + length of upstream branch)/2)
//        ln_prob += RbStatistics::Normal::lnPdf(ln_mean, stand_dev, ln_node_value) - ln_node_value;
//        ln_prob += RbStatistics::Normal::lnPdf(ln_mean, stand_dev, ln_node_value);
        ln_prob += RbStatistics::Lognormal::lnPdf(ln_mean, stand_dev, node_value);

    }
    
//    if ( node.isTip() == false )
//    {
//        double node_value = 0.0;
//        if ( node.isRoot() == true )
//        {
//            node_value = root_state->getValue();
//        }
//        else
//        {
//            node_value = (*this->value)[index];
//        }
//        double ln_node_value = log(node_value);
//        
//        std::vector<double> ln_mean = std::vector<double>(2, ln_node_value);
//        
//        std::vector<double> ln_child = std::vector<double>(2, 0);
//        ln_child[0] = log( (*this->value)[ node.getChild(0).getIndex() ] );
//        ln_child[1] = log( (*this->value)[ node.getChild(1).getIndex() ] );
//        double parent_branch_length  = 0.0;
//        if ( node.isRoot() == false )
//        {
//            parent_branch_length = node.getBranchLength()/2.0
//        }
//        double child_0_branch_length = node.getChild(0).getBranchLength()/2.0;
//        double child_1_branch_length = node.getChild(1).getBranchLength()/2.0;
//        MatrixReal cov = MatrixReal(2);
//        cov[0][0] = parent_branch_length + child_0_branch_length;
//        cov[0][1] = parent_branch_length;
//        cov[1][0] = parent_branch_length;
//        cov[1][1] = parent_branch_length + child_1_branch_length;
//        ln_prob += RbStatistics::MultivariateNormal::lnPdfCovariance(ln_mean, cov, ln_child, sigma->getValue() * sigma->getValue());
//        ln_prob -= ln_child[0];
//        ln_prob -= ln_child[1];
//
//    }
    
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
        double ln_parent_value = log( parent_values[parent_index] );
        double stand_dev = sigma->getValue() * sqrt(node.getBranchLength());
        double ln_mean = ln_parent_value; // + drift->getValue() * node.getBranchLength();
        
        // simulate the new Val
        RandomNumberGenerator* rng = GLOBAL_RNG;
        parent_values[index] = exp(RbStatistics::Normal::rv( ln_mean, stand_dev, *rng));
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
