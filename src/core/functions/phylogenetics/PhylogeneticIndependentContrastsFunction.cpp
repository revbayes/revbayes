#include "PhylogeneticIndependentContrastsFunction.h"

#include <cmath>

#include "ContinuousCharacterData.h"
#include "ContinuousTaxonData.h"
#include "RbVectorImpl.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PhylogeneticIndependentContrastsFunction::PhylogeneticIndependentContrastsFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t> *s, bool n ) : TypedFunction< RbVector<double> >( new RbVector<double>() ),
    tau( t ),
    data( d ),
    site( s ),
    normalized( n )
{
    // add the lambda parameter as a parent
    addParameter( tau );
    addParameter( data );
    addParameter( s );
    
    resetContrasts();
    
    update();
}


PhylogeneticIndependentContrastsFunction::~PhylogeneticIndependentContrastsFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PhylogeneticIndependentContrastsFunction* PhylogeneticIndependentContrastsFunction::clone( void ) const
{
    return new PhylogeneticIndependentContrastsFunction( *this );
}


void PhylogeneticIndependentContrastsFunction::update( void )
{
    
    RbVector<double> &v = *value;
    const Tree &tree = tau->getValue();
    
    size_t num_tips             = tree.getNumberOfTips();
    size_t num_internal_nodes   = num_tips-1;
    
    if ( v.size() != num_internal_nodes )
    {
        v.resize( num_internal_nodes );
    }
    
    const TopologyNode &root = tree.getRoot();
    recursiveComputeLnProbability(root, root.getIndex());
    
    for (size_t i=0; i<num_internal_nodes; ++i)
    {
        v[i] = node_means[ num_tips+i ];
        if ( normalized == true )
        {
            v[i] /= standard_deviations[ num_tips+i ];
        }
    }
    
}


void PhylogeneticIndependentContrastsFunction::recursiveComputeLnProbability( const TopologyNode &node, size_t node_index )
{
    
    // check for recomputation
    if ( node.isTip() == false )
    {
        
        double &mu_node  = this->contrasts[node_index];
        double &stdev  = this->standard_deviations[node_index];
        double &contrast  = this->node_means[node_index];
        
        // get the number of children
        size_t num_children = node.getNumberOfChildren();
        
        for (size_t j = 1; j < num_children; ++j)
        {
            
            size_t left_index = node_index;
            const TopologyNode *left = &node;
            if ( j == 1 )
            {
                left = &node.getChild(0);
                left_index = left->getIndex();
                recursiveComputeLnProbability( *left, left_index );
            }
            
            const TopologyNode &right = node.getChild(j);
            size_t right_index = right.getIndex();
            recursiveComputeLnProbability( right, right_index );
            
            
            // get the per node and site contrasts
            double &mu_left  = this->contrasts[left_index];
            double &mu_right = this->contrasts[right_index];
            
            // get the propagated uncertainties
            double delta_left  = this->contrast_uncertainty[left_index];
            double delta_right = this->contrast_uncertainty[right_index];
            
            // get the scaled branch lengths
            double v_left  = 0;
            if ( j == 1 )
            {
                v_left = left->getBranchLength();
            }
            double v_right = right.getBranchLength();
            
            // add the propagated uncertainty to the branch lengths
            double t_left  = v_left  + delta_left;
            double t_right = v_right + delta_right;
            
            // set delta_node = (t_l*t_r)/(t_l+t_r);
            this->contrast_uncertainty[node_index] = (t_left*t_right) / (t_left+t_right);
            
            stdev = sqrt(t_left+t_right);
            
            mu_node = (mu_left*t_right + mu_right*t_left) / (t_left+t_right);
                
            // compute the contrasts for this node
            contrast = mu_left - mu_right;
            
        } // end for-loop over all children
        
    } // end if we need to compute something for this node.
    
}


void PhylogeneticIndependentContrastsFunction::resetContrasts( void )
{
    size_t num_nodes = tau->getValue().getNumberOfNodes();
    
    // check if the vectors need to be resized
    contrasts                   = std::vector<double>(num_nodes, 0);
    contrast_uncertainty        = std::vector<double>(num_nodes, 0);
    standard_deviations         = std::vector<double>(num_nodes, 0);
    node_means                  = std::vector<double>(num_nodes, 0);
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    size_t site_index = site->getValue()-1;
    
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        if ( (*it)->isTip() )
        {
            const ContinuousTaxonData& taxon = data->getValue().getTaxonData( (*it)->getName() );
            double c = taxon.getCharacter(site_index);
            contrasts[(*it)->getIndex()] = c;
            contrast_uncertainty[(*it)->getIndex()] = 0;
        }
    }
    
}

void PhylogeneticIndependentContrastsFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == tau)
    {
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
    }
    else if (oldP == data)
    {
        data = static_cast<const TypedDagNode<ContinuousCharacterData>* >( newP );
    }
    else if (oldP == site)
    {
        site = static_cast<const TypedDagNode<std::int64_t>* >( newP );
    }
}



