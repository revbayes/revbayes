#include "UniformTopologyBranchLengthDistribution.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>

#include "Clade.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "RbBitSet.h"
#include "TimeInterval.h"
#include "TreeChangeEventMessage.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

UniformTopologyBranchLengthDistribution::UniformTopologyBranchLengthDistribution(const std::vector<Taxon> &tn, const Clade &og, TypedDistribution<double>* blp, bool rt) : TypedDistribution<Tree>( new Tree() ),
    branch_length_prior( blp ),
    num_taxa( tn.size() ),
    taxa( tn ),
    log_tree_topology_prob( RbConstants::Double::nan ),
    outgroup( og ),
    outgroup_provided( og.size() > 0 ),
    rooted( rt ),
    dirty_topology( true )
{
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_length_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    double branchLnFact = 0.0;
    double nodeLnFact = 0.0;
    for (size_t i = 2; i <= 2*num_taxa - 3 - 2*(!rooted); i++)
    {
        branchLnFact += std::log(i);
        if (i <= num_taxa - 2 - !rooted)
        {
            nodeLnFact += std::log(i);
        }
    }
    
    log_tree_topology_prob = (num_taxa - 2 - !rooted) * RbConstants::LN2 + nodeLnFact - branchLnFact;
    
    std::vector<Taxon> unordered_taxa = taxa;
    std::vector<std::string> ordered_taxa;
    for (size_t i = 0; i < unordered_taxa.size(); ++i)
    {
        ordered_taxa.push_back(unordered_taxa[i].getName());
    }
    
    // order taxon names
    std::sort(ordered_taxa.begin(), ordered_taxa.end());
    
    // pick an outgroup if this topology is unrooted and no outgroup was specified
    if ( outgroup.size() == 0 && rooted == false )
    {
        
        outgroup.addTaxon( ordered_taxa[0] );
        outgroup_provided = false;
        
    }
    
    std::map<Taxon, size_t> taxon_bitset_map;
    // add taxa to bitset map
    for (size_t i = 0; i < ordered_taxa.size(); ++i)
    {
        taxon_bitset_map[ordered_taxa[i]] = i;
    }
    
    RbBitSet b( num_taxa );
    for (size_t j = 0; j < outgroup.size(); j++)
    {
        size_t k = taxon_bitset_map[ outgroup.getTaxonName(j) ];
        
        b.set(k);
    }
    
    outgroup.setBitRepresentation( b );
    
    simulateTree();
    
    value->getTreeChangeEventHandler().addListener( this );
}


UniformTopologyBranchLengthDistribution::UniformTopologyBranchLengthDistribution(const UniformTopologyBranchLengthDistribution &d) : TypedDistribution<Tree>( d ),
    branch_length_prior( d.branch_length_prior->clone() ),
    num_taxa( d.num_taxa ),
    taxa( d.taxa ),
    log_tree_topology_prob( d.log_tree_topology_prob ),
    outgroup( d.outgroup ),
    outgroup_provided( d.outgroup_provided ),
    rooted( d.rooted ),
    dirty_topology( d.dirty_topology )
{
    
    // add the parameters of the distribution
    const std::vector<const DagNode*>& pars = branch_length_prior->getParameters();
    for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
    {
        this->addParameter( *it );
    }
    
    value->getTreeChangeEventHandler().addListener( this );
}



UniformTopologyBranchLengthDistribution::~UniformTopologyBranchLengthDistribution()
{
    
    delete branch_length_prior;
    // the tree will be deleted automatically by the base class
    
    value->getTreeChangeEventHandler().removeListener( this );
}


UniformTopologyBranchLengthDistribution& UniformTopologyBranchLengthDistribution::operator=(const UniformTopologyBranchLengthDistribution &d)
{
    
    if ( this != &d )
    {
        TypedDistribution<Tree>::operator=( d );
        
        // remove the old branch-length-prior parameters
        const std::vector<const DagNode*>& old_pars = branch_length_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = old_pars.begin(); it != old_pars.end(); ++it)
        {
            this->removeParameter( *it );
        }
        delete branch_length_prior;
        
        branch_length_prior         = d.branch_length_prior->clone();
        num_taxa                    = d.num_taxa;
        taxa                        = d.taxa;
        log_tree_topology_prob      = d.log_tree_topology_prob;
        outgroup                    = d.outgroup;
        outgroup_provided           = d.outgroup_provided;
        rooted                      = d.rooted;
        dirty_topology              = d.dirty_topology;

        // add the parameters of the distribution
        const std::vector<const DagNode*>& pars = branch_length_prior->getParameters();
        for (std::vector<const DagNode*>::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            this->addParameter( *it );
        }
    
    }

    return *this;
}


UniformTopologyBranchLengthDistribution* UniformTopologyBranchLengthDistribution::clone( void ) const
{
    
    return new UniformTopologyBranchLengthDistribution( *this );
}


double UniformTopologyBranchLengthDistribution::computeLnProbability( void )
{
    
    if ( dirty_topology == true )
    {
        size_t num_outgroup_species = outgroup.size();
        if ( num_outgroup_species > 0 )
        {
            // now we check that the outgroup is correct
            const TopologyNode &root = value->getRoot();
            const std::vector<TopologyNode*> &children = root.getChildren();
            bool contains_outgroup = false;
            for (size_t i=0; i<children.size(); ++i)
            {
                const TopologyNode &child = *(children[i]);
                size_t num_species_in_clade = child.getNumberOfNodesInSubtree(false) + 1;
                if ( num_outgroup_species == num_species_in_clade )
                {
                    Clade c = child.getClade();
                    if ( c == outgroup )
                    {
                        contains_outgroup = true;
                        break;
                    }
                
                }
            
            }
        
            if ( contains_outgroup == false )
            {
                return RbConstants::Double::neginf;
            }
        }
  
        dirty_topology = false;
        
    }
    
    const std::vector<TopologyNode*> &nodes = this->value->getNodes();
    size_t num_nodes = nodes.size();

    double ln_branch_length_prob = 0.0;
    for (size_t i=0; i<num_nodes; ++i)
    {
        const TopologyNode *node = nodes[i];
        if ( node->isRoot() == false )
        {
            double *bl = new double( value->getBranchLengthForNode(*node) );
            branch_length_prior->setValue( bl );
            ln_branch_length_prob += branch_length_prior->computeLnProbability();
        }
    }
    
    return ln_branch_length_prob + log_tree_topology_prob;
}


void UniformTopologyBranchLengthDistribution::fireTreeChangeEvent(const TopologyNode &n, const unsigned& m)
{
    
    if ( m == TreeChangeEventMessage::DEFAULT || m == TreeChangeEventMessage::TOPOLOGY )
    {
        
        dirty_topology = true;
    }
    
}


const std::vector<Taxon>& UniformTopologyBranchLengthDistribution::getTaxa( void ) const
{
    
    return taxa;
}


void UniformTopologyBranchLengthDistribution::redrawValue( void )
{
    value->getTreeChangeEventHandler().removeListener( this );
    simulateTree();
    value->getTreeChangeEventHandler().addListener( this );
    
    dirty_topology = true;
}


void UniformTopologyBranchLengthDistribution::setValue(RevBayesCore::Tree *v, bool force)
{

    value->getTreeChangeEventHandler().removeListener( this );
    
    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );
    
    value->getTreeChangeEventHandler().addListener( this );
    
    // Check that this isn't an artifact of arbitrary outgroup choice + clamping
    if ( rooted == false && outgroup_provided == false )
    {
        std::vector<Taxon> taxa;
        v->getRoot().getChild( 0 ).getTaxa( taxa );
        outgroup = Clade( taxa );
        
    }
    
    dirty_topology = true;
}


/**
 *
 */
void UniformTopologyBranchLengthDistribution::simulateClade(std::vector<TopologyNode *> &n )
{
    
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    while ( n.size() >= 2 )
    {
        
        // get all the nodes before the current age
        std::vector<TopologyNode*> active_nodes = n;
        
        // randomly pick two nodes
        size_t index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
        TopologyNode* left_child = active_nodes[index_left];
        active_nodes.erase(active_nodes.begin()+long(index_left));
        size_t index_right = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
        TopologyNode* right_right = active_nodes[index_right];
        active_nodes.erase(active_nodes.begin()+long(index_right));
        
        // erase the nodes also from the origin nodes vector
        n.erase(std::remove(n.begin(), n.end(), left_child), n.end());
        n.erase(std::remove(n.begin(), n.end(), right_right), n.end());
        
        
        // create a parent for the two
        TopologyNode *parent = new TopologyNode();
        parent->addChild( left_child );
        parent->addChild( right_right );
        left_child->setParent( parent );
        right_right->setParent( parent );
        
        // draw a branch length
        branch_length_prior->redrawValue();
        double bl = branch_length_prior->getValue();
        parent->setBranchLength( bl );
        
        // insert the parent to our list
        n.push_back( parent );
        
    }
    
    
}


void UniformTopologyBranchLengthDistribution::simulateTree( void )
{
    
    // the tree object
    //    Tree *psi = new Tree();
    Tree *psi = value;
    
    // internally we treat unrooted topologies the same as rooted
    psi->setRooted( rooted );
    
    // create the tip nodes
    std::vector<TopologyNode*> ingroup_nodes;
    std::vector<TopologyNode*> outgroup_nodes;
    std::vector<Taxon> ingroup_taxa;
    for (size_t i=0; i<num_taxa; ++i)
    {
        
        // create the i-th taxon
        TopologyNode* node = new TopologyNode( taxa[i], i );
        
        // draw a branch length
        branch_length_prior->redrawValue();
        double bl = branch_length_prior->getValue();
        node->setBranchLength( bl );
        
        bool is_outgroup = false;
        for (size_t j = 0; j < outgroup.size(); ++j)
        {
            if ( outgroup.getTaxon(j) == taxa[i] )
            {
                is_outgroup = true;
                break;
            }
        }
        
        // add the new node to the list
        if ( is_outgroup == true )
        {
            outgroup_nodes.push_back( node );
        }
        else
        {
            ingroup_nodes.push_back( node );
            ingroup_taxa.push_back( taxa[i] );
        }
        
    }
    
    if ( outgroup.size() != outgroup_nodes.size() )
    {
        throw RbException("Problem with the outgroup nodes. The names do not match.");
    }
    
    TopologyNode *root;
    
    // we need to simulate the ingroup
    simulateClade( ingroup_nodes );
    
    // check if using an outgroup
    if ( outgroup.size() > 0 )
    {
        // we need to simulate the outgroup
        simulateClade( outgroup_nodes );
        
        // now we create the root
        // we will take the ingroup root and the outgroup root
        TopologyNode *ingroup_root = ingroup_nodes[0];
        TopologyNode *outgroup_root = outgroup_nodes[0];
        
        if ( rooted == false )
        {
            // fix the root by connecting the ingroup with the outgroup
            root = ingroup_root;
            root->addChild( outgroup_root );
            outgroup_root->setParent( root );
        }
        else
        {
            // fix the root by connecting the ingroup and the outgroup to a new node
            root = new TopologyNode();
            root->addChild( ingroup_root );
            root->addChild( outgroup_root );
            ingroup_root->setParent( root );
            outgroup_root->setParent( root );
        }
    }
    // not using an outgroup
    else
    {
        if ( rooted == false )
        {
            // root the tree at the first non-tip child descending from the root
            TopologyNode *ingroup_root = ingroup_nodes[0];
            TopologyNode *left_child   = &(ingroup_root->getChild(0));
            TopologyNode *right_child  = &(ingroup_root->getChild(1));
            
            TopologyNode *new_child;
            if ( left_child->isTip() == true )
            {
                root = right_child;
                new_child = left_child;
            }
            else
            {
                root = left_child;
                new_child = right_child;
            }
            
            ingroup_root->removeChild(left_child);
            ingroup_root->removeChild(right_child);
            
            root->addChild(new_child);
            new_child->setParent(root);
            
            root->setParent(NULL);
            delete ingroup_root;
        }
        else
        {
            // simply set the root as the ingroup root
            root = ingroup_nodes[0];
        }
    }
    
    root->setBranchLength(0.0);
    
    // initialize the topology by setting the root
    psi->setRoot(root, true);
    
    // re-couple tip node names with tip indices
    // this is necessary because otherwise tip names get scrambled across replicates
    for (size_t i=0; i<num_taxa; i++)
    {
        psi->getTipNodeWithName(taxa[i].getName()).setIndex(i);
    }
    
    psi->orderNodesByIndex();
    
}


/** Swap a parameter of the distribution */
void UniformTopologyBranchLengthDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    if ( branch_length_prior != NULL )
    {
        branch_length_prior->swapParameter(oldP,newP);
    }
}

