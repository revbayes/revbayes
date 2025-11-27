#include "UniformTopologyDistribution.h"

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

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

UniformTopologyDistribution::UniformTopologyDistribution(const std::vector<Taxon> &tn, const Clade &og, const std::vector<Clade> &c, bool rt) : TypedDistribution<Tree>( new Tree() ),
    num_taxa( tn.size() ),
    taxa( tn ),
	constraints( c ),
    logTreeTopologyProb( RbConstants::Double::nan ),
    outgroup( og ),
    outgroup_provided( og.size() > 0 ),
	rooted( rt )
{
    
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
    
    logTreeTopologyProb = (num_taxa - 2 - !rooted) * RbConstants::LN2 + nodeLnFact - branchLnFact;
    
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

    for (size_t i = 0; i < constraints.size(); ++i)
    {
        RbBitSet b( num_taxa );
        for (size_t j = 0; j < constraints[i].size(); ++j)
        {
            size_t k = taxon_bitset_map[ constraints[i].getTaxonName(j) ];
            b.set(k);
            
        }

        constraints[i].setBitRepresentation( b );
    }

    simulateTree();
    
}


UniformTopologyDistribution::~UniformTopologyDistribution()
{
    // the tree will be deleted automatically by the base class
    
}


UniformTopologyDistribution* UniformTopologyDistribution::clone( void ) const
{

    return new UniformTopologyDistribution( *this );
}


double UniformTopologyDistribution::computeLnProbability( void )
{
    
	// first check if the current tree matches the clade constraints
    if ( matchesConstraints() == false )
    {
        return RbConstants::Double::neginf;
    }
    
    if ( outgroup.size() > 0 )
    {
		// now we check that the outgroup is correct
		const TopologyNode &root = value->getRoot();
		const std::vector<TopologyNode*> &children = root.getChildren();
		bool contains_outgroup = false;
		for (size_t i=0; i<children.size(); ++i)
		{
			const TopologyNode &child = *(children[i]);
			Clade c = child.getClade();
			if ( c == outgroup )
			{
				contains_outgroup = true;
				break;
			}
		}
        
		if ( contains_outgroup == false )
		{
			return RbConstants::Double::neginf;
		}
	}
    
    return logTreeTopologyProb;
}


void UniformTopologyDistribution::redrawValue( void )
{
    simulateTree();
}


void UniformTopologyDistribution::setValue(RevBayesCore::Tree *v, bool force)
{
    
    // delegate to super class
    TypedDistribution<Tree>::setValue( v, force );
    
    // Check that this isn't an artifact of arbitrary outgroup choice + clamping
    if ( rooted == false && outgroup_provided == false )
    {
        std::vector<Taxon> taxa;
        v->getRoot().getChild( 0 ).getTaxa( taxa );
        outgroup = Clade( taxa );
        
    }
    
}


/**
 *
 */
void UniformTopologyDistribution::simulateClade(std::vector<TopologyNode *> &n )
{
    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // The children contained in vector n have no parent.
    // Create new internal nodes to join 2 members of n together until n contains only 1 node.
    // Then use the single remaining node as the parent.
    while ( n.size() >= 2 )
    {
        
        // get all the nodes before the current age
        std::vector<TopologyNode*> active_nodes = n;
        
        // randomly pick two nodes
        size_t index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
        TopologyNode* left_child = active_nodes[index_left];
        active_nodes.erase(active_nodes.begin()+std::int64_t(index_left));
        size_t index_right = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
        TopologyNode* right_right = active_nodes[index_right];
        active_nodes.erase(active_nodes.begin()+std::int64_t(index_right));
                
        // erase the nodes also from the origin nodes vector
        n.erase(std::remove(n.begin(), n.end(), left_child), n.end());
        n.erase(std::remove(n.begin(), n.end(), right_right), n.end());
                
                
        // create a parent for the two
        TopologyNode *parent = new TopologyNode();
        parent->addChild( left_child );
        parent->addChild( right_right );
        left_child->setParent( parent );
        right_right->setParent( parent );
                
        // insert the parent to our list
        n.push_back( parent );
        
    }
    
    
}


void UniformTopologyDistribution::simulateTree( void )
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
    
    // we need a sorted vector of constraints, starting with the smallest
    std::vector<Clade> sorted_clades = {};
    for(int i = 0; i < constraints.size(); i++) //lyr changed
        if(constraints[i].isNegativeConstraint() == false)
            sorted_clades.push_back(constraints[i]);
    
    // create a clade that contains all species
    Clade all_ingroup_species = Clade(ingroup_taxa);
    sorted_clades.push_back(all_ingroup_species);
    
    // next sort the clades
    std::sort(sorted_clades.begin(),sorted_clades.end());
    
    // remove duplicates
    std::vector<Clade> tmp;
    tmp.push_back( sorted_clades[0] );
    for (size_t i = 1; i < sorted_clades.size(); ++i)
    {
        Clade &a = tmp[tmp.size()-1];
        Clade &b = sorted_clades[i];
        
        if ( a.size() != b.size() )
        {
            tmp.push_back( sorted_clades[i] );
        }
        else
        {
            bool equal = true;
            for (size_t j = 0; j < a.size(); ++j)
            {
                if ( a.getTaxon(j) != b.getTaxon(j) )
                {
                    equal = false;
                    break;
                }
            }
            if ( equal == false )
            {
                tmp.push_back( sorted_clades[i] );
            }
        }
        
    }
    sorted_clades = tmp;
    
    std::vector<Clade> virtual_taxa;
    for (size_t i = 0; i < sorted_clades.size(); ++i)
    {
        
        Clade &c = sorted_clades[i];
        std::vector<Taxon> taxa = c.getTaxa();
        std::vector<Clade> clades;
        
        for (int j = int(i)-1; j >= 0; --j)
        {
            const Clade &c_nested = sorted_clades[j];
            const std::vector<Taxon> &taxa_nested = c_nested.getTaxa();
            
            bool found_all = true;
            bool found_some = false;
            for (size_t k = 0; k < taxa_nested.size(); ++k)
            {
                std::vector<Taxon>::iterator it = std::find(taxa.begin(), taxa.end(), taxa_nested[k]);
                if ( it != taxa.end() )
                {
                    taxa.erase( it );
                    found_some = true;
                }
                else
                {
                    found_all = false;
                }
                
            }
            
            if ( found_all == true )
            {
                //                c.addTaxon( virtual_taxa[j] );
                //                taxa.push_back( virtual_taxa[j] );
                clades.push_back( virtual_taxa[j] );
            }
            
            if ( found_all == false && found_some == true )
            {
                throw RbException("Cannot simulate tree: conflicting monophyletic clade constraints. Check that all clade constraints are properly nested.");
            }
            
        }
        
        std::vector<TopologyNode*> nodes_in_clade;
        
        
        for (size_t k = 0; k < taxa.size(); ++k)
        {
            clades.push_back( Clade( taxa[k] ) );
        }
        
        for (size_t k = 0; k < clades.size(); ++k)
        {
            for (size_t j = 0; j < ingroup_nodes.size(); ++j)
            {
                if (ingroup_nodes[j]->getClade() == clades[k])
                {
                    nodes_in_clade.push_back( ingroup_nodes[j] );
                    ingroup_nodes.erase( ingroup_nodes.begin()+j );
                    break;
                }
                
            }
            
        }
        
        
        simulateClade( nodes_in_clade );
        ingroup_nodes.push_back( nodes_in_clade[0] );
        
        std::vector<Taxon> v_taxa;
        nodes_in_clade[0]->getTaxa(v_taxa);
        virtual_taxa.push_back( Clade(v_taxa) );
        
    }
    
    TopologyNode *root;
    
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
void UniformTopologyDistribution::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
}

/**
 * We check here if all the constraints are satisfied.
 * These are hard constraints, that is, the clades must be monophyletic.
 *
 * \return     True if the constraints are matched, false otherwise.
 */
bool UniformTopologyDistribution::matchesConstraints( void ) 
{
    
    if ( constraints.empty() == true )
    {
		return true;
	}
    else
    {
		
		const TopologyNode &root = value->getRoot();
		for (std::vector<Clade>::iterator it = constraints.begin(); it != constraints.end(); ++it) 
		{
            if ( root.containsClade( *it, true ) == false && it->isNegativeConstraint() == false ) //LYR changed
			{
				return false;
			}
            else if ( root.containsClade( *it, true ) == true && it->isNegativeConstraint() == true ) //LYR changed
            {
                return false;
            }
		}
        
		return true;
	}
    
}
