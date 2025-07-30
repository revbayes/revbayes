#include <cmath>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NewickConverter.h"
#include "RbException.h"
#include "Tree.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"
#include "TreeUtilities.h"
#include "Clade.h"
#include "DagNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbBitSet.h"
#include "RbBoolean.h"
#include "RbFileManager.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "TaxonMap.h"
#include "TreeChangeEventHandler.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }
namespace RevBayesCore { class TreeChangeEventListener; }

using namespace RevBayesCore;

/* Copy constructor */
Tree::Tree(const Tree& t) :
    changeEventHandler( ),
    root( NULL ),
    rooted( t.rooted ),
    is_negative_constraint( t.is_negative_constraint ),
    num_tips( t.num_tips ),
    num_nodes( t.num_nodes ),
    taxon_bitset_map( t.taxon_bitset_map )
{

    // need to perform a deep copy of the BranchLengthTree nodes
    if (t.root != NULL)
    {
        TopologyNode * newRoot = t.getRoot().clone();

        // set the root. This will also set the nodes vector.
        // do not reorder node indices when copying (WP)
        setRoot(newRoot, false);
    }
    
}


/* Copy constructor */
Tree::Tree(Tree&& t)
{
    operator=( std::move(t) );
}


/* Destructor */
Tree::~Tree(void)
{
    
    nodes.clear();
    std::set<TreeChangeEventListener*> l = changeEventHandler.getListeners();
    for ( std::set<TreeChangeEventListener*>::iterator it = l.begin(); it != l.end(); ++it )
    {
        changeEventHandler.removeListener( *it );
    }

    delete root;

}


Tree& Tree::operator=(const Tree &t)
{

    if (this != &t)
    {
        // nothing really to do here, should be done in the derived classes
        // @todo: Find a better solution - Sebastian
        // Problem: If we redraw the tree because the initial states are invalid,
        // then we somehow need to remember the tree event change listeners.
        // But it is not nice if the tree distribution needs to remember this!!!
//        changeEventHandler = t.changeEventHandler;



        nodes.clear();
        delete root;
        root = NULL;
        num_tips               = t.num_tips;
        num_nodes              = t.num_nodes;
        rooted                 = t.rooted;
        is_negative_constraint = t.is_negative_constraint;

        TopologyNode* newRoot = t.root->clone();

        // set the root. This will also set the nodes vector
        // do not reorder node indices when copying (WP)
        setRoot(newRoot, false);

    }

    return *this;
}

Tree& Tree::operator=(Tree&& t)
{
    // If there are listeners watching t, then maybe we shouldn't steal its nodes.
    // (However, since t is a temporary variable, then it probably shouldn't have listeners.)
    if (not t.changeEventHandler.getListeners().empty())
    {
        // Use the copy constructor in this case.
        return operator=(t);
    }

    // Otherwise steal nodes from t instead of copying them.
    // By swapping with t, we let t do the work of destroying our nodes.
    if (this != &t)
    {
        std::swap( rooted                , t.rooted                );
        std::swap( num_tips              , t.num_tips              );
        std::swap( num_nodes             , t.num_nodes             );
        std::swap( is_negative_constraint, t.is_negative_constraint);
        std::swap( root                  , t.root                  );
        std::swap( nodes                 , t.nodes                 );
        std::swap( taxon_bitset_map      , t.taxon_bitset_map      );

        // This loop is maybe a reason to NOT record a tree pointer on the nodes...
        for(auto& node: nodes)
            node->setTree(this);

        for(auto& node: t.nodes)
            node->setTree(&t);
    }

    return *this;
}

bool Tree::operator==(const Tree &t) const
{

    return getNewickRepresentation() == t.getNewickRepresentation();
}


bool Tree::operator!=(const Tree &t) const
{

    return !operator==(t);
}


bool Tree::operator<(const Tree &t) const
{

    return getNewickRepresentation() < t.getNewickRepresentation();
}


bool Tree::operator<=(const Tree &t) const
{

    return operator<(t) || operator==(t);
}


void Tree::addBranchParameter(std::string const &name, const std::vector<double> &parameters, bool internalOnly)
{

    getRoot().addBranchParameters(name,parameters,internalOnly);

}


void Tree::addNodeParameter(std::string const &name, const std::vector<double> &parameters, bool internalOnly)
{

    getRoot().addNodeParameters(name,parameters,internalOnly);

}

void Tree::addNodeParameter(std::string const &name, const std::vector<std::string> &parameters, bool internalOnly)
{

    getRoot().addNodeParameters(name,parameters,internalOnly);

}


void Tree::clearParameters( void )
{

    clearNodeParameters();
    clearBranchParameters();

}

void Tree::clearBranchParameters( void )
{

    getRoot().clearBranchParameters();

}


void Tree::clearNodeParameters( void )
{

    getRoot().clearNodeParameters();

}


/* Clone function */
Tree* Tree::clone(void) const
{

    return new Tree(*this);
}


void Tree::collapseNegativeBranchLengths(double l)
{
    double age = 0.0;
    double parent_age = 0.0;
    for (size_t i = 0; i < nodes.size(); i++)
    {
        if ( nodes[i]->isRoot() == false )
        {
            age = nodes[i]->getAge();
            parent_age = nodes[i]->getParent().getAge();

            if ( parent_age < age )
            {
                nodes[i]->setAge( parent_age - l, false );
            }
        }
    }

}

/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
bool Tree::containsClade(const TopologyNode &n, bool unrooted) const
{
    RbBitSet your_taxa = RbBitSet( getNumberOfTips() );
    n.getTaxa( your_taxa );

    bool contains = root->containsClade( your_taxa, true );

    if ( contains == false && unrooted == true )
    {
        your_taxa.flip();
        contains = root->containsClade( your_taxa, true );
    }

    return contains;
}


/**
 * Drop the tip node with the given name.
 * The name should correspond to the taxon name, not the species name.
 * This will throw an error if the name doesn't exist.
 */
void Tree::dropTipNodeWithName( const std::string &n )
{
    // get the index of this name
    size_t index = getTipIndex( n );
    dropTipNode( index );
}



/**
 * Drop the tip node with the given name.
 * The name should correspond to the taxon name, not the species name.
 * This will throw an error if the name doesn't exist.
 */
void Tree::dropTipNode( size_t index )
{
    // get the index of this name
    TopologyNode &node = getTipNode( index );
    if (node.isRoot() == true && nodes.size() == 1)
    {
        // there is nothing left to prune
        node.setName("");
        return;
    }
    
    TopologyNode &parent        = node.getParent();
    TopologyNode &grand_parent  = parent.getParent();
    
    /* If the parent of the current node has more than two children, we can just delete this child without affecting
     * the rest of the tree (i.e., we do not have to remove any internal nodes, nor adjust character history).
     */
    if (parent.getNumberOfChildren() > 2)
    {
        parent.removeChild( &node );
    }
    else
    {
        if (parent.isRoot() == false)
        {
            TopologyNode *sibling = &parent.getChild( 0 );
            if ( sibling == &node )
            {
                sibling = &parent.getChild( 1 );
            }
            grand_parent.removeChild( &parent );
            parent.removeChild( sibling );
            grand_parent.addChild( sibling );
            sibling->setParent( &grand_parent );

            // update character history
            if (parent.getTimeInStates().size() > 0 && sibling->getTimeInStates().size() > 0)
            {
                std::vector<double> sibling_state_times = sibling->getTimeInStates();
                for (size_t i = 0; i < parent.getTimeInStates().size(); i++)
                {
                    sibling_state_times[i] += parent.getTimeInStates()[i];
                }
                sibling->setTimeInStates(sibling_state_times);
                sibling->setNumberOfShiftEvents( sibling->getNumberOfShiftEvents() + parent.getNumberOfShiftEvents() );
            }
        }
        else
        {
            TopologyNode *sibling = &root->getChild( 0 );
            if ( sibling == &node )
            {
                sibling = &root->getChild( 1 );
            }
            root->removeChild(&node);
            sibling->setParent(NULL);
            root = sibling;
            if (root->getTimeInStates().size() > 0)
            {
                root->setTimeInStates(std::vector<double>(root->getTimeInStates().size(), 0.0));
                root->setNumberOfShiftEvents( 0 );
            }
        }
    }

    reindexNodes();

    // count the number of tips
    num_tips = 0;
    for (size_t i = 0; i < num_nodes; ++i)
    {
        if ( nodes[i] == NULL )
        {
            throw RbException("Problem while reading in tree.");
        }
        num_tips += ( nodes[i]->isTip() ? 1 : 0);
    }

}



void Tree::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, double &rv) const
{

    if ( n == "rootAge" )
    {
        rv = getRoot().getAge();
    }
    else if ( n == "branchLength" )
    {
        int index = (int)static_cast<const TypedDagNode<std::int64_t> *>( args[0] )->getValue()-1;
        rv = getNode( index ).getBranchLength();
    }
    else if ( n == "nodeAge" )
    {
        int index = (int)static_cast<const TypedDagNode<std::int64_t> *>( args[0] )->getValue()-1;
        rv = getNode( index ).getAge();
    }
    else if ( n == "treeLength" )
    {
        rv = getTreeLength();
    }
    else if (n == "gammaStatistic")
    {
        rv = RevBayesCore::TreeUtilities::getGammaStatistic( *this );
    }
    else if ( n == "meanInverseES" )
    {
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >* c = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( args[0] );
        size_t state_index = (size_t)static_cast<const TypedDagNode<int> *>( args[1] )->getValue();
        rv = RevBayesCore::TreeUtilities::getMeanInverseES( *this, c->getValue(), state_index );
    }
    else if ( n == "calculateMPD" )
    {
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >* c = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( args[0] );
        size_t state_index = (size_t)static_cast<const TypedDagNode<int> *>( args[1] )->getValue();
        size_t site_index = (size_t)static_cast<const TypedDagNode<int> *>( args[2] )->getValue() - 1;
        // why doesn't this work?
        //bool zscore = static_cast<const TypedDagNode<bool> *>( args[3] )->getValue();
        bool zscore = false;
        if (args[3]->getValueAsString() == "TRUE")
        {
            zscore = true;
        }
        bool branch_lengths = false;
        if (args[4]->getValueAsString() == "TRUE")
        {
            branch_lengths = true;
        }
        size_t reps = (size_t)static_cast<const TypedDagNode<int> *>( args[5] )->getValue();
        rv = RevBayesCore::TreeUtilities::calculateMPD(*this, c->getValue(), site_index, state_index, zscore, branch_lengths, reps);
    }
    else if ( n == "calculateMNTD" )
    {
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >* c = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( args[0] );
        size_t state_index = (size_t)static_cast<const TypedDagNode<int> *>( args[1] )->getValue();
        size_t site_index = (size_t)static_cast<const TypedDagNode<int> *>( args[2] )->getValue() - 1;
        // why doesn't this work?
        //bool zscore = static_cast<const TypedDagNode<bool> *>( args[3] )->getValue();
        bool zscore = false;
        if (args[3]->getValueAsString() == "TRUE")
        {
            zscore = true;
        }
        bool branch_lengths = false;
        if (args[4]->getValueAsString() == "TRUE")
        {
            branch_lengths = true;
        }
        size_t reps = (size_t)static_cast<const TypedDagNode<int> *>( args[5] )->getValue();
        rv = RevBayesCore::TreeUtilities::calculateMNTD(*this, c->getValue(), site_index, state_index, zscore, branch_lengths, reps);
    }
    else
    {
        throw RbException() << "A tree object does not have a member method called '" << n << "'.";
    }

}


void Tree::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, std::int64_t &rv) const
{

    if ( n == "parent" )
    {
        std::int64_t index = static_cast<const TypedDagNode<std::int64_t> *>( args[0] )->getValue()-1;
        if (getNode( index ).isRoot() == true)
        {
            throw RbException("Root has no parent.");
        }
        rv = std::int64_t( getNode( index ).getParent().getIndex() )+1;
    }
    else if ( n == "child" )
    {
        int parent_index = static_cast<const TypedDagNode<int> *>( args[0] )->getValue()-1;
        if (parent_index < num_tips)
        {
            throw RbException("Tips have no children");
        }
        int child_index = static_cast<const TypedDagNode<int> *>( args[1] )->getValue()-1;
        if (child_index < 0)
        {
            throw RbException("Child indices begin at 1.");
        }
        if (child_index > getNode( parent_index ).getNumberOfChildren() - 1)
        {
            throw RbException("Node does not have child with this index.");
        }
        rv = int( getNode( parent_index ).getChild(child_index).getIndex() )+1;
    }

    else if ( n == "numSampledAncestors")
    {
        rv = 0;
        for (size_t i=0; i< num_tips; i++)
        {
            rv += nodes[i]->isSampledAncestorTip();
        }
    }
    else if (n == "colless")
    {
        int s = 0;

        rv = RevBayesCore::TreeUtilities::getCollessMetric( getRoot(), s);
    }
    else if (n == "nnodes")
    {
        rv = nodes.size();
    }
    else if (n == "nbranches")
    {
        rv = nodes.size() - 1;
    }
    else if (n == "ntips")
    {
        rv = num_tips;
    }
    else if ( n == "fitchScore" )
    {
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >* c = static_cast<const TypedDagNode< AbstractHomologousDiscreteCharacterData >* >( args[0] );
        rv = RevBayesCore::TreeUtilities::getFitchScore( *this, c->getValue() );
    }
    else
    {
        throw RbException() << "A tree object does not have a member method called '" << n << "'.";
    }

}


void Tree::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, Boolean &rv) const
{

    if ( n == "isContainedInClade" )
    {
        int index = (int)static_cast<const TypedDagNode<std::int64_t> *>( args[0] )->getValue()-1;
        Clade clade = static_cast<const TypedDagNode<Clade> *>( args[1] )->getValue();
        clade.resetTaxonBitset( getTaxonBitSetMap() );

        if ( index < 0 || index >= nodes.size() )
        {
            std::stringstream s;
            s << "The index of the node must be between 1 and " << int(nodes.size()) << ".";
            throw RbException( s.str() );
        }


        size_t clade_index = 0;
        bool found = false;
        size_t min_clade_size = nodes.size() + 2;
        size_t taxa_count = clade.size();

        for (size_t i = getNumberOfTips(); i < nodes.size(); ++i)
        {

            TopologyNode *node = nodes[i];
            size_t clade_size = size_t( (node->getNumberOfNodesInSubtree(true) + 1) / 2);
            if ( clade_size < min_clade_size && clade_size >= taxa_count && node->containsClade( clade, false ) )
            {
                found = true;
                clade_index = node->getIndex();
                min_clade_size = clade_size;
                if ( taxa_count == clade_size )
                {
                    break;
                }

            }

        }

        if ( found == true )
        {
            while ( index != clade_index && nodes[index]->isRoot() == false )
            {
                index = int( nodes[index]->getParent().getIndex() );
            }

        }

        rv = Boolean( index == clade_index );
    }
    else if ( n == "containsClade" )
    {
        Clade clade = static_cast<const TypedDagNode<Clade> *>( args[0] )->getValue();
        clade.resetTaxonBitset( getTaxonBitSetMap() );

        bool contained = root->containsClade(clade, true);

        rv = Boolean( contained );
    }
    else if ( n == "hasSameTopology" )
    {
        const Tree& t = static_cast<const TypedDagNode<Tree> *>( args[0] )->getValue();

        rv = Boolean( hasSameTopology(t) );
    }
    else
    {
        throw RbException() << "A tree object does not have a member method called '" << n << "'.";
    }

}

/* fill the nodes vector by a phylogenetic traversal recursively starting with this node.
 * The tips fill the slots 0,...,n-1 followed by the internal nodes and then the root.
 */
void Tree::fillNodesByPhylogeneticTraversal(TopologyNode* node)
{
    // now call this function recursively for all your children
    for (size_t i=0; i<node->getNumberOfChildren(); i++)
    {
        fillNodesByPhylogeneticTraversal(&node->getChild(i));
    }

    if (node->isTip())
    {
        // all the tips go to the beginning
        nodes.insert(nodes.begin(), node);
    }
    else
    {
        // this is phylogenetic ordering so the internal nodes come last
        nodes.push_back(node);
    }
}

const std::vector<std::vector<double> > Tree::getAdjacencyMatrix(void) const
{
    std::vector<std::vector<double> > adjacency(num_nodes, std::vector<double>(num_nodes, 0.0));

    for (size_t i = 0; i < nodes.size(); i++)
    {
        const TopologyNode* nd = nodes[i];
        std::vector<TopologyNode*> children = nd->getChildren();
        for (std::vector<TopologyNode*>::iterator ch = children.begin(); ch != children.end(); ch++)
        {
            adjacency[nd->getIndex()][(*ch)->getIndex()] = (*ch)->getBranchLength();
        }
        if (!nd->isRoot())
        {
            const TopologyNode* pa = &nd->getParent();
            adjacency[nd->getIndex()][pa->getIndex()] = pa->getBranchLength();
        }
    }

    return adjacency;
}


std::vector<Taxon> Tree::getFossilTaxa() const
{
    std::vector< Taxon > taxa;
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        if ( n.isFossil() == true )
        {
            taxa.push_back( n.getTaxon() );
        }

    }
    std::sort(taxa.begin(), taxa.end());

    return taxa;
}


/** We provide this function to allow a caller to randomly pick one of the interior nodes.
 This version assumes that the root is always the last and the tips the first in the nodes vector. */
const TopologyNode& Tree::getInteriorNode( size_t indx ) const
{

    size_t n = getNumberOfTips();
    if ( indx > (n-2) )
    {
        throw RbException() << "Cannot acces interior node '" << StringUtilities::to_string(indx) << "' for a tree with " << StringUtilities::to_string(n) << " tips.";
    }
    return *nodes[ indx + n ];
}


TopologyNode& Tree::getMrca(const Clade &c)
{

    return *(root->getMrca( c ));
}


const TopologyNode& Tree::getMrca(const Clade &c) const
{

    return *(root->getMrca( c ));
}

const TopologyNode& Tree::getMrca(const Clade &c, bool strict) const
{

    return *(root->getMrca( c,strict ));
}


const TopologyNode& Tree::getMrca(const TopologyNode &n) const
{

    return *(root->getMrca( n ));
}


std::string Tree::getNewickRepresentation(bool round ) const
{

    if ( root == NULL )
    {
        return "";
    }
    else
    {
        return root->computeNewick( round );
    }

}



TopologyNode& Tree::getNode(size_t idx)
{

    if ( idx >= nodes.size() )
    {
        throw RbException("Index out of bounds in getNode.");
    }

    return *nodes[idx];
}


const TopologyNode& Tree::getNode(size_t idx) const
{

    if ( idx >= nodes.size() )
    {
        throw RbException("Index out of bounds in getNode.");
    }

    return *nodes[idx];
}


const std::vector<TopologyNode*>& Tree::getNodes(void) const
{

    return nodes;
}


std::vector<RbBitSet>* Tree::getNodesAsBitset(void) const
{

    std::vector<RbBitSet>* bs = new std::vector<RbBitSet>();

    for ( size_t i=0; i<nodes.size(); ++i )
    {
        TopologyNode *n = nodes[i];
        if ( n->isTip() == false )
        {
            RbBitSet taxa_this_node = RbBitSet(num_tips);
            n->getTaxa(taxa_this_node);
            bs->push_back( taxa_this_node );
        }
    }

    return bs;
}


std::map<RbBitSet, TopologyNode*> Tree::getBitsetToNodeMap(void) const
{

    std::map<RbBitSet, TopologyNode*> bstn;

    for ( size_t i=0; i<nodes.size(); ++i )
    {
        TopologyNode *n = nodes[i];
        if ( n->isTip() == false )
        {
            RbBitSet taxa_this_node = RbBitSet(num_tips);
            n->getTaxa(taxa_this_node);
            bstn[taxa_this_node] = n;
        }
    }

    return bstn;
}


std::vector<std::int64_t> Tree::getNodeIndices(void) const
{
    std::vector<std::int64_t> indices;

    for ( size_t i=0; i<nodes.size(); ++i )
    {
        indices.push_back(nodes[i]->getIndex() );
    }

    return indices;
}




/**
 * Calculate the number of interior nodes in the BranchLengthTree by deducing the number of
 * tips from number of nodes, and then subtract 1 more if the BranchLengthTree is rooted.
 */
size_t Tree::getNumberOfInteriorNodes( void ) const
{

    size_t preliminaryNumIntNodes = getNumberOfNodes() - getNumberOfTips();

    if ( isRooted() )
    {
	if (preliminaryNumIntNodes > 1)
	    return preliminaryNumIntNodes - 1;
	else
	    return 0;
    }
    else
    {
        return preliminaryNumIntNodes;
    }

}


size_t Tree::getNumberOfNodes(void) const
{

    return num_nodes;
}


/**
 * return the number of tips.
 */
size_t Tree::getNumberOfTips( void ) const
{

    return num_tips;
}


/**
 * return the number of tips.
 */
size_t Tree::getNumberOfExtinctTips( void ) const
{
    size_t num_extinct = 0;
    for (size_t i = 0; i < num_tips; i++)
    {
        num_extinct += nodes[i]->isFossil();
    }

    return num_extinct;
}


/**
 * return the number of tips.
 */
size_t Tree::getNumberOfExtantTips( void ) const
{
    return num_tips - getNumberOfExtinctTips();
}


/**
 * return the number of tips.
 */
size_t Tree::getNumberOfSampledAncestors( void ) const
{
    size_t num_sa = 0;
    for (size_t i = 0; i < num_tips; i++)
    {
        num_sa += nodes[i]->isSampledAncestorTip();
    }

    return num_sa;
}


std::string Tree::getPlainNewickRepresentation() const
{

    return root->computePlainNewick();
}


TopologyNode& Tree::getRoot(void)
{
    return *root;
}


const TopologyNode& Tree::getRoot(void) const
{
    return *root;
}


/**
 * Get the tree and character history in newick format
 * compatible with SIMMAP and phytools
 */
std::string Tree::getSimmapNewickRepresentation(bool round) const
{
    return root->computeSimmapNewick(round);
}


/**
 * Get all the species names for this topology.
 * This might include duplicates.
 */
std::vector<std::string> Tree::getSpeciesNames() const
{
    std::vector< std::string > snames;
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        snames.push_back( n.getTaxon().getSpeciesName() );
    }

    return snames;
}

std::vector<Taxon> Tree::getTaxa() const
{
    std::vector< Taxon > taxa;
    for (auto& node: nodes)
    {
        if (node->isTip() or node->getNumberOfChildren() == 1)
        {
            Taxon taxon = node->getTaxon();
            taxon.setAge(node->getAge());
            taxa.push_back( taxon );
        }
    }

    return taxa;
}


/**
 * Check that all tip ages are within their specified taxon age ranges
 * If forceAdjust is true (default), will adjust misspecified ages, otherwise throws an exception
 */
void RevBayesCore::Tree::checkTaxonAges(bool forceAdjust)
{
    for (auto& node : nodes)
    {
        if(!node->isTip()) continue;
        Taxon taxon = node->getTaxon();
        if(node->getAge() < taxon.getMinAge()) {
            if(forceAdjust) {
                node->setAge(taxon.getMinAge());
                RBOUT("Age of taxon " + taxon.getName() + " was below the specified minimum and has been adjusted.");
            } else throw RbException() << "Age of taxon " << taxon.getName() << " is below the specified minimum.";
        }
        if(node->getAge() > taxon.getMaxAge()) {
            if(forceAdjust) {
                node->setAge(taxon.getMaxAge());
                RBOUT("Age of taxon " + taxon.getName() + " was above the specified maximum and has been adjusted.");
            } else throw RbException() << "Age of taxon " << taxon.getName() << " is above the specified maximum.";
        }
    }
}


/**
 * Returns a map of the taxa to their BitSet indices.
 * The taxa are ordered alphabetically in the BitSet.
 * Eventually this should be refactored with the TaxonMap class.
 */
const std::map<std::string, size_t> &Tree::getTaxonBitSetMap(void) const
{
    if (taxon_bitset_map.size() == 0)
    {
        // get all taxon names
        std::vector<Taxon> unordered_taxa = getTaxa();
        std::vector<std::string> ordered_taxa;
        for (size_t i = 0; i < unordered_taxa.size(); ++i)
        {
            ordered_taxa.push_back(unordered_taxa[i].getName());
        }

        // order taxon names
        std::sort(ordered_taxa.begin(), ordered_taxa.end());

        // add taxa to bitset map
        for (size_t i = 0; i < ordered_taxa.size(); ++i)
        {
            taxon_bitset_map[ordered_taxa[i]] = i;
        }
    }
    
    return taxon_bitset_map;
}

/**
 * Get the tip index for this name.
 */
size_t Tree::getTipIndex( const std::string &name ) const
{
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        if ( name == n.getName() )
        {
            return n.getIndex();
        }
    }

    // if name not found
    throw RbException() << "Could not find tip node with name '" << name << "' in tree." ;
}


std::vector<std::string> Tree::getTipNames() const
{

    std::vector<std::string> names;
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        const TopologyNode& n = getTipNode( i );
        names.push_back( n.getName() );
    }

    return names;
}


/**
 * We provide this function to allow a caller to randomly pick one of the tip nodes.
 * This version assumes that the tips are first in the nodes vector.
 */
TopologyNode& Tree::getTipNode( size_t index )
{

//
//    if ( index >= getNumberOfTips() )
//    {
//        throw RbException("Index out of bounds in getTipNode()!");
//    }
//    if (!nodes[ index ]->isTip())
//    {
//        throw RbException("Node at index is not a tip but should have been!");
//    }

    return *nodes[ index ];
}


const TopologyNode& Tree::getTipNode(size_t index) const
{

    //
    //    if ( index >= getNumberOfTips() )
    //    {
    //        throw RbException("Index out of bounds in getTipNode()!");
    //    }
    //    if (!nodes[ index ]->isTip())
    //    {
    //        throw RbException("Node at index is not a tip but should have been!");
    //    }

    return *nodes[index];
}


/**
 * Get the tip node with the given name.
 * The name should correspond to the taxon name, not the species name.
 * This will throw an error if the name doesn't exist.
 */
TopologyNode& Tree::getTipNodeWithName( const std::string &n )
{
    // get the index of this name
    size_t index = getTipIndex( n );

    return *nodes[ index ];
}


/**
 * Get the tip node with the given name.
 * The name should correspond to the taxon name, not the species name.
 * This will throw an error if the name doesn't exist.
 */
const TopologyNode& Tree::getTipNodeWithName( const std::string &n ) const
{
    // get the index of this name
    size_t index = getTipIndex( n );

    return *nodes[ index ];
}



/**
 * Get all the tip nodes with this species name.
 * If there is none, then we return an empty vector.
 */
std::vector<TopologyNode*> Tree::getTipNodesWithSpeciesName( const std::string &name )
{
    // create the vector of the tip nodes with this species name
    std::vector<TopologyNode*> tipNodes;

    // loop over all tips
    for (size_t i = 0; i < getNumberOfTips(); ++i)
    {
        // get the i-th tip
        TopologyNode& n = getTipNode( i );

        // test if the species name matches
        if ( name == n.getSpeciesName() )
        {
            // add this tip node to our list
            tipNodes.push_back( &n );
        }
    }

    // return the vector
    return tipNodes;
}


double Tree::getTmrca(const Clade &c)
{

    return root->getTmrca( c );
}


double Tree::getTmrca(const TopologyNode &n)
{

    return root->getTmrca( n );
}


double Tree::getTmrca(const std::vector<Taxon> &t)
{

    return root->getTmrca( t );
}


TreeChangeEventHandler& Tree::getTreeChangeEventHandler( void ) const
{

    return changeEventHandler;
}


double Tree::getTreeLength( void ) const
{

    double tl = 0.0;
    // loop over all nodes
    for (size_t i = 0; i < num_nodes; ++i)
    {
        // get the i-th node
        const TopologyNode& n = *nodes[i];

        if ( n.isRoot() == false )
        {
            // add the branch length
            tl += n.getBranchLength();
        }

    }

    return tl;
}


bool Tree::hasSameTopology(const Tree &t) const
{

    std::string a = getPlainNewickRepresentation();
    std::string b = t.getPlainNewickRepresentation();

    return a == b;
}


// Serialize (resurrect) the object from a file
void Tree::initFromFile( const path &dir, const std::string &fn )
{
    path filename = dir / (fn + ".newick");

    // open the stream to the file
    std::ifstream inStream( filename.string() );

    std::string s = "";
    while ( inStream.good() )
    {
        // Read a line
        std::string line;
        safeGetline( inStream, line );

        // append
        s += line;

    }

    return initFromString( s );
}


void Tree::initFromString(const std::string &s)
{
    NewickConverter converter;
    Tree* bl_tree = converter.convertFromNewick( s );
    if (isTimeTree())
    {
        auto time_tree = TreeUtilities::convertTree(*bl_tree, false);
        delete bl_tree;
        *this = *time_tree;
        delete time_tree;
    }
    else
    {
        *this = *bl_tree;
        delete bl_tree;
    }
}


bool Tree::isBinary(void) const
{
    for (size_t i = 0; i < getNumberOfInteriorNodes(); ++i)
    {
        const TopologyNode &n = getInteriorNode( i );
        if (n.getNumberOfChildren() != 2)
        {
            return false;
            break;
        }
    }
    return true;
}


bool Tree::isBroken( void ) const
{

    for (size_t i = 0; i < getNumberOfInteriorNodes(); ++i)
    {

        const TopologyNode &n = getInteriorNode( i );
        double age = n.getAge();

        for (size_t j = 0; j < n.getNumberOfChildren(); ++j)
        {
            const TopologyNode &child = n.getChild( j );

            double est_age = child.getAge() + child.getBranchLength();

            if ( std::fabs(age-est_age) > 1E-4 )
            {
                return true;
            }

        }

    }

    for (size_t i = 0; i < getNumberOfNodes(); ++i)
    {

        const TopologyNode &n = getNode( i );

        if ( n.isRoot() == false )
        {
            double my_age = n.getAge();
            double my_parents_age = n.getParent().getAge();

            if ( std::fabs( my_parents_age - my_age - n.getBranchLength() ) > 1E-4 )
            {
                return true;
            }

        }

    }


    return false;
}

bool Tree::isNegativeConstraint(void) const
{

    return is_negative_constraint;
}

bool Tree::isRooted(void) const
{

    return rooted;
}

bool Tree::isTimeTree(void) const
{
    bool any_node_times = false;
    bool all_node_times = true;

    for(auto& node: nodes)
    {
        if (not std::isnan(node->getAge()))
            any_node_times = true;
        else
            all_node_times = false;
    }

    if (any_node_times and not all_node_times)
        throw RbException()<<"isTimeTree: only some nodes have node times!";

    return all_node_times;
}


void Tree::makeRootBifurcating(const Clade& outgroup)
{
    size_t num_root_children = root->getNumberOfChildren();
    if ( num_root_children == 3 )
    {
        // which one of the 3 children do we want as outgroup?
        size_t good_outgroup = 0;
        double good_bl = 0.0;
        double half_bl = 0.0;
        std::vector<TopologyNode*> nodes_to_move = root->getChildren() ;
        for (size_t i = 0; i < num_root_children; ++i)
        {
            if ( root->getChild(i).getClade() == outgroup )
            {
                good_outgroup = i;
                good_bl = root->getChild(i).getBranchLength();
                half_bl = good_bl / 2.0;
                root->getChild(i).setBranchLength(half_bl);
                //  root->getChild(i).setAge(root->getChild(i).getAge() + half_bl);
                break;
            }
        }
        nodes_to_move.erase(nodes_to_move.begin() + good_outgroup);
          
        TopologyNode *new_child = new TopologyNode();
        double parent_age = root->getAge();
  
        double new_age = parent_age - half_bl;
    
        for (size_t i = 0; i < nodes_to_move.size(); ++i)
        {
            TopologyNode* tmp = nodes_to_move[i];
            //  double newBl = tmp->getBranchLength() - halfMini;
            root->removeChild(tmp);
            new_child->addChild(tmp);
            tmp->setParent (new_child);
            //  tmp->setBranchLength(newBl) ;
        }

        root->addChild(new_child);
        new_child->setParent( root );
        new_child->setAge(new_age);
        new_child->setBranchLength( half_bl );
    }
    else if (num_root_children > 2)
    {
        throw RbException() << "Cannot reroot on '" << outgroup.toString() << "'; the root has more than 3 children.";
    } // end-if the root node has 3 children
    

    // we need to reset the root so that the vector of nodes get filled again with the new number of nodes
    // it only makes sense to reindex the nodes because we have more nodes now!
    bool reindex = true;
    setRoot( &getRoot(), reindex );

    // clear the taxon bitset map
    // the next time someone call getTaxonBitset() it will be rebuilt
    taxon_bitset_map.clear();

} 


bool Tree::tryReadIndicesFromParameters(bool remove)
{
    // This routine assumes that the nodes array already contains the tree nodes.

    // If we want to allow for nodes NOT having valid indices, we could do
    // `setupIndices( not *has_indices )` instead of conditionally doing orderNodesByIndex().

    std::optional<bool> has_indices;

    std::vector<TopologyNode*> nodes2(nodes.size(), nullptr);

    std::optional<std::string> failure;

    for (auto& node: nodes)
    {
        // 1. Get the index as a string, optionally erasing it from the parameters.
        std::optional<std::string> value;
        if (remove)
            value = node->eraseNodeParameter("index");
        else
            value = node->getNodeParameter("index");

        try
        {
            // 2. Does this node have an index?
            bool has_index = value.has_value();

            // 3. Complain if some nodes have indices and some do not.
            if (not has_indices)
                has_indices = has_index;

            if (has_index != *has_indices)
                throw RbException()<<"readIndexFromParameter: Some nodes have an index parameter, and some do not!";

            // 4. Set the index if there is one.
            if (value)
            {
		int index = -1;
		try
		{
		    index = stoi(*value) - 1;
		}
		catch (...)
		{
		    throw RbException()<<"tryReadIndicesFromParameters: index '"<<*value<<"' is not an integer";
		}

                if (index < 0)
                    throw RbException()<<"tryReadIndicesFromParameters: Negative node index "<<index;
                else if (index >= nodes2.size())
                    throw RbException()<<"tryReadIndicesFromParameters: Tree only has "<<nodes.size()<<" nodes, but got index "<<index;
                else if (nodes2[index] != nullptr)
                    throw RbException()<<"tryReadIndicesFromParameters: Node index "<<index<<" used twice!";

                nodes2[index] = node;
            }
        }
        // NB - If remove == true, then we can't quit.  We have to finish erasing the "index" parameter from all nodes.
        catch (RbException& e)
        {
            failure = e.getMessage();
        }
	catch (std::exception& e)
	{
	    failure = e.what();
	}
	catch (...)
	{
	    failure = "tryReadIndicesFromParameters: unhandled exception";
	}
    }

    // 5. Report failure message.
    if (failure)
    {
        RBOUT("Warning: " + *failure + ".  Ignoring indices.");

        has_indices = false;
    }

    // 6.If the tree is empty, set has_indices to false.
    if (not has_indices.has_value())
        has_indices = false;


    // 7. If we set the indices
    if (*has_indices)
    {
	// If we get here, we know that every node had a unique index.

        nodes = nodes2;

        for(int i=0;i < nodes.size(); i++)
            nodes[i]->setIndex(i);
    }
    else
    {
        // We assume the nodes have valid indices that we can fall back on if reading indices fails.
        for(int i=0;i<nodes.size();i++)
            assert(nodes[i]->getIndex() == i);
    }

    // 8. Did we successfully set indices from parameters?
    return *has_indices;
}


void Tree::writeIndicesToParameters()
{
    for (auto& node: nodes)
    {
        node->setNodeParameter("index",std::to_string(node->getIndex()+1));
    }
}


// method to order nodes by their existing index
// used when reading in tree with existing node indexes we need to keep
void Tree::orderNodesByIndex( void )
{


    std::vector<TopologyNode*> nodes_copy = std::vector<TopologyNode*>(nodes.size());
    std::vector<bool> used = std::vector<bool>(nodes.size(),false);
    for (int i = 0; i < nodes.size(); i++)
    {
        if ( nodes[i]->getIndex() >= nodes.size() )
        {
            throw RbException() << "Problem while working with tree: Node had bad index. Index was '" << nodes[i]->getIndex() << "' while there are only '" << nodes.size() << "' nodes in the tree.";
        }
        else if ( used[nodes[i]->getIndex()] == true )
        {
            throw RbException() << "Problem while working with tree: Node had bad index. Two nodes had same index of '" << nodes[i]->getIndex() << "'.";
        }
        else
        {
            used[nodes[i]->getIndex()] = true;
        }
        nodes_copy[ nodes[i]->getIndex() ] = nodes[i];
    }

    nodes = nodes_copy;

}

// Prints tree for user (rounding)
void Tree::printForUser( std::ostream &o, const std::string &sep, int l, bool left ) const {
    std::int64_t previousPrecision = o.precision();
    std::ios_base::fmtflags previousFlags = o.flags();

    std::fixed( o );
    o.precision( 3 );

    o << *this;
    
    o.setf( previousFlags );
    o.precision( previousPrecision );
}

// Prints tree for simple storing (rounding)
void Tree::printForSimpleStoring( std::ostream &o, const std::string &sep, int l, bool left, bool flatten ) const {
    std::stringstream ss;
    ss << *this;
    std::string s = ss.str();
    if ( l > 0 ) 
    {
        StringUtilities::fillWithSpaces(s, l, left);
    }
    o << s;
}

// Prints tree for complex storing (no rounding; mostly checkpointing)
void Tree::printForComplexStoring ( std::ostream &o, const std::string &sep, int l, bool left, bool flatten ) const {
    std::stringstream ss;

    // call getNewickRepresentation with round == FALSE
    ss << this->getNewickRepresentation( false );
    std::string s = ss.str();
    if (l > 0) {
        StringUtilities::fillWithSpaces(s, l, left);
    }
    o << s;
}

json Tree::toJSON() const
{
    return this->getNewickRepresentation(false);
}

void Tree::collapseSampledAncestors()
{
    for(auto& node: nodes)
    {
        if (node->isSampledAncestorTip())
        {
            node->getParent().setTaxon(node->getTaxon());
            node->getParent().removeChild(node);
        }
    }
    setRoot(root, true);
}

// Pick a random node which is not the root, a tip, or the parent of a sampled ancestor.
// First try doing so at random; if that does not work, use the more computationally demanding strategy of finding all eligible nodes.
TopologyNode* Tree::pickRandomInternalNode(RandomNumberGenerator* rng) const
{
    TopologyNode* node;
    double u = rng->uniform01();
    
    for (size_t i = 0; i < 10; i++)
    {
        size_t node_idx = size_t( std::floor(getNumberOfNodes() * u) );
        node = (TopologyNode*)&getNode(node_idx);
        if ( !node->isRoot() && !node->isTip() && !node->isSampledAncestorParent() ) return node;
    }
    
    // check that there is at least one node which is not the root, a tip, or the parent of a SA
    std::vector<TopologyNode*> eligible_nodes;
        
    for (auto& to_check: getNodes())
    {
        if ( !to_check->isRoot() && !to_check->isTip() && !to_check->isSampledAncestorParent() )
        {
            eligible_nodes.push_back( to_check );
        }
    }
        
    if (eligible_nodes.size() == 0)
    {
        node = NULL;
    }
    else
    {
        size_t node_idx = size_t( std::floor(eligible_nodes.size() * u) );
        node = eligible_nodes[node_idx];
    }
    
    return node;
}

void Tree::pruneTaxa(const RbBitSet& prune_map )
{
    nodes.clear();

    // bootstrap all nodes from the root and add the in a pre-order traversal
    recursivelyPruneTaxa(root, prune_map);

    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->setIndex(i);
    }
    
    // we also need to reset our internal variables
    num_nodes = nodes.size();
    num_tips = (num_nodes+1) / 2;
}

bool Tree::recursivelyPruneTaxa( TopologyNode* n, const RbBitSet& prune_map )
{
    if( n->isTip()  )
    {
        bool prune = prune_map[n->getIndex()];

        if( prune == false )
        {
            nodes.insert(nodes.begin(), n);
        }

        return prune;
    }

    std::vector<TopologyNode*> children = n->getChildren();

    //std::vector<TopologyNode*> retained_children;
    std::vector<TopologyNode*> pruned_children;
    for(size_t i = 0; i < children.size(); i++)
    {
        if( recursivelyPruneTaxa(children[i], prune_map) )
        {
            pruned_children.push_back(children[i]);
        }
    }

    // if we don't prune any children, then add this node to the list of keepers
    if( pruned_children.empty() )
    {
        nodes.push_back(n);
    }

    // if we prune all or zero children, then just continue up the tree
    if( pruned_children.empty() || pruned_children.size() == children.size() )
    {
        return pruned_children.size() == children.size();
    }

    // if we prune some, but not all children, then patch over this node

    // first prune the dead branches
    for(size_t i = 0; i < pruned_children.size(); i++)
    {
        n->removeChild(pruned_children[i]);
        pruned_children[i]->setParent(NULL);
        delete pruned_children[i];
    }

    // if there are still enough retained root children, then return
    if( n->isRoot() && rooted == false && children.size() - pruned_children.size() < 3 )
    {
        // there are not enough retained children
        // patch up the node
        std::vector<TopologyNode*> retained_children = n->getChildren();

        for(size_t i = 0; i < retained_children.size(); i++)
        {
            n->removeChild(retained_children[i]);
            retained_children[i]->setParent(NULL);
        }

        delete n;

        if( retained_children.size() == 1 )
        {
            n = retained_children.back();

            retained_children = n->getChildren();
            for(size_t i = 0; i < retained_children.size(); i++)
            {
                n->removeChild(retained_children[i]);
                retained_children[i]->setParent(NULL);
            }

            delete n;
        }

        root = retained_children.back();
        if(root->isTip()) { // cannot use tip as the new root, so pick the other one
            root = retained_children.front();
            root->addChild(retained_children.back());
            retained_children.back()->setParent(root);
        } else { // otherwise, proceed normally
            root->addChild(retained_children.front());
            retained_children.front()->setParent(root);
        }
    }
    // if there are still at least 2 retained children, then return
    else if( children.size() - pruned_children.size() < 2 )
    {
        // there is only one retained child
        // patch up the 2-degree node
        TopologyNode* retained_child = &n->getChild(0);

        n->removeChild(retained_child);

        if( n->isRoot() )
        {
            retained_child->setParent(NULL);
            root = retained_child;
        }
        else
        {
            TopologyNode* parent = &n->getParent();

            parent->removeChild(n);
            parent->addChild(retained_child);
            retained_child->setParent(parent);
            n->setParent(NULL);
        }

        delete n;
    }

    return false;
}


void Tree::reindexNodes()
{
    nodes.clear();

    // bootstrap all nodes from the root and add the in a pre-order traversal
    fillNodesByPhylogeneticTraversal(root);
    for (unsigned int i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->setIndex(i);
    }

    num_nodes = nodes.size();
}

void Tree::removeDegree2Node(TopologyNode* n)
{
    // 1. Check that the node is in fact degree 2
    if (n->getDegree() != 2)
        throw RbException()<<"removeDegree2Node: node has degree "<<n->getDegree()<<"!";

    // 2. If n is the root, pivot so that it is a child of the root
    TopologyNode* new_root = nullptr;
    if (n->isRoot())
    {
        // Pick the child with the highest degree for the new root
        TopologyNode* new_root = &n->getChild(0);
        for(auto& child: n->getChildren())
            if (child->getDegree() > new_root->getDegree())
                new_root = child;
        reverseParentChild( *new_root );

        root = new_root;
        root->setTree(this);

        assert(root->isRoot());
    }


    // 3. Save references from n.
    //    Now we have p -> n -> c
    assert(not n->isRoot());
    auto& p = n->getParent();
    auto& c = n->getChild(0);
    double combined_branch_length = n->getBranchLength() + c.getBranchLength();

    // 4. Remove references from n
    n->setParent(nullptr);
    n->removeChild(&c);

    // 5. Link parent and child and set new branch length
    p.removeChild(n);
    p.addChild(&c);
    c.setParent(&p);
    c.setBranchLength(combined_branch_length);

    // 6. We need to reindex nodes because we removed an index.
    reindexNodes();

    // 7. Destroy the removed node.
    delete n;
}

bool Tree::removeNodeIfDegree2(TopologyNode& n)
{
    if (n.getDegree() == 2)
    {
        removeDegree2Node( &n );
        return true;
    }
    else
        return false;
}

bool Tree::removeRootIfDegree2()
{
    return removeNodeIfDegree2(getRoot());
}

void Tree::removeDuplicateTaxa( void )
{

    bool removed_replicate = true;
    while ( removed_replicate == true )
    {
        removed_replicate = false;
        for ( size_t i=0; i<(num_tips-1); ++i )
        {
            const std::string &name_a = nodes[ i ]->getName();
            for ( size_t j=i+1; j<num_tips; ++j )
            {
                const std::string &name_b = nodes[ j ]->getName();
                if ( name_a == name_b )
                {
                    removed_replicate = true;
                    dropTipNode( j );
                    break;
                }

            }

            if ( removed_replicate == true )
            {
                break;
            }

        }

    }

}


void Tree::renameNodeParameter(const std::string &old_name, const std::string &new_name)
{
    getRoot().renameNodeParameter(old_name, new_name);
}


/**
 * Change node ids to be as in reference
 * @param const reference tree whose node ids will be used
 * @return tree with node ids as in reference but topology and branch lengths as in target
 */
void Tree::renumberNodes(const Tree &reference)
{
    // First, internal nodes
    std::map<RbBitSet, TopologyNode*> ref = reference.getBitsetToNodeMap();
    std::map<RbBitSet, TopologyNode*> toret = getBitsetToNodeMap();
    for (std::map<RbBitSet,TopologyNode*>::iterator it=ref.begin(); it!=ref.end(); ++it) {
      toret[it->first]->setIndex( it->second->getIndex());
    }
    // Second, the tip nodes
    std::vector<std::string> tipNames = getTipNames();
    for (size_t i = 0; i<tipNames.size(); ++i)
    {
      getTipNodeWithName(tipNames[i]).setIndex(reference.getTipNodeWithName(tipNames[i]).getIndex());
    }
    
}


void Tree::reroot(const Clade &o, bool make_bifurcating, bool reindex)
{
    bool strict = true;

    // for safety we reset the bitrepresentation of the clade
    Clade outgroup = o;
    outgroup.resetTaxonBitset( getTaxonBitSetMap() );

    if ( root->containsClade(outgroup, strict ) == false )
    {

        // check for the inverted clade
        RbBitSet b = outgroup.getBitRepresentation();
        b.flip();
        outgroup.setBitRepresentation(b);

        if ( root->containsClade(outgroup, strict ) == false )
        {
            throw RbException() << "Cannot reroot the tree because we could not find an outgroup with name '" << outgroup.toString() << "'.";
        }

    }

    // get the node representing the outgroup
    TopologyNode *outgroup_node = root->getNode(outgroup, strict);

    // check that we properly got a node
    if ( outgroup_node == NULL )
    {
        throw RbException() << "Cannot reroot the tree because we could not find an outgroup with name '" << outgroup.toString() << "'.";
    }

    // only reroot if the node is not the root node
    if ( outgroup_node->isRoot() == false )
    {
        // reset parent/child relationships
        reroot(*outgroup_node, make_bifurcating, reindex);
    }

}


void Tree::reroot(const std::string &outgroup, bool make_bifurcating, bool reindex)
{
    std::vector<std::string> tip_names = getTipNames();
    size_t outgroup_index = tip_names.size();
    for (size_t i=0; i<tip_names.size(); ++i)
    {
        if ( tip_names[i] == outgroup )
        {
            outgroup_index = i;
            break;
        }
    }

    if ( outgroup_index == tip_names.size() )
    {
        throw RbException() << "Cannot reroot the tree because we could not find an outgroup with name '" << outgroup << "'.";
    }
    
    TopologyNode& outgroup_node = getTipNode( outgroup_index );
    reroot(outgroup_node, make_bifurcating, reindex);

}


void Tree::reroot(TopologyNode &n, bool make_bifurcating, bool reindex)
{
    // reset parent/child relationships
    reverseParentChild( n.getParent() );
    n.getParent().setParent( NULL );
    
    // set the new root
    setRoot( &n.getParent(), reindex );
    
    // do we want to make the tree bifurcating?
    if ( make_bifurcating == true )
    {
        // first, if we have a trifurcation at the root, we need to change it into a bifurcation
        makeRootBifurcating(n.getClade());

        // second, we suppress any internal nodes of outdegree 1 ("knuckles" / sampled ancestors)
        suppressOutdegreeOneNodes(true);
        
    }
}


/**
 * Resets a map of the taxa to their BitSet indices.
 * The taxa are ordered alphabetically in the BitSet.
 */
void Tree::resetTaxonBitSetMap( void )
{
    taxon_bitset_map.clear();
    
    // get all taxon names
    std::vector<Taxon> unordered_taxa = getTaxa();
    std::vector<std::string> ordered_taxa;
    for (size_t i = 0; i < unordered_taxa.size(); ++i)
    {
        ordered_taxa.push_back(unordered_taxa[i].getName());
    }

    // order taxon names
    std::sort(ordered_taxa.begin(), ordered_taxa.end());

    // add taxa to bitset map
    for (size_t i = 0; i < ordered_taxa.size(); ++i)
    {
        taxon_bitset_map[ordered_taxa[i]] = i;
    }
    
}


void Tree::resolveMultifurcations( bool resolve_root )
{
    // Collect the multifurcation nodes.
    std::vector<TopologyNode*> multifurcating_nodes;
    for (auto node: nodes)
        if (node->getNumberOfChildren() > 2)
            multifurcating_nodes.push_back(node);

    // set the 'use_ages' attribute for all nodes
    root->setUseAges( isTimeTree(), true );

    // Resolve the multifurcations.
    for(auto node: multifurcating_nodes)
        node->resolveMultifurcation( resolve_root );

    // We need to reindex nodes because we added new ones
    reindexNodes();
}


TopologyNode& Tree::reverseParentChild(TopologyNode &n)
{
    // This routine makes n the new root by
    // * first making n.getParent() the new root, and then
    // * making n.getParent() into a child of n.
    // It returns the original root.

    TopologyNode* ret = &n;

    if ( !n.isRoot() )
    {
        TopologyNode &p = n.getParent();
        ret = &(reverseParentChild(p));

        // we need to re-orient the branches/indices so that
        // nodes remain associated with the same parameters
        p.setIndex(n.getIndex());
        p.setBranchLength(n.getBranchLength());

        p.removeChild( &n );
        n.setParent( nullptr ); // Avoid loop in parent edges from n -> p -> n
        p.setParent( &n );
        n.addChild( &p );
    }

    return *ret;
}


void Tree::setNegativeConstraint(bool tf)
{
    is_negative_constraint = tf;
}


void Tree::setRooted(bool tf)
{
    rooted = tf;
}


void Tree::setRoot( TopologyNode* r, bool reindex )
{

    // delete the old root if it's not in this tree
    bool found = false;

    TopologyNode* old_root = root;

    // set the root
    root = r;

    nodes.clear();

    // bootstrap all nodes from the root and add the in a pre-order traversal
    fillNodesByPhylogeneticTraversal(r);

    if ( reindex == true )
    {
        for (unsigned int i = 0; i < nodes.size(); ++i)
        {
            nodes[i]->setIndex(i);
        }
    }
    else
    {
        orderNodesByIndex();
    }

    num_nodes = nodes.size();

    // count the number of tips
    num_tips = 0;
    for (size_t i = 0; i < num_nodes; ++i)
    {
        if ( nodes[i] == NULL )
        {
            throw RbException("Problem while setting the root of the tree.");
        }
        if ( nodes[i] == old_root)
        {
            found = true;
        }
        num_tips += ( nodes[i]->isTip() ? 1 : 0);
    }


    root->setTree( this );

    if ( found == false )
    {
        delete old_root;
    }


}


//!< Set the indices of the taxa from the taxon map
void Tree::setTaxonIndices(const TaxonMap &tm)
{
    assert(root);

    // 1. Initialize the number of times that we have seen each map taxon to 0
    std::map<Taxon,int> use_count;
    for(int i=0;i<tm.size();i++)
        use_count.insert({tm.getTaxon(i), 0});

    // 2. Initialize the nodes to their index.
    int next_non_taxon_index = use_count.size();
    for(auto& node: nodes)
    {
        if (node->isTip() or node->getNumberOfChildren() == 1)
        {
            // 3a. If this is a tip or sampled ancester, check that we have an index for it.
            auto& taxon = node->getTaxon();
            if (not tm.hasTaxon(taxon))
                throw RbException()<<"setTaxonIndices: tree has unexpected "<<taxon<<" with no index";

            // 3b. Record that we used this taxon
            auto& count = use_count.at(taxon);
            count++;

            // 3c. Set the index for the node.
            node->setIndex( tm.getTaxonIndex(taxon) );
        }
        else
        {
            // 4. Otherwise set an index starting after all the taxon indices.
            node->setIndex( next_non_taxon_index++ );
        }
    }

    // 5. Check that all the taxon labels are used exactly once.
    for(auto& [taxon, count]: use_count)
    {
        if (count < 1)
        {
            for(auto& node: nodes)
            {
                if (node->getTaxon() == taxon)
                    throw RbException()<<"setTaxonIndices: taxon "<<taxon<<" with "<<node->getNumberOfChildren()<<" children is neither a tip nor a sampled ancestor";
            }
            throw RbException()<<"setTaxonIndices: tree is missing the taxon \'"<<taxon<<"\'";
        }
        if (count > 1)
            throw RbException()<<"setTaxonIndices: tree has multiple copies of the taxon "<<taxon;
    }

    // 6. Sort the nodes by their index.
    orderNodesByIndex();
}


/**
 * Change the name of a taxon
 *
 * \param[in] current_name    self explanatory.
 * \param[in] new_name         self explanatory.
 */
void Tree::setTaxonName(const std::string& current_name, const std::string& new_name)
{

    TopologyNode& node = getTipNodeWithName( current_name );
    Taxon& t = node.getTaxon();
    t.setName( new_name );
    
    // clear the taxon bitset map
    // the next time someone call getTaxonBitset() it will be rebuilt
    taxon_bitset_map.clear();
}


/**
 * Change the taxon object with specified name
 *
 * \param[in] current_name   self explanatory.
 * \param[in] new_taxon      self explanatory.
 */
void Tree::setTaxonObject(const std::string& current_name, const Taxon& new_taxon)
{

    const std::string &new_name = new_taxon.getName();

    TopologyNode& node = getTipNodeWithName( current_name );
    node.setTaxon( new_taxon );

    // clear the taxon bitset map
    // the next time someone call getTaxonBitset() it will be rebuilt
    taxon_bitset_map.clear();

}


void Tree::suppressOutdegreeOneNodes(bool replace)
{
    // if we are replacing rather than removing outdegree-1 nodes, we can let recursion take care of the
    // entire tree
    if (replace)
    {
        getRoot().suppressOutdegreeOneNodes( true );
    }
    // if we want them removed, we need to iterate over all nodes (potentially multiple times)
    else
    {
        for (size_t i = 0; i < nodes.size(); i++)
        {
            if (nodes[i]->isRoot() and nodes[i]->getNumberOfChildren() == 1)
            {
                TopologyNode* new_root = &nodes[i]->getChild(0);
                double new_bl = nodes[i]->getBranchLength() + new_root->getBranchLength();
                nodes[i]->removeChild( new_root );
                new_root->setParent( NULL );
                new_root->setBranchLength( new_bl );

                root = new_root;
                root->setTree(this);
            }
            else
            {
                // this will only have effect if nodes[i]->getNumberOfChildren() == 1
                nodes[i]->suppressOutdegreeOneNodes( false );
                // we need to reindex nodes because we removed an index
                reindexNodes();
            }
        }
            
        int outdegreeOneNodes = 0;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            if (nodes[i]->getNumberOfChildren() == 1)
            {
                outdegreeOneNodes++;
            }
        }
            
        // we have eliminated cases of -->A--> but not of -->A-->B-->; call myself recursively to handle
        // the latter as well
        if (outdegreeOneNodes > 0) suppressOutdegreeOneNodes( false );
    }

    // reindex here, in case suppressOutdegreeOneNodes
    // * replaced outdegree-1 nodes ("knuckles") by new tips with no index (replace = true)
    // * removed them entirely                                             (replace = false)
    
    // we need to reset the root so that the vector of nodes get filled again with the new number of nodes
    setRoot( &getRoot(), true );

    // clear the taxon bitset map
    // the next time someone call getTaxonBitset() it will be rebuilt
    taxon_bitset_map.clear();
}


void Tree::unroot( void )
{

    if ( isRooted() == true )
    {

        // get the root node because we need to make this tree unrooted (for topology comparison)
        TopologyNode *old_root = &getRoot();

        // make the tree use branch lengths instead of ages
        old_root->setUseAges(false, true);

        size_t child_index = 0;
        if ( old_root->getChild(child_index).isTip() == true )
        {
            child_index = 1;
        }
        TopologyNode *new_root = &old_root->getChild( child_index );
        TopologyNode *second_child = &old_root->getChild( (child_index == 0 ? 1 : 0) );

        double bl_first = new_root->getBranchLength();
        double bl_second = second_child->getBranchLength();

        old_root->removeChild( new_root );
        old_root->removeChild( second_child );
        new_root->setParent( NULL );
        new_root->addChild( second_child );
        second_child->setParent( new_root );

        second_child->setBranchLength( bl_first + bl_second );

        // finally we need to set the new root to our tree copy
        setRooted( false );
        setRoot( new_root, true);

    }

}


// Write this object into a file in its default format.
void Tree::writeToFile( const path &dir, const std::string &fn ) const
{
    
    // do not write a file if the tree is invalid
    if (this->getNumberOfTips() > 1)
    {
        path filename = dir / (fn + ".newick");
        create_directories(dir);

        // open the stream to the file
        std::ofstream outStream( filename.string() );

        // write the value of the node
        outStream << getNewickRepresentation();
        outStream << std::endl;

        // close the stream
        outStream.close();
    }
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const Tree& x)
{

    o << x.getNewickRepresentation();

    return o;
}
