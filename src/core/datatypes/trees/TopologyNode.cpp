#include <stdio.h>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <map>
#include <math.h>
#include <set>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <limits>

#include "Clade.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RbSettings.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeChangeEventMessage.h"
#include "RbBitSet.h"
#include "TaxonMap.h"
#include "TreeChangeEventHandler.h"
#include "RbConstants.h" // IWYU pragma: keep

using namespace RevBayesCore;

using boost::optional;

/** Default constructor (interior node, no name). Give the node an optional index ID */
TopologyNode::TopologyNode() {}


/** Default constructor (interior node, no name). Give the node an optional index ID */
TopologyNode::TopologyNode(size_t indx)
    : index(indx)
{

}


/** Constructor of node with name. Give the node an optional index ID */
TopologyNode::TopologyNode(const Taxon& t, const optional<size_t>& indx) :
    taxon(t),
    index(indx)
{

}


/** Constructor of node with name. Give the node an optional index ID */
TopologyNode::TopologyNode(const std::string& n, const optional<size_t>& indx) :
    taxon(n),
    index(indx)
{

}

/** Copy constructor. We use a shallow copy. */
TopologyNode::TopologyNode(const TopologyNode &n) :
    use_ages( n.use_ages ),
    age( n.age ),
    branch_length( n.branch_length ),
    parent( n.parent ),
    tree( NULL ),
    taxon( n.taxon ),
    index( n.index ),
    interior_node( n.interior_node ),
    root_node( n.root_node ),
    tip_node( n.tip_node ),
    sampled_ancestor( n.sampled_ancestor ),
    node_comments( n.node_comments ),
    branch_comments( n.branch_comments ),
    time_in_states( n.time_in_states ),
    num_shift_events( n.num_shift_events ),
    burst_speciation( n.burst_speciation ),
    sampling_event( n.sampling_event ),
    serial_sampling( n.serial_sampling ),
    serial_speciation( n.serial_speciation )
{

    // copy the children
    for (std::vector<TopologyNode*>::const_iterator it = n.children.begin(); it != n.children.end(); it++)
    {
        TopologyNode* the_node = *it;
        TopologyNode* theClone = the_node->clone();
        children.push_back( theClone );
        theClone->setParent(this);
    }


}


/** Destructor */
TopologyNode::~TopologyNode(void)
{
    // we do not own the parent so we do not delete it

    // free memory of children
    removeAllChildren();

    // make sure that I was removed from my parent
    if (parent != NULL)
    {
        parent->removeChild(this);
    }
    
}


TopologyNode& TopologyNode::operator=(const TopologyNode &n)
{

    if (this == &n)
    {

        removeAllChildren();

        // copy the members
        age                     = n.age;
        branch_comments         = n.branch_comments;
        branch_length           = n.branch_length;
        burst_speciation        = n.burst_speciation;
        index                   = n.index;
        interior_node           = n.interior_node;
        node_comments           = n.node_comments;
        parent                  = n.parent;
        root_node               = n.root_node;
        sampled_ancestor        = n.sampled_ancestor;
        sampling_event          = n.sampling_event;
        serial_sampling         = n.serial_sampling;
        serial_speciation       = n.serial_speciation;
        taxon                   = n.taxon;
        time_in_states          = n.time_in_states;
        num_shift_events        = n.num_shift_events;
        tip_node                = n.tip_node;
        use_ages                = n.use_ages;

        // copy the children
        for (std::vector<TopologyNode*>::const_iterator it = n.children.begin(); it != n.children.end(); it++)
        {
            children.push_back( (*it)->clone() );
        }

        // add myself as a new child to the parent node
        parent->addChild(this);

    }

    return *this;
}


void TopologyNode::addBranchParameter(const std::string &n, double p)
{

    if ( n == "index" || n == "species" )
    {
        std::cerr << "Illegal branch parameter with name '" << n << "'.\n";
    }

    std::stringstream o;
    char s[32];
    snprintf(s, sizeof(s), "%f",p);
    o << n << "=" << s;
    std::string comment = o.str();
    branch_comments.push_back( comment );

}


void TopologyNode::addBranchParameter(const std::string &n, const std::string &p)
{

    if ( n == "index" || n == "species" )
    {
        std::cerr << "Illegal branch parameter with name '" << n << "'.\n";
    }

    std::string comment = n + "=" + p;
    branch_comments.push_back( comment );

}



void TopologyNode::addBranchParameters(std::string const &n, const std::vector<double> &p, bool internalOnly) {

    if ( !internalOnly || !isTip()  )
    {
        std::stringstream o;
        o << n << "=" << p[ getIndex() ];
        std::string comment = o.str();
        branch_comments.push_back( comment );

        for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
        {
            (*it)->addBranchParameters(n, p, internalOnly);
        }

    }

}

void TopologyNode::addBranchParameters(std::string const &n, const std::vector<std::string> &p, bool internalOnly) {

    if ( !internalOnly || !isTip()  )
    {
        std::string comment = n + "=" + p[ getIndex() ];
        branch_comments.push_back( comment );

        for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
        {
            (*it)->addBranchParameters(n, p, internalOnly);
        }

    }

}


/** Add a child node. We own it from here on. */
void TopologyNode::addChild(TopologyNode* c, size_t pos )
{
    // add child to beginning if pos is out of bounds
    if( pos > children.size() )
    {
        throw RbException("Child position index out of bounds");
    }

    // add the child at pos offset from the end
    children.insert((children.rbegin() + pos).base(), c);

    // fire tree change event
    if ( tree != NULL )
    {
        tree->getTreeChangeEventHandler().fire( *c, RevBayesCore::TreeChangeEventMessage::TOPOLOGY );
    }

    tip_node = false;
    interior_node = true;

}


void TopologyNode::addNodeParameter(const std::string &n, double p)
{

    if ( n == "index" || n == "species" )
    {
        std::cerr << "Illegal node parameter with name '" << n << "' and value "<< p <<".\n";
    }

    std::stringstream o;
    char s[32];
    snprintf(s, sizeof(s), "%f",p);
    o << n << "=" << s; //SK
    std::string comment = o.str();
    node_comments.push_back( comment );

}


void TopologyNode::addNodeParameter(const std::string &n, const std::string &p)
{

    if ( n == "index" || n == "species" )
    {
        std::cerr << "Illegal node parameter with name '" << n << "' and value "<< p <<".\n";
    }

    std::string comment = n + "=" + p;
    node_comments.push_back( comment );

}


void TopologyNode::addNodeParameters(std::string const &n, const std::vector<double> &p, bool internalOnly)
{

    if ( internalOnly == false || isTip() == false  )
    {
        std::stringstream o;
        char s[32];
        size_t num_tip_nodes = 0;
        if ( internalOnly == true && tree != NULL )
        {
            num_tip_nodes = tree->getNumberOfTips();
        }
        snprintf(s, sizeof(s), "%f",p[ getIndex() - num_tip_nodes]);
        o << n << "=" << s; //SK
        std::string comment = o.str();
        node_comments.push_back( comment );

        for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
        {
            (*it)->addNodeParameters(n, p, internalOnly);
        }

    }

}

void TopologyNode::addNodeParameters(std::string const &n, const std::vector<std::string> &p, bool internal_only)
{

    if ( !internal_only || !isTip()  )
    {
        std::stringstream o;
        o << n << "=" << p[ getIndex() ];
        std::string comment = o.str();
        node_comments.push_back( comment );

        for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
        {
            (*it)->addNodeParameters(n, p, internal_only);
        }
    }
}


/*
 * Build newick string.
 * If simmap = true build a newick string compatible with SIMMAP and phytools.
 */
std::string TopologyNode::buildNewickString( bool simmap = false, bool round = true )
{

    // create the newick string
    std::stringstream o;

    std::fixed(o);
    // depending on the value of round, get standard precision or maximum
    if (round)
    {
        o.precision( 6 );
    }
    else
    {
        o.precision( std::numeric_limits<double>::digits10 );
    }

    std::vector<std::string> fossil_comments;

    // ensure we have an updated copy of branch_length variables
    if ( isRoot() == false )
    {
        recomputeBranchLength();
    }

    // test whether this is a internal or external node
    if ( tip_node == true )
    {
        // this is a tip so we just return the name of the node
        o << taxon.getName();

    }
    else
    {
        std::string fossil_name = "";

        o << "(";
        size_t j = 0;
        for (size_t i=0; i< children.size(); i++)
        {
            if ( RbSettings::userSettings().getCollapseSampledAncestors()
                    && children[i]->isSampledAncestor()
                    && (children[i]->getName() < fossil_name || fossil_name == "" ) )
            {
                fossil_name = children[i]->getName();
                fossil_comments = children[i]->getNodeParameters();
            }
            else
            {
                if (j > 0)
                {
                    o << ",";
                }
                j++;
                o << children[i]->buildNewickString( simmap, round );
            }
        }

        o << ")" << fossil_name;

    }

    if ( ( node_comments.size() + fossil_comments.size() > 0 || RbSettings::userSettings().getPrintNodeIndex() == true ) && simmap == false )
    {
        o << "[&";

        bool needsComma = false;

        // first let us print the node index, we must increment by 1 to match RevLanguage indexing
        if ( RbSettings::userSettings().getPrintNodeIndex() == true )
        {
            o << "index=" << getIndex() + 1;
            needsComma = true;
        }

        for (size_t i = 0; i < node_comments.size(); ++i)
        {
            if ( needsComma == true )
            {
                o << ",";
            }
            o << node_comments[i];
            needsComma = true;
        }

        for (size_t i = 0; i < fossil_comments.size(); ++i)
        {
            if ( needsComma == true )
            {
                o << ",";
            }
            o << fossil_comments[i];
            needsComma = true;
        }

        //Finally let's print the species name (always)
//        if ( needsComma == true )
//        {
//            o << ",";
//        }
//        o << "&species=" << getSpeciesName();

        o << "]";
    }

    if ( simmap == false )
    {
        double br = getBranchLength();

        if( RevBayesCore::RbMath::isNan(br) == false )
        {
            o << ":" << br;
        }
    }
    else
    {
        if ( isRoot() == false )
        {
            bool found = false;
            for (size_t i = 0; i < node_comments.size(); ++i)
            {
                if ( node_comments[i].substr(0, 18) == "character_history=" )
                {
                    o << ":" << node_comments[i].substr(18, node_comments[i].length());
                    found = true;
                    break;
                }
            }
            if ( found == false )
            {
                throw RbException("Error while writing SIMMAP newick string: no character history found for node.");
            }
        }
    }

    if ( branch_comments.size() > 0 && simmap == false )
    {
        o << "[&";
        for (size_t i = 0; i < branch_comments.size(); ++i)
        {
            if ( i > 0 )
            {
                o << ",";
            }
            o << branch_comments[i];
        }
        o << "]";
    }

    if ( root_node == true )
    {
        o << ";";
    }

    return o.str();
}



void TopologyNode::clearParameters(void)
{
    clearBranchParameters();
    clearNodeParameters();
}


void TopologyNode::clearBranchParameters( void )
{

    branch_comments.clear();
    if ( !isTip()  )
    {

        for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
        {
            (*it)->clearBranchParameters();
        }
    }
}


void TopologyNode::clearNodeParameters( void )
{

    node_comments.clear();
    if ( !isTip()  )
    {

        for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
        {
            (*it)->clearNodeParameters();
        }
    }
}


/** Clone function */
TopologyNode* TopologyNode::clone(void) const
{

    return new TopologyNode(*this);
}



std::string TopologyNode::computeNewick( bool round )
{

    return buildNewickString(false, round);
}


/* Build newick string */
std::string TopologyNode::computePlainNewick( void ) const
{

    // test whether this is a internal or external node
    if ( tip_node == true )
    {
        // this is a tip so we just return the name of the node
        return taxon.getName();
    }
    else
    {
        std::string fossil = "";
        std::string newick = "(";
        std::vector<std::string> child_newick;
        for (size_t i = 0; i < getNumberOfChildren(); ++i)
        {
            const TopologyNode& child = getChild( i );
            if ( RbSettings::userSettings().getCollapseSampledAncestors()
                    && child.isSampledAncestor()
                    && (child.getName() < fossil || fossil == "") )
            {
                fossil = child.getName();
            }
            else
            {
                child_newick.push_back( child.computePlainNewick() );
            }
        }
        sort(child_newick.begin(), child_newick.end());
        for (std::vector<std::string>::iterator it = child_newick.begin(); it != child_newick.end(); ++it)
        {
            if ( it != child_newick.begin() )
            {
                newick += ",";
            }
            newick += *it;
        }
        newick += ")";
        newick += fossil;

        return newick;
    }

}


std::string TopologyNode::computeSimmapNewick( bool round )
{
    return buildNewickString( true, round );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
bool TopologyNode::containsClade(const TopologyNode *c, bool strict) const
{
    RbBitSet your_taxa = RbBitSet( tree->getNumberOfTips() );
    c->getTaxa( your_taxa );

    return containsClade( your_taxa, strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
bool TopologyNode::containsClade(const Clade &c, bool strict) const
{

    return containsClade( c.getBitRepresentation(), strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
bool TopologyNode::containsClade(const RbBitSet &your_taxa, bool strict) const
{
    size_t n = tree->getNumberOfTips();
    RbBitSet my_taxa   = RbBitSet( n );
    getTaxa( my_taxa );

    if ( your_taxa.size() != my_taxa.size() )
    {
        throw RbException("Cannot check if the clade is contained within a node because of a problem in bit representation of clades.");
    }

    // this node needs to have at least as many taxa to contain the other clade
    if ( your_taxa.getNumberSetBits() > my_taxa.getNumberSetBits() )
    {
        // quick negative abort to safe computational time
        return false;
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.isSet(i) == true && my_taxa.isSet(i) == false )
        {
            return false;
        }

    }

    // now check, if required, that the contained clade is monophyletic in the containing clade.
    if ( strict == true )
    {
        // we already know from our check above that all taxa from the contained clade are present in this clade.
        // so we just need to check if there are additional taxa in this clade
        // and if so, then we need to check that the contained clade is contained in one of my children.
        if ( your_taxa.getNumberSetBits() < my_taxa.getNumberSetBits() )
        {

            // loop over all children
            for (std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
            {
                // check if the clade is contained in this child
                bool is_contained_in_child = (*it)->containsClade( your_taxa, strict );
                if ( is_contained_in_child == true )
                {
                    // yeah, so we can abort and return true
                    return true;
                }
            }

            return false;
        }

    }

    return true;
}



bool TopologyNode::doesUseAges( void ) const
{
    return use_ages;
}



bool TopologyNode::equals(const TopologyNode& node) const
{

    if (this == &node)
    {
        return true;
    }

    // test if the name is the same
    if (taxon != node.taxon)
    {
        return false;
    }

    // test if the index is the same
    if (index != node.index)
    {
        return false;
    }

    // test if the parent is the same
    if (parent != node.parent)
    {
        return false;
    }

    // test if the size of the children differs
    if (children.size() != node.children.size())
    {
        return false;
    }

    // test if all children are the same
    for (size_t i=0; i<children.size(); i++)
    {
        if (children[i]->equals(*node.children[i]))
        {
            return false;
        }
    }

    return true;
}


/*
 * Fill this map recursively with all clade indices.
 */
std::string TopologyNode::fillCladeIndices(std::map<std::string,size_t> &clade_index_map) const
{

    std::string newick = "";

    if ( isTip() == true )
    {
        newick = taxon.getName();
    }
    else
    {
        std::string fossil = "";
        newick = "(";
        std::vector<std::string> child_newick;
        for (size_t i = 0; i < getNumberOfChildren(); ++i)
        {
            const TopologyNode& child = getChild( i );
            if ( RbSettings::userSettings().getCollapseSampledAncestors()
                && child.isSampledAncestor()
                && (child.getName() < fossil || fossil == "") )
            {
                fossil = child.getName();
            }
            else
            {
                child_newick.push_back( child.fillCladeIndices(clade_index_map) );
            }
        }
        sort(child_newick.begin(), child_newick.end());
        for (std::vector<std::string>::iterator it = child_newick.begin(); it != child_newick.end(); ++it)
        {
            if ( it != child_newick.begin() )
            {
                newick += ",";
            }
            newick += *it;
        }
        newick += ")";
        newick += fossil;

    }

    // now insert the newick string for this node/clade with the index of this node
    clade_index_map.insert( std::pair<std::string,size_t>(newick, getIndex()) );


    // finally return my newick string so that my parents can use it
    return newick;
}


void TopologyNode::fireTreeChangeEvent( const unsigned& m )
{

    // fire tree change event
    if ( tree != NULL )
    {
        tree->getTreeChangeEventHandler().fire( *this, m );
    }

}

/*
 * Get the Age.
 * We internally store the age so can return it. However, if we invalidated the age ( age = Inf ),
 * then we need to compute the age from the time.
 */
double TopologyNode::getAge( void ) const
{

    return age;
}


RbBitSet TopologyNode::getAllClades(std::vector<RbBitSet> &all_clades, size_t num_tips, bool internal_only) const
{

    RbBitSet this_bs = RbBitSet(num_tips);
    if ( isTip() == true )
    {
//        taxa.set( index );
        // We can't use indices for tree comparison because different trees may
        // have different indices for the same taxon.
        // Instead make the BitSet ordered by taxon names.
        // Eventually this should be refactored with the TaxonMap class.
        std::map<std::string, size_t> taxon_bitset_map = tree->getTaxonBitSetMap();
        this_bs.set( taxon_bitset_map[taxon.getName()] );
        
        if ( internal_only == false )
        {
            all_clades.push_back( this_bs );
        }
    }
    else
    {
        for ( std::vector<TopologyNode* >::const_iterator i=children.begin(); i!=children.end(); i++ )
        {
            RbBitSet child_bs = (*i)->getAllClades(all_clades, num_tips, internal_only);
            this_bs |= child_bs;
        }
        all_clades.push_back( this_bs );

    }

    return this_bs;
}


/*
 * Get the branch length.
 * We compute the difference of my time and my parents time.
 */
double TopologyNode::getBranchLength( void ) const
{

    return branch_length;
}


/*
 * Get the branch parameters.
 */
const std::vector<std::string>& TopologyNode::getBranchParameters( void ) const
{

    return branch_comments;
}


/**
 * Get the index of a clade
 */
size_t TopologyNode::getCladeIndex(const TopologyNode *c) const
{

    size_t n = tree->getNumberOfTips();
    RbBitSet my_taxa   = RbBitSet( n );
    RbBitSet your_taxa = RbBitSet( n );
    getTaxa( my_taxa );
    c->getTaxa( your_taxa );

    // sanity check
    if ( your_taxa.size() != my_taxa.size() )
    {
        throw RbException("Cannot compute the clade index because of a problem in bit representation of clades.");
    }

    // this node needs to have at least as many taxa to contain the other clade
    if ( your_taxa.getNumberSetBits() > my_taxa.getNumberSetBits() )
    {
        // quick negative abort to safe computational time
        throw RbException("Node does not have at least as many taxa as input clade.");
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.isSet(i) == true && my_taxa.isSet(i) == false )
        {
            throw RbException("Node does not contain any taxa in clade.");
        }

    }

    // we already know from our check above that all taxa from the contained clade are present in this clade.
    // so we just need to check if there are additional taxa in this clade
    // and if so, then we need to check that the contained clade is contained in one of my children.
    if ( your_taxa.getNumberSetBits() < my_taxa.getNumberSetBits() )
    {

        // loop over all children
        for (std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
        {

            // check if the clade is contained in this child
            try
            {
                return (*it)->getCladeIndex( c );
            }
            catch(RbException&)
            {
                continue;
            }

        }

        // the clade is not one of my children, and we require strict identity
        throw RbException("Input clade is not a child node.");

    }

    // finally return my index
    return getIndex();
}


/** Get child at index i */
const TopologyNode& TopologyNode::getChild(size_t i) const
{

    return *children[i];
}


/** Get child at index i */
TopologyNode& TopologyNode::getChild(size_t i)
{

    return *children[i];
}


const std::vector<TopologyNode*>& TopologyNode::getChildren( void ) const
{
    return children;
}


/** Loop over children and get their indices */
std::vector<int> TopologyNode::getChildrenIndices() const
{

    std::vector<int> temp;

    for ( std::vector<TopologyNode* >::const_iterator i=children.begin(); i!=children.end(); i++ )
    {
        temp.push_back( int( (*i)->getIndex() ) );
    }

    return temp;
}


Clade TopologyNode::getClade( void ) const
{
    Clade c;

    // get the clade taxa
    std::vector<Taxon> taxa;

    if ( tree != NULL )
    {
        // initialize the clade bitset
        RbBitSet bitset( tree->getNumberOfTips() );
        getTaxa(taxa, bitset);

        c = Clade(taxa, bitset);
    }
    else
    {
        getTaxa(taxa);
        c = Clade(taxa);
    }

    c.setAge( age );

    std::set<Taxon> mrca;

    if( isTip() )
    {
        if( isSampledAncestor() )
        {
            mrca.insert( getTaxon() );
        }
    }
    else
    {
        // if a child is a sampled ancestor, its taxon is a mrca
        for (size_t i = 0; i < children.size(); i++)
        {
            if ( children[i]->isSampledAncestor() )
            {
                mrca.insert( children[i]->getTaxon() );
            }
        }
    }

    c.setMrca( mrca );

    return c;
}

bool TopologyNode::hasIndex( void) const
{
    return (bool)index;
}

size_t TopologyNode::getIndex( void ) const
{
    if (not index)
        throw RbException()<<"Problem while working with tree: Node index was never set.";

    return *index;
}

/**
 * Get the indices of nodes contained in the subtree starting with this node as the root.
 * This either returns 1 if this is a tip node (or 0 if we do not count tipes)
 * or computes recursively the number of nodes in both children plus one for this node.
 *
 * \param[in]   countTips   Shall we count tips?
 * \return                  Subtree size.
 */
void TopologyNode::getIndicesOfNodesInSubtree( bool countTips, std::vector<size_t>* indices ) const
{

    if ( tip_node )
    {
        if (countTips)
        {
            indices->push_back( getIndex() );
        }
    }
    else
    {
        indices->push_back( getIndex() );
        // now call this function recursively for all your children
        children[0]->getIndicesOfNodesInSubtree(countTips, indices);
        children[1]->getIndicesOfNodesInSubtree(countTips, indices);
    }

}


/**
 * Get the maximal depth starting from this node.
 * The depth here mean the maximal path length along the branches until a terminal node (tip) is reached.
 * For ultrametric trees all path lengths are equivalent, but for serial sampled trees not.
 * Hence, we compute the maximal depths by recursively exploring each path along the branches to the children.
 *
 * \return    The maximal depth (path length) from this node to the most recent tip.
 */
double TopologyNode::getMaxDepth( void ) const
{

    // iterate over the childen
    double max = 0.0;
    for (std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        double m = 0.0;
        TopologyNode& node = *(*it);
        if ( node.isTip() )
        {
            m = node.getBranchLength();
        }
        else
        {
            m = node.getBranchLength() + node.getMaxDepth();
        }

        if ( m > max )
        {
            max = m;
        }
    }

    return max;
}


const std::string& TopologyNode::getName( void ) const
{

    return getTaxon().getName();
}


/**
 * Is the argument clade contained in the clade descending from this node?
 */
TopologyNode* TopologyNode::getMrca(const Clade &c)
{

    return getNode( c, false );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 */
const TopologyNode* TopologyNode::getMrca(const Clade &c) const
{

    return getNode( c, false );
}

const TopologyNode* TopologyNode::getMrca(const Clade &c, bool strict) const
{

    return getNode( c, strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 */
const TopologyNode* TopologyNode::getMrca(const TopologyNode &n) const
{

    return getNode( n, false );
}

/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
TopologyNode* TopologyNode::getNode(const TopologyNode &n, bool strict)
{

    RbBitSet your_taxa = RbBitSet( tree->getNumberOfTips() );
    n.getTaxa( your_taxa );

    return getNode( your_taxa, strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
const TopologyNode* TopologyNode::getNode(const TopologyNode &n, bool strict) const
{

    RbBitSet your_taxa = RbBitSet( tree->getNumberOfTips() );
    n.getTaxa( your_taxa );

    return getNode( your_taxa, strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
TopologyNode* TopologyNode::getNode(const Clade &c, bool strict)
{

    return getNode( c.getBitRepresentation(), strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
TopologyNode* TopologyNode::getNode(const RbBitSet &your_taxa, bool strict)
{

    size_t n = tree->getNumberOfTips();
    RbBitSet my_taxa   = RbBitSet( n );
    getTaxa( my_taxa );

    if ( your_taxa.size() != my_taxa.size() )
    {
        throw RbException("Cannot retrieve a node because of a problem in bit representation of clades.");
    }

    // this node needs to have at least as many taxa to contain the other clade
    if ( your_taxa.getNumberSetBits() > my_taxa.getNumberSetBits() )
    {
        // quick negative abort to safe computational time
        return NULL;
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.isSet(i) == true && my_taxa.isSet(i) == false )
        {
            return NULL;
        }

    }


    // we already know from our check above that all taxa from the contained clade are present in this clade.
    // so we just need to check if there are additional taxa in this clade
    // and if so, then we need to check that the contained clade is contained in one of my children.
    if ( your_taxa.getNumberSetBits() < my_taxa.getNumberSetBits() )
    {

        // loop over all children
        for (std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
        {
            // check if the clade is contained in this child
            TopologyNode *is_contained_in_child = (*it)->getNode( your_taxa, strict );
            if ( is_contained_in_child != NULL )
            {
                // yeah, so we can abort and return true
                return is_contained_in_child;
            }
        }

        // now check, if required, that the contained clade is monophyletic in the containing clade.
        // this will only be done if we haven't found the clade within one of our children
        if ( strict == true )
        {
            return NULL;
        }

    }


    return this;
}



/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
const TopologyNode* TopologyNode::getNode(const Clade &c, bool strict) const
{

    return getNode( c.getBitRepresentation(), strict );
}


/**
 * Is the argument clade contained in the clade descending from this node?
 * By strict we mean that the contained clade has to be monophyletic in the containing clade.
 */
const TopologyNode* TopologyNode::getNode(const RbBitSet &your_taxa, bool strict) const
{
    size_t n = tree->getNumberOfTips();
    RbBitSet my_taxa   = RbBitSet( n );
    getTaxa( my_taxa );

    if ( your_taxa.size() != my_taxa.size() )
    {
        throw RbException("Cannot retrieve a (const) node because of a problem in bit representation of clades.");
    }

    // this node needs to have at least as many taxa to contain the other clade
    if ( your_taxa.getNumberSetBits() > my_taxa.getNumberSetBits() )
    {
        // quick negative abort to safe computational time
        return NULL;
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.isSet(i) == true && my_taxa.isSet(i) == false )
        {
            return NULL;
        }

    }

    // we already know from our check above that all taxa from the contained clade are present in this clade.
    // so we just need to check if there are additional taxa in this clade
    // and if so, then we need to check that the contained clade is contained in one of my children.
    if ( your_taxa.getNumberSetBits() < my_taxa.getNumberSetBits() )
    {

        // loop over all children
        for (std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
        {
            // check if the clade is contained in this child
            TopologyNode *is_contained_in_child = (*it)->getNode( your_taxa, strict );
            if ( is_contained_in_child != NULL )
            {
                // yeah, so we can abort and return true
                return is_contained_in_child;
            }
        }

        // now check, if required, that the contained clade is monophyletic in the containing clade.
        // this will only be done if we haven't found the clade within one of our children
        if ( strict == true )
        {
            return NULL;
        }

    }

    return this;
}


/*
 * Get the node parameters.
 */
const std::vector<std::string>& TopologyNode::getNodeParameters( void ) const
{

    return node_comments;
}



size_t TopologyNode::getNumberOfChildren( void ) const
{

    return children.size();
}


/**
 * Get the number of nodes contained in the subtree starting with this node as the root.
 * This either returns 1 if this is a tip node (or 0 if we do not count tipes)
 * or computes recursively the number of nodes in both children plus one for this node.
 *
 * \param[in]   countTips   Shall we count tips?
 * \return                  Subtree size.
 */
size_t TopologyNode::getNumberOfNodesInSubtree( bool countTips ) const
{

    if ( tip_node )
    {
        return (countTips ? 1 : 0);
    }
    else
    {
        return children[0]->getNumberOfNodesInSubtree(countTips) + children[1]->getNumberOfNodesInSubtree(countTips) + 1;
    }

}

TopologyNode& TopologyNode::getParent(void)
{

    return *parent;

}

const TopologyNode& TopologyNode::getParent(void) const
{

    return *parent;
}


std::string TopologyNode::getIndividualName() const
{
    std::string name = taxon.getSpeciesName();
    return name;
}


size_t TopologyNode::getNumberOfShiftEvents( void ) const
{
    return num_shift_events;
}


std::string TopologyNode::getSpeciesName() const
{
    std::string name = taxon.getSpeciesName();
    return name;
}


void TopologyNode::getTaxa(std::vector<Taxon> &taxa) const
{

    if ( isTip() )
    {
        taxa.push_back( taxon );
    }
    else
    {
        for ( std::vector<TopologyNode* >::const_iterator i=children.begin(); i!=children.end(); i++ )
        {
            (*i)->getTaxa( taxa );
        }
    }


}


void TopologyNode::getTaxa(RbBitSet &taxa) const
{

    if ( isTip() == true )
    {
//        taxa.set( index );
        // We can't use indices for tree comparison because different trees may
        // have different indices for the same taxon.
        // Instead make the BitSet ordered by taxon names.
        // Eventually this should be refactored with the TaxonMap class.
        std::map<std::string, size_t> taxon_bitset_map = tree->getTaxonBitSetMap();
        taxa.set( taxon_bitset_map[taxon.getName()] );
    }
    else
    {
        for ( std::vector<TopologyNode* >::const_iterator i=children.begin(); i!=children.end(); i++ )
        {
            (*i)->getTaxa( taxa );
        }
    }


}


void TopologyNode::getTaxa(std::vector<Taxon> &taxa, RbBitSet &bitset) const
{

    if ( isTip() )
    {
        taxa.push_back( taxon );
        std::map<std::string, size_t> taxon_bitset_map = tree->getTaxonBitSetMap();
        bitset.set( taxon_bitset_map[taxon.getName()] );
    }
    else
    {
        for ( std::vector<TopologyNode* >::const_iterator i=children.begin(); i!=children.end(); i++ )
        {
            (*i)->getTaxa( taxa, bitset );
        }
    }


}


const Taxon& TopologyNode::getTaxon( void ) const
{
    return taxon;
}


Taxon& TopologyNode::getTaxon( void )
{
    return taxon;
}


std::vector<double> TopologyNode::getTimeInStates()
{
    return time_in_states;
}


double TopologyNode::getTmrca(const Clade &c) const
{
    const std::vector<Taxon>& yourTaxa = c.getTaxa();

    return getTmrca( yourTaxa );
}


double TopologyNode::getTmrca(const TopologyNode &n) const
{
    std::vector<Taxon> yourTaxa;
    n.getTaxa( yourTaxa );

    return getTmrca( yourTaxa );
}

double TopologyNode::getTmrca(const std::vector<Taxon> &yourTaxa) const
{

    std::vector<Taxon> myTaxa;
    getTaxa( myTaxa );

    if ( myTaxa.size() < yourTaxa.size() )
    {
        return -1;
    }

    for (std::vector<Taxon>::const_iterator y_it = yourTaxa.begin(); y_it != yourTaxa.end(); ++y_it)
    {
        bool found = false;
        for (std::vector<Taxon>::const_iterator it = myTaxa.begin(); it != myTaxa.end(); ++it)
        {
            if ( *y_it == *it )
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            return -1;
        }
    }

    if ( myTaxa.size() == yourTaxa.size() )
    {
        return getAge();
    }
    else
    {
        double tmrca = getAge();
        bool contains = false;
        for (std::vector<TopologyNode*>::const_iterator it = children.begin(); it != children.end(); ++it)
        {
            double child_tmrca = (*it)->getTmrca( yourTaxa );
            contains |= ( child_tmrca >= 0.0 );
            if ( contains == true )
            {
                tmrca = child_tmrca;
                break;
            }
        }
        return tmrca;
    }
}


bool TopologyNode::isFossil( void ) const
{

    return age > 0.0 && tip_node;
}


bool TopologyNode::isInternal( void ) const
{

    return interior_node;
}


bool TopologyNode::isRoot( void ) const
{

    return parent == NULL;
}


bool TopologyNode::isSampledAncestor(  bool propagate ) const
{

    bool sa = sampled_ancestor;
    if( propagate == true )
    {
        for(size_t i = 0; i < children.size(); i++)
        {
            sa = sa || children[i]->isSampledAncestor(false);
        }
    }

    return sa;
}


bool TopologyNode::isTip( void ) const
{

    return tip_node;
}


bool TopologyNode::isUltrametric(double &depth) const
{
    // if we are simply a tip node, then this subtree must be ultrametric
    if ( isTip() == true )
    {
        depth = 0;
        return true;
    }
    
    double my_depth = 0.0;
    bool am_ultrametric = children[0]->isUltrametric(my_depth);
    my_depth += children[0]->getBranchLength();
    
    for ( size_t i=1; i<children.size(); ++i )
    {
        double child_depth = 0.0;
        bool child_is_ultrametric = children[i]->isUltrametric(child_depth);
        child_depth += children[i]->getBranchLength();
        
        am_ultrametric = (am_ultrametric && child_is_ultrametric);
        am_ultrametric = (am_ultrametric && fabs(child_depth - my_depth) < 1E-4);
    }
    
    return am_ultrametric;
}


/**
 * Make this node an all its children bifurcating.
 * The root will not be changed. We throw an error if this node
 * has more than 2 children. If this node has only one child,
 * then we insert a dummy child.
 * This function is called recursively.
 */
void TopologyNode::makeBifurcating( bool as_fossils )
{

    if ( isTip() == false )
    {

        // we need to be able to bifurcate sampled ancestor root nodes
        //if ( isRoot() == false )
        //{

            if ( getNumberOfChildren() == 1 )
            {

                // should we remove the node or add it as a fossil?
                if ( as_fossils == true )
                {
                    TopologyNode *new_fossil = new TopologyNode( getTaxon() );
                    taxon = Taxon("");

                    // connect to the old fossil
                    addChild( new_fossil );
                    new_fossil->setParent( this );

                    // set the fossil flags
                    setSampledAncestor( false );
                    new_fossil->setSampledAncestor( true );

                    // set the age and branch-length of the fossil
                    new_fossil->setAge( age );
                    new_fossil->setBranchLength( 0.0 );
                }
                else
                {
                    // we are going to delete myself by connect my parent and my child
                    TopologyNode& parent = getParent();
                    TopologyNode& child = getChild(0);
                    
                    // the new branch length needs to be the branch length of the parent and child
                    double summ = getBranchLength() + child.getBranchLength();
                    
                    // now remove myself from the parent
                    parent.removeChild( this );
                    
                    // and my child from me
                    removeChild( &child );
                    
                    // and stich my parent and my child together
                    parent.addChild( &child );
                    child.setParent( &parent );
                    
                    // finally, adapt the branch lengths
                    child.setBranchLength(summ);
                    
                }


            }

        //}

        // call this function recursively for all its children
        for (size_t i=0; i<getNumberOfChildren(); ++i)
        {
            getChild( i ).makeBifurcating( as_fossils );
        }

    }

}


void TopologyNode::recomputeAge( bool recursive )
{

    if ( children.size() == 0 )
    {
        age = 0.0;
    }
    else 
    {
        if ( recursive == true )
        {
            for (size_t i=0; i<children.size(); ++i)
            {
                children[i]->recomputeAge(recursive);
            }
        }
        age = children[0]->getBranchLength() + children[0]->getAge();
    }

}


void TopologyNode::recomputeBranchLength( void )
{

    if ( parent == NULL )
    {
        branch_length = 0.0;
    }
    else if ( RbMath::isFinite( age ) == true )
    {
        branch_length = parent->getAge() - age;
    }

}


/** Remove all children. We need to call intelligently the destructor here. */
void TopologyNode::removeAllChildren(void)
{

    // empty the children vector
    while (children.size() > 0)
    {
        TopologyNode* the_node = children[0];
        // free the memory
        delete the_node;
    }

    taxon = Taxon("");

    tip_node = true;
    interior_node = false;
}




/** Remove a child from the vector of children */
size_t TopologyNode::removeChild(TopologyNode* c)
{

    std::vector<TopologyNode* >::iterator it = find(children.begin(), children.end(), c);
    size_t pos = 0;
    if ( it != children.end() )
    {
        // get offset from the end
        pos = std::distance(it, children.end())-1;
        children.erase(it);
    }
    else
    {
        throw RbException("Cannot find node in list of children nodes");
    }

    // update the flags
    tip_node      = (children.size() == 0);
    interior_node = (children.size()  > 0);

    // fire tree change event
    if ( tree != NULL )
    {
        tree->getTreeChangeEventHandler().fire( *c, RevBayesCore::TreeChangeEventMessage::TOPOLOGY );
        tree->getTreeChangeEventHandler().fire( *this, RevBayesCore::TreeChangeEventMessage::TOPOLOGY );
    }

    return pos;
}


void TopologyNode::removeTree(Tree *t)
{

    // only remove the tree if we had a pointer stored to it
    if ( tree == t )
    {
        tree = NULL;
    }

    for (std::vector<TopologyNode *>::iterator i = children.begin(); i != children.end(); ++i)
    {
        (*i)->removeTree( t );
    }

}


void TopologyNode::renameNodeParameter(const std::string &old_name, const std::string &new_name)
{
    for (size_t i = 0; i < node_comments.size(); i++)
    {
        size_t equal_sign = node_comments[i].find("=");
        std::string param_name = node_comments[i].substr(0, equal_sign);
        if (param_name.compare(old_name) == 0)
        {
            node_comments[i] = new_name + node_comments[i].substr(equal_sign);
            break;
        }
    }

    for (std::vector<TopologyNode*>::iterator it = children.begin(); it != children.end(); ++it)
    {
        (*it)->renameNodeParameter(old_name, new_name);
    }
}


void TopologyNode::setAge(double a, bool propagate)
{
    if ( sampled_ancestor == true && propagate == true )
    {
        parent->setAge(a);
        return;
    }

    age = a;
    
//    // we should also update the taxon age if this is a tip node
//    if ( isTip() == true )
//    {
//        getTaxon().setAge( a );
//    }

    // we need to recompute my branch-length
    recomputeBranchLength();

    // we also need to recompute the branch lengths of my children
    for (std::vector<TopologyNode *>::iterator it = children.begin(); it != children.end(); ++it)
    {
        TopologyNode *child = *it;
        if ( child->isSampledAncestor() )
        {
            child->setAge(a, false);
        }
        else
        {
            child->recomputeBranchLength();
        }

        // fire tree change event
        if ( tree != NULL )
        {
            tree->getTreeChangeEventHandler().fire( *child, RevBayesCore::TreeChangeEventMessage::BRANCH_LENGTH );
        }
    }

}


void TopologyNode::setBranchLength(double b, bool flag_dirty)
{

    branch_length = b;
    
    // we need to mark that we are not using ages but branch lengths
    // otherwise, this code would need to reformat the branch length into ages!!!
    // Sebastian (20200819): For some forgotten reason this has been commented out.
    // I think there are issues if not all ages or branch lengths are set yet.
    // So the caller needs to make sure of what type of trees are used (unfortunately).
//    setUseAges( false, false);
    
    
    // fire tree change event
    if ( flag_dirty == true && tree != NULL )
    {
        tree->getTreeChangeEventHandler().fire( *this, RevBayesCore::TreeChangeEventMessage::BRANCH_LENGTH );
    }

}


void TopologyNode::setIndex( size_t idx )
{

    index = idx;

}


void TopologyNode::setName(std::string const &n)
{

    taxon.setName( n );
    taxon.setSpeciesName( n );

}


void TopologyNode::setNumberOfShiftEvents(size_t n)
{
    num_shift_events = n;
}

//SK
void TopologyNode::setNodeType(bool tip, bool root, bool interior)
{

	tip_node = tip;
	root_node = root;
	interior_node = interior;

}


void TopologyNode::setParent(TopologyNode* p)
{

    // we only do something if this isn't already our parent
    if (p != parent)
    {
        // we do not own the parent so we do not have to delete it
        parent = p;

        // we need to recompute our branch length
        recomputeBranchLength();

        // fire tree change event
        if ( tree != NULL )
        {
            tree->getTreeChangeEventHandler().fire( *this, RevBayesCore::TreeChangeEventMessage::DEFAULT );
        }

    }

    root_node = (parent == NULL);
}

void TopologyNode::setUseAges(bool tf, bool recursive)
{

    // if this node did use ages before but not we do not anymore
    if ( use_ages == true && tf == false )
    {

        // check if we need to call the recursion
        if ( recursive == true )
        {

            // call all our children
            for (size_t i=0; i<children.size(); ++i)
            {
                children[i]->setUseAges(tf, recursive);
            }

        }

        // now we need to compute the branch lengths
        recomputeBranchLength();

        // make the age not usable to be safe
        age = RbConstants::Double::nan;

    }
    // finally set our internal flag
    use_ages = tf;

}


void TopologyNode::setSampledAncestor(bool tf)
{

    sampled_ancestor = tf;

}


void TopologyNode::setSpeciesName(std::string const &n)
{

    taxon.setSpeciesName( n );

}


void TopologyNode::setTaxon(Taxon const &t)
{

    taxon = t;

}


//!< Set the indices of the taxa from the taxon map
void TopologyNode::setTaxonIndices(const TaxonMap &tm)
{

    if ( isTip() == true )
    {
        size_t idx = tm.getTaxonIndex( taxon );
        index = idx;
    }
    else
    {
        for ( std::vector<TopologyNode* >::const_iterator i=children.begin(); i!=children.end(); i++ )
        {
            (*i)->setTaxonIndices( tm );
        }
    }


}


void TopologyNode::setTimeInStates(std::vector<double> t)
{
    time_in_states = t;
}


void TopologyNode::setTree(Tree *t)
{

    tree = t;
    for (std::vector<TopologyNode *>::iterator i = children.begin(); i != children.end(); ++i)
    {
        (*i)->setTree( t );
    }

}
