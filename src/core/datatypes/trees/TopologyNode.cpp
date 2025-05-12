#include <cstdio>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <map>
#include <cmath>
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
#include "StringUtilities.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"

using namespace RevBayesCore;

using boost::optional;

using std::string;
using std::vector;

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
    sampled_ancestor_tip( n.sampled_ancestor_tip ),
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
        node_comments           = n.node_comments;
        parent                  = n.parent;
        sampled_ancestor_tip    = n.sampled_ancestor_tip;
        sampling_event          = n.sampling_event;
        serial_sampling         = n.serial_sampling;
        serial_speciation       = n.serial_speciation;
        taxon                   = n.taxon;
        time_in_states          = n.time_in_states;
        num_shift_events        = n.num_shift_events;
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
        throw RbException() << "Illegal branch parameter with name '" << n << "'.";
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
        throw RbException() << "Illegal branch parameter with name '" << n << "'.\n";
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


bool TopologyNode::hasNodeComment(const std::string& comment) const
{
    for(auto& node_comment: node_comments)
        if (node_comment == comment)
            return true;
    return false;
}

std::optional<std::string> TopologyNode::getNodeParameter(const std::string& name) const
{
    for(auto& node_comment: node_comments)
    {
        if (node_comment.substr(0,name.size()) == name and node_comment.size() > name.size() and node_comment[name.size()] == '=')
        {
            return node_comment.substr(name.size()+1);
        }
    }

    // Not found
    return {};
}

// If the parameter is already set, modify the existing comment instead of appending a new one.
bool TopologyNode::setNodeParameter(const std::string& name, const std::string& value)
{
    // value probably better not have any commas in it.

    for(auto& node_comment: node_comments)
    {
        // If this command starts with '<name>=', then set the value here.
        if (node_comment.substr(0,name.size()) == name and node_comment.size() > name.size() and node_comment[name.size()] == '=')
        {
            node_comment = name + "=" + value;
            return true;
        }
    }

    // Otherwise append a new comment to the end.
    node_comments.push_back(name + "=" + value);
    return false;
}

std::optional<std::string> TopologyNode::eraseNodeParameter(const std::string& name)
{
    // 1. Find the index of the parameter, if it exists.
    std::optional<int> found_index;
    for(int i=0;i<node_comments.size();i++)
    {
        auto& node_comment = node_comments[i];
        if (node_comment.substr(0,name.size()) == name and node_comment.size() > name.size() and node_comment[name.size()] == '=')
        {
            found_index = i;
            break;
        }
    }

    // 2. If it doesn't exist, then we're done.
    if (not found_index) return {};

    // 3. Save the parameter value.
    string value = node_comments[*found_index].substr(name.size()+1);

    // 4. Move the comment to the end of the array and pop it.
    if (*found_index < node_comments.size()-1)
    {
        std::swap(node_comments[*found_index], node_comments.back());
    }
    node_comments.pop_back();

    // 5. Return the saved value.
    return value;
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
        tree->getTreeChangeEventHandler().fire( *this, RevBayesCore::TreeChangeEventMessage::TOPOLOGY );
    }
}


void TopologyNode::addNodeParameter(const std::string &n, double p)
{

    if ( n == "index" || n == "species" )
    {
        throw RbException() << "Illegal node parameter with name '" << n << "' and value "<< p <<".\n";
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
        throw RbException() << "Illegal node parameter with name '" << n << "' and value "<< p <<".\n";
    }

    addNodeParameter_(n,p);
}


void TopologyNode::addNodeParameter_(const std::string &n, const std::string &p)
{
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


std::ostream& TopologyNode::buildNewick( std::ostream& o, bool simmap = false)
{
    // ensure we have an updated copy of branch_length variables
    if ( not isRoot() )
    {
        recomputeBranchLength();
    }

    // 1. Write out the child newicks if there are any.
    if (not isTip())
    {
        o << "(";
        for (size_t i=0; i< children.size(); i++)
        {
            if (i > 0)
            {
                o << ",";
            }
            children[i]->buildNewick(o, simmap);
        }
        o << ")";
    }

    // 2. Write out the node name is there is any.
    if (children.size() < 2)
        o << taxon.getName();

    // 3. Write out node comments if there are any and (simmap == false)
    if ( ( node_comments.size() > 0 or RbSettings::userSettings().getPrintNodeIndex() == true ) && simmap == false )
    {
        o << "[&";

        // first let us print the node index, we must increment by 1 to match RevLanguage indexing
        if ( RbSettings::userSettings().getPrintNodeIndex() == true )
        {
            o << "index=" << getIndex() + 1;
            if (node_comments.size() > 0)
                o<<",";
        }

        StringUtilities::join(o, node_comments, ",");

        o << "]";
    }

    // 4a. Write ":" + branch length + branch_comments if (simmap == false)
    if ( simmap == false )
    {
        double br = getBranchLength();

        if( RevBayesCore::RbMath::isNan(br) == false )
        {
            o << ":" << br;
        }

        if ( branch_comments.size() > 0 )
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
    }
    // 4b. Write ":" + simmap comment if (simmap == true)
    else
    {
        if ( isRoot() == false )
        {
            // Find the simmap comment
            optional<string> simmap_comment;
            for (auto& node_comment: node_comments)
            {
                if ( node_comment.substr(0, 18) == "character_history=" )
                {
                    simmap_comment = node_comment.substr(18);
                    break;
                }
            }

            // Print the simmap comment
            if ( simmap_comment )
                o << ":" << simmap_comment.value();
            else
                throw RbException("Error while writing SIMMAP newick string: no character history found for node.");
        }
    }

    // 5. Write ";" if we're done with the tree.
    if ( isRoot() )
    {
        // FIXME - move to caller?
        o << ";";
    }

    return o;
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

    buildNewick(o, simmap);

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
    /* NOTE: Representing a topology as with WITH NO ANNOTATIONS
     * means that we have to represent sampled ancestors as
     * outdegree-1 nodes.
     *
     * If you want to build a tree object from that, you may
     * need to call tree->suppressOutdegreeOneNodes(true, true).
     */

    // test whether this is a internal or external node
    if ( isTip() )
    {
        // this is a tip so we just return the name of the node
        return taxon.getName();
    }
    else
    {
	// If this is non-empty, then there is a taxon here.  Right?
        string node_name = taxon.getName();

        std::vector<std::string> child_newicks;
        for (size_t i = 0; i < getNumberOfChildren(); ++i)
	{
	    if (getChild(i).isSampledAncestorTip())
	    {
		assert(node_name == "");
		node_name = getChild(i).taxon.getName();
	    }
	    else
		child_newicks.push_back( getChild(i).computePlainNewick() );
	}

        sort(child_newicks.begin(), child_newicks.end());

        return "(" + StringUtilities::join(child_newicks, ",") + ")" + node_name;
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
    if ( your_taxa.count() > my_taxa.count() )
    {
        // quick negative abort to safe computational time
        return false;
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.test(i) == true && my_taxa.test(i) == false )
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
        if ( your_taxa.count() < my_taxa.count() )
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
        std::vector<std::string> child_newicks;
        for (size_t i = 0; i < getNumberOfChildren(); ++i)
            child_newicks.push_back( getChild(i).fillCladeIndices(clade_index_map) );

        sort(child_newicks.begin(), child_newicks.end());

        newick = "(" + StringUtilities::join(child_newicks, ",") + ")";
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
        int bit_index = tree->getTaxonBitSetMap().at(taxon.getName());

        this_bs.set( bit_index );
        
        if ( internal_only == false )
        {
            all_clades.push_back( this_bs );
        }
    }
    else
    {
        for ( auto& child: children )
            this_bs |= child->getAllClades(all_clades, num_tips, internal_only);
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
    if ( your_taxa.count() > my_taxa.count() )
    {
        // quick negative abort to safe computational time
        throw RbException("Node does not have at least as many taxa as input clade.");
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.test(i) == true && my_taxa.test(i) == false )
        {
            throw RbException("Node does not contain any taxa in clade.");
        }

    }

    // we already know from our check above that all taxa from the contained clade are present in this clade.
    // so we just need to check if there are additional taxa in this clade
    // and if so, then we need to check that the contained clade is contained in one of my children.
    if ( your_taxa.count() < my_taxa.count() )
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
    // Maybe this should be an assert, but not sure what the contract is here.
    if (i >= children.size())
        throw RbException()<<"Trying to get child "<<i<<", but node only has "<<children.size()<<" children.";

    return *children[i];
}


/** Get child at index i */
TopologyNode& TopologyNode::getChild(size_t i)
{
    // Maybe this should be an assert, but not sure what the contract is here.
    if (i >= children.size())
        throw RbException()<<"Trying to get child "<<i<<", but node only has "<<children.size()<<" children.";

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
        if( isSampledAncestorTip() )
        {
            mrca.insert( getTaxon() );
        }
    }
    else
    {
        // if a child is a sampled ancestor, its taxon is a mrca
        for (size_t i = 0; i < children.size(); i++)
        {
            if ( children[i]->isSampledAncestorTip() )
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
 * This either returns 1 if this is a tip node (or 0 if we do not count tips)
 * or computes recursively the number of nodes in both children plus one for this node.
 *
 * \param[in]   countTips   Shall we count tips?
 * \return                  Subtree size.
 */
void TopologyNode::getIndicesOfNodesInSubtree( bool countTips, std::vector<size_t>* indices ) const
{

    if ( isTip() )
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
    if ( your_taxa.count() > my_taxa.count() )
    {
        // quick negative abort to safe computational time
        return NULL;
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.test(i) == true && my_taxa.test(i) == false )
        {
            return NULL;
        }

    }


    // we already know from our check above that all taxa from the contained clade are present in this clade.
    // so we just need to check if there are additional taxa in this clade
    // and if so, then we need to check that the contained clade is contained in one of my children.
    if ( your_taxa.count() < my_taxa.count() )
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
    if ( your_taxa.count() > my_taxa.count() )
    {
        // quick negative abort to safe computational time
        return NULL;
    }

    // check that every taxon of the clade is in this subtree
    for (size_t i=0; i<n; ++i)
    {

        // if I don't have any of your taxa then I cannot contain you.
        if ( your_taxa.test(i) == true && my_taxa.test(i) == false )
        {
            return NULL;
        }

    }

    // we already know from our check above that all taxa from the contained clade are present in this clade.
    // so we just need to check if there are additional taxa in this clade
    // and if so, then we need to check that the contained clade is contained in one of my children.
    if ( your_taxa.count() < my_taxa.count() )
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

size_t TopologyNode::getDegree( void ) const
{
    int degree = children.size();
    if (not isRoot())
        degree++;
    return degree;
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

    if ( isTip() )
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
        int bit_index = tree->getTaxonBitSetMap().at(taxon.getName());

        taxa.set( bit_index );
    }
    else
    {
        for ( auto& child: children )
        {
            child->getTaxa( taxa );
        }
    }


}


void TopologyNode::getTaxa(std::vector<Taxon> &taxa, RbBitSet &bitset) const
{

    if ( isTip() )
    {
        taxa.push_back( taxon );
        int bit_index = tree->getTaxonBitSetMap().at(taxon.getName());
        bitset.set( bit_index );
    }
    else
    {
        for ( auto& child : children)
        {
            child->getTaxa( taxa, bitset );
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

    return age > 0.0 && isTip();
}


bool TopologyNode::isInternal( void ) const
{

    return not isTip();
}


bool TopologyNode::isRoot( void ) const
{

    return parent == NULL;
}


bool TopologyNode::isSampledAncestorTip() const
{
    // Only tips can have the sampled_ancestor_tip flag set.
    assert(not sampled_ancestor_tip or isTip());

    return sampled_ancestor_tip;
}


bool TopologyNode::isSampledAncestorParent() const
{
    // Only tips can have the sampled_ancestor_tip flag set.
    assert(not sampled_ancestor_tip or isTip());

    for(auto child: children)
	if (child->isSampledAncestorTip())
	    return true;

    return false;
}


bool TopologyNode::isSampledAncestorTipOrParent() const
{
    return isSampledAncestorTip() or isSampledAncestorParent();
}


// This function exists partly to distinguish "real" sampled ancestors
// from the parent/tip "fake" sampled ancestors.
bool TopologyNode::isSampledAncestorKnuckle() const
{
    // Only tips can have the sampled_ancestor_tip flag set.
    assert(not sampled_ancestor_tip or isTip());

    // This function intentionally does NOT query the sampled_ancestor_tip flag.
    // One question is whether we should require there to be a Taxon present here to be
    //   considered a sampled ancestor.
    
    return getNumberOfChildren() == 1;
}


bool TopologyNode::isTip( void ) const
{

    return children.empty();
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

    // fire tree change event
    if ( tree != NULL )
    {
        // tree->getTreeChangeEventHandler().fire( *c, RevBayesCore::TreeChangeEventMessage::TOPOLOGY );
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


void TopologyNode::resolveMultifurcation(bool resolve_root)
{
    if (children.size() < 3) return;
    
    // What if the root had 4 children?
    if (isRoot() and not resolve_root) return;

    std::cerr<<"\n\nresolveMultifurcation:  children.size() = "<<children.size()<<"\n";
    RandomNumberGenerator* rng = GLOBAL_RNG;
    // "active" children are those that are younger than the child currently under consideration
    std::vector<TopologyNode*> active_children;
            
    if (use_ages)
    {
        // The following is adapted from UniformSerialSampledTimeTreeDistribution::simulateCoalescentAges()
        std::vector<double> coalescence_times;

        std::vector<TopologyNode*> sorted_children = children;
        std::sort(sorted_children.begin(), sorted_children.end(), [](auto x, auto y) {return x->getAge() < y->getAge();});

        // for each tip, simulate an age between max age and tip age
        double max_age = getAge();

        // Skip the youngest child, since it could pick an age that is older than any sibling to coalesce with.
        // Skip the oldest child, since we only want n-2 coalescence events.
        for(int i = 1; i < sorted_children.size()-1; ++i)
        {
            // get the age of the tip
            double a = sorted_children[i]->getAge();

            // simulate the age of a node
            double new_age = a + rng->uniform01() * (max_age - a);

            // add the age to the vector of coalescence_times
            coalescence_times.push_back(new_age);
        }

        // sort the coalescence_times (from youngest to oldest)
        std::sort(coalescence_times.begin(), coalescence_times.end());

        active_children.push_back(sorted_children[0]);
        std::cerr<<"    child ages: ";
        for(auto& child: sorted_children)
            std::cerr<<child->getAge()<<"  ";
        std::cerr<<"\n";
        std::cerr<<"    parent age: "<<getAge()<<"\n";

        // The sorted_children from index `used_children` to the end have not been seen yet.
        int used_children = 1;
                
        std::cerr<<"    event=LEAF  time = "<<sorted_children[0]->getAge()<<"  active_children.size() = "<<active_children.size()<<"\n";

        // Apply the coalescence events
        for (int i=0; i<coalescence_times.size();i++)
        {
            std::cerr<<"i = "<<i<<"/"<<coalescence_times.size()<<"  active_children.size() = "<<active_children.size()<<"   children.size() = "<<children.size()<<"\n";
            // get the age of the current child
            double next_coal_time = coalescence_times[i];

            // If the next event is adding a leaf instead of a coalescent, then add the next leaf.
            for(;used_children < sorted_children.size() and sorted_children[used_children]->getAge() < next_coal_time; used_children++)
            {
                active_children.push_back( sorted_children[used_children] );
                std::cerr<<"    event=LEAF  time = "<<sorted_children[used_children]->getAge()<<"  active_children.size() = "<<active_children.size()<<"\n";
            }

            assert(active_children.size() >= 2);

            // randomly draw one child (arbitrarily called left) node from the list of active children
            size_t left = static_cast<size_t>( floor( rng->uniform01() * active_children.size() ) );
            TopologyNode* leftChild = active_children.at(left);

            // remove the randomly drawn node from the list
            active_children.erase( active_children.begin() + std::int64_t(left) );

            // randomly draw one child (arbitrarily called right) node from the list of active children
            size_t right = static_cast<size_t>( floor( rng->uniform01() * active_children.size() ) );
            TopologyNode* rightChild = active_children.at(right);

            // remove the randomly drawn node from the list
            active_children.erase( active_children.begin() + std::int64_t(right) );

            // remove the two also from the list of the children of the current node
            int old_size = children.size();
            children.erase( std::remove(children.begin(), children.end(), leftChild), children.end() );
            children.erase( std::remove(children.begin(), children.end(), rightChild), children.end() );

            // create a parent for the two
            TopologyNode* prnt = new TopologyNode(); // leave the new node without index
            prnt->setAge( next_coal_time );
            prnt->addChild( leftChild );
            prnt->addChild( rightChild );
            leftChild->setParent( prnt );
            rightChild->setParent( prnt );
            active_children.push_back( prnt );

            // add the newly created parent to the list of the children of the current node
            addChild( prnt );
            prnt->setParent( this );

            std::cerr<<"    event=COAL  time = "<<next_coal_time<<"  active_children.size() = "<<active_children.size()<<"   children.size() = "<<children.size()<<"\n";
            std::cerr<<"          left = "<<leftChild->getAge()<<"  right = "<<rightChild->getAge()<<"   parent = "<<next_coal_time<<"\n";

            assert(children.size()+1 == old_size);
            assert(prnt->getAge() > leftChild->getAge());
            assert(prnt->getAge() > rightChild->getAge());
        }
        assert(children.size() == 2);
    }
    else // use_ages == false
    {
        double brlen = getBranchLength();
                
        // The following is adapted from UniformTopologyDistribution::simulateClade()
                
        while ( children.size() >= 2 )
        {
            active_children = children;
            std::cerr<<"    branch-length tree:  active_children.size() = "<<active_children.size()<<"\n";
                    
            // randomly draw one child (arbitrarily called left) node from the list of active children
            size_t left = static_cast<size_t>( floor( rng->uniform01() * active_children.size() ) );
            TopologyNode* leftChild = active_children.at(left);
                    
            // remove the randomly drawn node from the list
            active_children.erase( active_children.begin() + std::int64_t(left) );
                    
            // randomly draw one child (arbitrarily called left) node from the list of active children
            size_t right = static_cast<size_t>( floor( rng->uniform01() * active_children.size() ) );
            TopologyNode* rightChild = active_children.at(right);
                    
            // remove the randomly drawn node from the list
            active_children.erase( active_children.begin() + std::int64_t(right) );
                    
            // remove the two also from the list of the children of the current node
            children.erase( std::remove(children.begin(), children.end(), leftChild), children.end() );
            children.erase( std::remove(children.begin(), children.end(), rightChild), children.end() );
                    
            // create a parent for the two
            TopologyNode* prnt = new TopologyNode(); // leave the new node without index
            prnt->setBranchLength(0.0);              // set the length of the branch subtending it to zero
            prnt->addChild( leftChild );
            prnt->addChild( rightChild );
            leftChild->setParent( prnt );
            rightChild->setParent( prnt );
            // we don't need active_children.push_back( prnt ) here because of the first line inside of this loop
                    
            // add the newly created parent to the list of the children of the current node
            addChild( prnt );
            prnt->setParent( this );
        }
                
        // Give my only child my old branch length, and set my new branch length to 0
        // This is to make sure everything is handled properly when we call suppressOutdegreeOneNodes() on myself
        children[0]->setBranchLength(brlen);
        setBranchLength(0.0);
    }
    std::cerr<<"DONE:  children.size() = "<<children.size()<<"\n";
}


void TopologyNode::scaleAgesFromTaxonAgesMBL(double minbl)
{
    //    1. get all my children
    //    2. if (a child is a tip)
    //           set its age based on the taxon it contains
    //       else
    //           call myself on the child
    //    3. collect the ages of all my children
    //    4. set my age to the age of the oldest of my children + min br. len.
    //    5. recompute the branch lengths of all my children
    
    const std::vector<TopologyNode*>& children = getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        if ( children[i]->isTip() )
        {
            double tip_age = ( children[i]->getTaxon().getMinAge() + children[i]->getTaxon().getMaxAge() ) / 2;
            children[i]->setAge( tip_age );
        }
        else
        {
            children[i]->scaleAgesFromTaxonAgesMBL(minbl);
        }
    }
    
    std::vector<double> ages;
    for (size_t i = 0; i < children.size(); i++)
    {
        ages.push_back( children[i]->getAge() );
    }
    
    double max_age = *std::max_element(ages.begin(), ages.end());
    setAge(max_age + minbl);
    
    // now we need to recompute the branch lengths of my children
    for (size_t i = 0; i < children.size(); i++)
    {
        children[i]->recomputeBranchLength();
    }
}


void TopologyNode::setAge(double a, bool propagate)
{
    if(getTaxon().getName() != "" && getTaxon().getMinAge() != getTaxon().getMaxAge()) {
        if(a < getTaxon().getMinAge() || a > getTaxon().getMaxAge()) {
            std::cerr << "Attempting to set new age of taxon " << getTaxon().getName() << " incompatible with age range" << std::endl;

            // NOTE: This code intentionally constructs trees with different ages for sampled-ancestors-parents and sampled-ancestor-tips.
	    //
	    // If we try to set the age of a sampled-ancestor parent to an age outside of the age rate, then this code will leave the
	    // sampled-ancestor-parent and sampled-ancestor-top with different ages.
	    //
            // If a proposal is setting the age, then the proposed should have a -Inf probability and be rejected.
	    //
            // However, if this occurs somewhere else (such as during undoProposal), then this can leave the tree in an inconsistent state
	    // So: AVOID situations where undoProposal tries to set the age of a sampled-ancestor-parent to an age outside of the age
	    //     range.

            return;
            //throw RbException() << "New age of taxa " << getTaxon().getName() << " incompatible with age range";
        }
    }
    if ( sampled_ancestor_tip == true && propagate == true )
    {
        parent->setAge(a);
        return;
    }

    age = a;

    // we need to recompute my branch-length
    recomputeBranchLength();
    
    // fire tree change event
    // we need to also flag this node as dirty (instead of only its children) as
    // 1) this node can be a tip, and
    // 2) not necessarily every flagging mechanism works recursively (e.g., flagging nodes for P matrices recomputation)
    if ( tree != NULL )
    {
        tree->getTreeChangeEventHandler().fire( *this, RevBayesCore::TreeChangeEventMessage::BRANCH_LENGTH );
    }

    // we also need to recompute the branch lengths of my children
    for (std::vector<TopologyNode *>::iterator it = children.begin(); it != children.end(); ++it)
    {
        TopologyNode *child = *it;
        if ( child->isSampledAncestorTip() )
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


void TopologyNode::setParent(TopologyNode* p, bool recompute_branch_length)
{

    // we only do something if this isn't already our parent
    if (p != parent)
    {
        // we do not own the parent so we do not have to delete it
        parent = p;
        
        if (recompute_branch_length == true)
        {
            // we need to recompute our branch length
            recomputeBranchLength();
            
            // fire tree change event
            if ( tree != NULL )
            {
                tree->getTreeChangeEventHandler().fire( *this, RevBayesCore::TreeChangeEventMessage::DEFAULT );
            }
        }
        else
        {
            // fire tree change event
            if ( tree != NULL && parent != NULL )
            {
                tree->getTreeChangeEventHandler().fire( *parent, RevBayesCore::TreeChangeEventMessage::TOPOLOGY );
            }
        }

    }
}


void TopologyNode::setSampledAncestor(bool tf)
{
    sampled_ancestor_tip = tf;

    // Only tips can have the sampled_ancestor_tip flag set.
    assert(not sampled_ancestor_tip or isTip());
}


void TopologyNode::setSpeciesName(std::string const &n)
{

    taxon.setSpeciesName( n );

}


void TopologyNode::setTaxon(Taxon const &t)
{

    taxon = t;

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


/**
 * If this node has an outdegree of 1 (i.e., only one descendant), either replace it
 * by a bifurcation where one child is subtended by a zero-length branch, and do
 * the same for all its children by calling this function recursively, or remove the
 * node entirely.
 */
void TopologyNode::suppressOutdegreeOneNodes( bool replace )
{
    
    if ( getNumberOfChildren() == 1 )
    {

        /* Should we remove the node or replace it by a bifurcation with a zero-length branch instead?
         * I.e.,
         *
         *    A                   A                 A
         *    |                   |                 |
         *    B      -->      B --+--\      or      |
         *    |                      |              |
         *    C                      C              C
         *
         *                 (replace = true)  (replace = false)
         */
        if (replace) // this solution is general enough to handle sampled ancestor root nodes
        {
            TopologyNode *new_fossil = new TopologyNode( getTaxon() );
            taxon = Taxon("");

            // connect to the old sampled ancestor
            addChild( new_fossil );
            new_fossil->setParent( this );

            // set the sampled ancestor flags
            setSampledAncestor( false );
            new_fossil->setSampledAncestor( true );

            // set the age and branch length of the newly added tip
            new_fossil->setAge( age );
            new_fossil->setBranchLength( 0.0 );
        }
        else // this only works for non-root nodes; the root is handled at the Tree level instead
        {
            // we are going to delete myself by connecting my parent and my child
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
    
    if (replace)
    {
        // call this function recursively for all children of this node
        for (size_t i = 0; i < getNumberOfChildren(); ++i)
        {
            getChild( i ).suppressOutdegreeOneNodes( true );
        }
    }

}


std::pair<double,double> getStartEndAge(const RevBayesCore::TopologyNode& node)
{
    double end_age = node.getAge();

    if (not RbMath::isFinite( end_age ))
    {
        // we assume by default that the end is at time 0
        end_age = 0;
    }

    double branch_length = node.getBranchLength();

    // From recursivelyDrawStochasticCharacterMap
    if (branch_length < 0.0)
    {
        branch_length = 1.0;
    }

    // This works for the root node.
    double start_age = end_age + branch_length;

    return {start_age, end_age};
}

