/**
 * @file
 * This file contains the declaration of a TopologyNode. Tree nodes are member objects and therefore can hold
 * all the variables one might to associate to the tree. The tree nodes are used to create the structure
 * of a tree. They provide access to their parent and children.
 *
 * The usage of tree nodes is to create and give easy access to parts of the dag, namely variables hold
 * by the nodes which need access to ancestral variables.
 *
 * We do not distinguish between branch parameter and node parameters. A branch parameter is simply set as
 * a node parameter of the descending node.
 *
 * A tree node can have a distribution associated with it, for instance in the situation when we condition
 * on a group of taxa being monophyletic.
 *
 *
 * NOTE: This class might be a temporary solution being the unified solution for all tree nodes. In the future
 * we might make this class abstract and implement at least the two types of tree nodes: bifurcating tree nodes
 * which are restricted two have exactly two descendants and multifurcating tree nodes which can have any number
 * of tree nodes. Perhaps there might be also tip nodes as a derived class and a root node, which is the
 * correct OO design approach but comes with a lot of overhead. A fast conversion method would be needed.
 *
 * @brief Declaration of TopologyNode
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-17 10:31:20 +0200 (Tue, 17 Jul 2012) $
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-12-04, version 1.0
 * @extends MemberObject
 *
 * $Id: TopologyNode.h 1682 2012-07-17 08:31:20Z hoehna $
 */

#ifndef TopologyNode_H
#define TopologyNode_H

#include "Clade.h"
#include "RbBitSet.h"
#include "RbConstants.h"

#include "TreeChangeEventMessage.h"
#include "Taxon.h"

#include <vector>
#include <map>
#include <string>
#include <boost/optional.hpp>

namespace RevBayesCore {
    
    class Tree;
    
    class TopologyNode  {
        
    public:
        TopologyNode();                                                                                                                 //!< Default constructor with no index
        TopologyNode(size_t indx);                                                                                                      //!< Default constructor with index
        TopologyNode(const std::string& n, const boost::optional<size_t>& indx = {});                                                   //!< Constructor with name and optional index
        TopologyNode(const Taxon& t, const boost::optional<size_t>& indx = {});                                                         //!< Constructor with taxon and optional index
        TopologyNode(const TopologyNode &n);                                                                                            //!< Copy constructor
        virtual                                    ~TopologyNode(void);                                                                 //!< Destructor
        TopologyNode&                               operator=(const TopologyNode& n);


        // Basic utility functions
        TopologyNode*                               clone(void) const;                                                                  //!< Clone object
        bool                                        equals(const TopologyNode& node) const;                                             //!< Test whether this is the same node

        // public methods
        void                                        addChild(TopologyNode* c, size_t pos = 0);                                          //!< Adds a child node

        void                                        addBranchParameter(const std::string &n, double p);
        void                                        addBranchParameter(const std::string &n, const std::string &p);
        void                                        addBranchParameters(const std::string &n, const std::vector<double> &p, bool io);
        void                                        addBranchParameters(const std::string &n, const std::vector<std::string> &p, bool io);

        void                                        addNodeParameter(const std::string &n, double p);
        void                                        addNodeParameter(const std::string &n, const std::string &p);
        void                                        addNodeParameter_(const std::string &n, const std::string &p);
        void                                        addNodeParameters(const std::string &n, const std::vector<double> &p, bool io);
        void                                        addNodeParameters(const std::string &n, const std::vector<std::string> &p, bool io);

        bool                                        hasNodeComment(const std::string& comment) const;                                   //!< Checks for a comment -- maybe not of the form key=value.
        bool                                        setNodeParameter(const std::string& name, const std::string& value);                //!< Adds OR REPLACES a node parmaeter, returning true if it was already present.
        std::optional<std::string>                  getNodeParameter(const std::string& name) const;                                    //!< Gets the value of a node parameter if present.
        std::optional<std::string>                  eraseNodeParameter(const std::string& name);                                        //!< Erases a node parameter if present, returning the value.


        void                                        clearParameters(void);                                                              //!< Clear the node and branch parameters
        void                                        clearBranchParameters(void);
        void                                        clearNodeParameters(void);
        virtual std::string                         computeNewick(bool round = true);                                                   //!< Compute the newick string for this clade
        std::string                                 computePlainNewick(void) const;                                                     //!< Compute the newick string for this clade as a plain string without branch length
        std::string                                 computeSimmapNewick(bool round = true);                                             //!< Compute the newick string compatible with SIMMAP and phytools
        bool                                        containsClade(const TopologyNode* c, bool strict) const;
        bool                                        containsClade(const Clade &c, bool strict) const;
        bool                                        containsClade(const RbBitSet &c, bool strict) const;
        bool                                        doesUseAges(void) const;                                                            //!< Does this node use ages or branch lengths?
        std::string                                 fillCladeIndices(std::map<std::string,size_t> &clade_index_map) const;              //!< Fill this map recursively with all clade indices.
        void                                        fireTreeChangeEvent(const unsigned& m = RevBayesCore::TreeChangeEventMessage::DEFAULT);
        double                                      getAge(void) const;                                                                 //!< Get the age (time ago from present) for this node
        RbBitSet                                    getAllClades(std::vector<RbBitSet> &taxa, size_t n, bool io) const;                 //!< Fill all the taxon bitset
        const std::vector<std::string>&             getBranchParameters(void) const;                                                    //!< Get the branch length leading towards this node
        double                                      getBranchLength(void) const;                                                        //!< Get the branch length leading towards this node
        size_t                                      getCladeIndex(const TopologyNode* c) const;
        const TopologyNode&                         getChild(size_t i) const;                                                           //!< Returns the i-th child
        TopologyNode&                               getChild(size_t i);                                                                 //!< Returns the i-th child (non-const to return non-const node)
        const std::vector<TopologyNode*>&           getChildren(void) const;
        std::vector<int>                            getChildrenIndices(void) const;                                                     //!< Return children indices
        Clade                                       getClade(void) const;                                                               //!< Get the clade this node represents
        bool                                        hasIndex(void) const;                                                               //!< Does the node have an index
        size_t                                      getIndex(void) const;                                                               //!< Get index of node
        void                                        getIndicesOfNodesInSubtree(bool countTips, std::vector<size_t>* indices) const;                                                               //!< Get index of node
        std::string                                 getIndividualName() const;                                                          //!< Get the species name for the node
        double                                      getMaxDepth(void) const;                                                            //!< Get the maximum depth from this node (time between this node and most recent tip)
        const std::string&                          getName() const;                                                                    //!< Get name of node
        TopologyNode*                               getMrca(const Clade &c);
        const TopologyNode*                         getMrca(const Clade &c) const;
        const TopologyNode*                         getMrca(const Clade &c, bool strict) const;
        const TopologyNode*                         getMrca(const TopologyNode &n) const;
//        const TopologyNode*                         getMrca(const std::vector<Taxon> &t) const;
        TopologyNode*                               getNode(const Clade &c, bool strict);
        TopologyNode*                               getNode(const RbBitSet &c, bool strict);
        TopologyNode*                               getNode(const TopologyNode &n, bool strict);
        const TopologyNode*                         getNode(const Clade &c, bool strict) const;
        const TopologyNode*                         getNode(const RbBitSet &c, bool strict) const;
        const TopologyNode*                         getNode(const TopologyNode &n, bool strict) const;
        const std::vector<std::string>&             getNodeParameters(void) const;                                                      //!< Get the branch length leading towards this node
        size_t                                      getNumberOfChildren(void) const;                                                    //!< Returns the number of children
        size_t                                      getNumberOfShiftEvents(void) const;
        size_t                                      getNumberOfNodesInSubtree(bool tips) const;
        size_t                                      getDegree() const;                                                                  //!< Returns the degree of the node
        TopologyNode&                               getParent(void);                                                                    //!< Returns the node's parent
        const TopologyNode&                         getParent(void) const;                                                              //!< Returns the node's parent
        std::string                                 getSpeciesName() const;                                                             //!< Get the species name for the node
        void                                        getTaxa(std::vector<Taxon> &taxa) const;                                            //!< Fill the vector of taxa
        void                                        getTaxa(RbBitSet &taxa) const;                                                      //!< Fill the taxon bitset
        void                                        getTaxa(std::vector<Taxon> &taxa, RbBitSet &bitset) const;                          //!< Fill the vector of taxa and the taxon bitset
        Taxon&                                      getTaxon();                                                                         //!< Get the taxon for this node
        const Taxon&                                getTaxon() const;                                                                   //!< Get the taxon for this node
        std::vector<double>                         getTimeInStates();
        double                                      getTmrca(const Clade &c) const;
        double                                      getTmrca(const TopologyNode &n) const;
        double                                      getTmrca(const std::vector<Taxon> &t) const;
        bool                                        isFossil(void) const;                                                               //!< Is node a fossil?
        bool                                        isInternal(void) const;                                                             //!< Is node internal?
        bool                                        isRoot(void) const;                                                                 //!< Is node root?
        bool                                        isSampledAncestorTip() const;                                                       //!< Is node a tip sampled ancestor?
        bool                                        isSampledAncestorParent() const;                                                    //!< Is child node a tip a sampled ancestor?
        bool                                        isSampledAncestorTipOrParent() const;                                               //!< Is node or child node a tip a sampled ancestor?
        bool                                        isSampledAncestorKnuckle() const;                                                   //!< Does this one have only one child?
        bool                                        isTip(void) const;                                                                  //!< Is node tip?
        void                                        makeBifurcating(bool as_fossils);                                                   //!< Make this and all its descendants bifurcating.
        void                                        recomputeAge(bool recursive);                                                       //!< Recompute the age of this node based on the childs age and the branch length leading to it.
        void                                        recomputeBranchLength(void);                                                        //!< Recompute the length of this branch based on the ages.
        void                                        renameNodeParameter(const std::string &old_name, const std::string &new_name);
        void                                        removeAllChildren(void);                                                            //!< Removes all of the children of the node
        size_t                                      removeChild(TopologyNode* c);                                                       //!< Removes a specific child
        void                                        removeTree(Tree *t);                                                                //!< Removes the tree pointer
        void                                        resolveMultifurcation(bool resolve_root = true);                                    //!< If node has more than 2 children, randomly resolve them into a bifurcating tree
        void                                        setAge(double a, bool propagate = true );                                           //!< Set the age of this node (should only be done for tips).
        void                                        setBranchLength(double b, bool flag_dirty=true);                                    //!< Set the length of the branch leading to this node.
        void                                        setIndex(size_t idx);                                                               //!< Set the index of the node

        void                                        setName(const std::string& n);                                                      //!< Set the name of this node
        void                                        setNumberOfShiftEvents(size_t n);                                                   //!< Set the number of shift events for stochastic character maps
        void                                        setParent(TopologyNode* p);                                                         //!< Sets the node's parent
        void                                        setSampledAncestor(bool tf);                                                        //!< Set if the node is a sampled ancestor
        void                                        setSpeciesName(std::string const &n);                                               //!< Set the species name of this node
        void                                        setTaxon(Taxon const &t);                                                           //!< Set the taxon of this node
        void                                        setTimeInStates(std::vector<double> t);
        void                                        setTree(Tree *t);                                                                   //!< Sets the tree pointer
        void                                        setUseAges(bool tf, bool recursive);

        // internal helper functions
        bool getBurstSpeciation(void) const { return burst_speciation; }
        bool getSamplingEvent(void) const { return sampling_event; }
        bool getSerialSampling(void) const { return serial_sampling; }
        bool getSerialSpeciation(void) const { return serial_speciation; }
        void setBurstSpeciation(bool tf) { burst_speciation = tf; }
        void setSamplingEvent(bool tf) { sampling_event = tf; }
        void setSerialSampling(bool tf) { serial_sampling = tf; }
        void setSerialSpeciation(bool tf) { serial_speciation = tf; }

    protected:

        bool burst_speciation = false;
        bool sampling_event = false;
        bool serial_sampling = false;
        bool serial_speciation = false;

        // helper methods
        std::ostream&                               buildNewick(std::ostream&, bool simmap);                                            //!< compute the newick string for a tree rooting at this node
        std::string                                 buildNewickString(bool simmap, bool round);                                         //!< compute the newick string for a tree rooting at this node

        // protected members
        bool                                        use_ages = true;
        double                                      age = RbConstants::Double::nan;
        double                                      branch_length = RbConstants::Double::nan;
        std::vector<TopologyNode*>                  children;                                                                           //!< Vector holding the node's children. Note that the parent owns the children but not the other way around.
        TopologyNode*                               parent = nullptr;                                                                   //!< Pointer to the parent of the node. It is a regular pointer instead of a super smart pointer to avoid loops in the reference counting.
        Tree*                                       tree = nullptr;                                                                     //!< A pointer to the tree for convinience access
        Taxon                                       taxon;                                                                              //!< Taxon of the node, i.e. identifier/taxon name, plus species it comes from

        boost::optional<size_t>                     index;                                                                              //!< Node index
        bool                                        sampled_ancestor_tip = false;

        // information for newick representation
        std::vector<std::string>                    node_comments;
        std::vector<std::string>                    branch_comments;

        // for stochastic maps
        std::vector<double>                         time_in_states;
        size_t                                      num_shift_events = 0;

//        RevLanguage::RevPtr<TaxonMap>               taxon_map;

     // std::map<std::string,std::string>           nodeFields;
     // std::map<std::string,std::string>           branchFields;
    };
}

std::pair<double,double> getStartEndAge(const RevBayesCore::TopologyNode& node);

#endif
