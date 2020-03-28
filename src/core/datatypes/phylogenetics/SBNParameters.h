/**
 * @file
 * This file contains the declaration of SBNParameters, which is
 * class that holds the parameters for an SBN, and methods for use
 * with SBNs.
 *
 * @brief Declaration of SBNParameters
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#ifndef SBNParameters_H
#define SBNParameters_H

#include "RbConstants.h"
#include "RbBitSet.h"
// #include "Split.h"
#include "Subsplit.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDistribution.h"

#include <vector>
#include <unordered_map>

namespace RevBayesCore {

    class SBNParameters {

    public:
                                           SBNParameters();  //!< Constructor
                                           SBNParameters(std::vector<Taxon> taxa, const std::string &branch_length_approximation, bool allow_unseen = false, double epsilon = 0.0);  //!< Constructor
                                           SBNParameters(const SBNParameters &sbn);
        virtual                            ~SBNParameters();


        // overloaded operators
        SBNParameters&        operator=(const SBNParameters& sbn);

        // Access to members
        const std::string &                                                           getBranchLengthApproximationMethod(void) const; // Tell someone how this SBN is approximating branch lengths
        std::unordered_map<RbBitSet,std::vector<double> >&                            getEdgeLengthDistributionParameters(void); // For each split in an unrooted tree, the parameters of the edge length distribution
        std::unordered_map<RbBitSet,std::vector<double> >&                            getNodeTimeDistributionParameters(void); // For each clade in a rooted tree, the parameters of the edge length distribution
        const size_t                                                                  getNumTaxa(void) const; // The number of taxa in the tree the SBN describes
        std::vector<Taxon>&                                                           getTaxa(void); // The taxa in the tree the SBN describes
        const std::vector<Taxon>&                                                     getTaxa(void) const; // The taxa in the tree the SBN describes

        // Helper functions for SBN distributions
        double                              computeEdgeLengthProbability( const RbBitSet &split, double length) const; // Compute probability of edge length for split
        double                              computeNodeTimeProbability( const RbBitSet &clade, double time) const; // Compute probability of node time (or node time proportion) for clade
        double                              computeLnProbabilityRootedTopology( const Tree &tree ) const;
        std::vector<std::pair<size_t,double> >  computeLnProbabilityTopologyAndRooting( const Tree &tree ) const;
        std::vector<std::pair<size_t,double> >  computeLnProbabilityTopologyAndRooting( const std::vector<std::vector<size_t> >& compact_tree ) const;
        double                              computeLnProbabilityUnrootedTopology( const Tree &tree ) const;
        double                              computeLnProbabilityUnrootedTopology( const std::vector<std::vector<size_t> >& compact_tree ) const;
        double                              computeRootSplitProbability( const size_t index ) const;
        double                              computeRootSplitProbability( const Subsplit &s ) const;
        double                              computeSubsplitTransitionProbability( const size_t parent_index, const size_t child_index ) const;
        double                              computeSubsplitTransitionProbability( const Subsplit &parent, const Subsplit &child ) const;
        double                              drawEdgeLength( const RbBitSet &split) const; // Draw a new edge length for a split
        double                              drawNodeTime( const RbBitSet &split) const; // Draw a new node time (or node time proportion) for a clade
        size_t                              drawRootSplit( void ) const;
        size_t                              drawSubsplitForY( size_t s ) const;
        size_t                              drawSubsplitForZ( size_t s ) const;
        size_t                              getIndex( const Subsplit &s ) const; // Find index of subsplit s in the master subsplit list
        Subsplit                            getSubsplit( size_t i ) const; // Find subsplit for index i in the master subsplit list
        const Subsplit&                     getSubsplitReference( size_t i ) const; // Find subsplit for index i in the master subsplit list

        // Functions for learning SBNs
        void                                fitBranchLengthDistributions(std::vector<Tree> &trees);
        void                                fitNodeTimeDistributions(std::vector<Tree> &trees);
        void                                learnTimeCalibratedSBN( std::vector<Tree>& trees );
        void                                learnUnconstrainedSBNSA( std::vector<Tree> &trees );
        void                                learnUnconstrainedSBNEM( std::vector<Tree> &trees, double &alpha, size_t min_loop_iter, size_t max_loop_iter );
        void                                makeCPDs(std::unordered_map<std::pair<size_t,size_t>,double>& parent_child_counts);
        void                                makeRootSplits(std::unordered_map<size_t,double>& root_split_counts);
        void                                normalizeCPDForSubsplit(size_t parent_index);

        // Helper functions for learning SBNs
        void                                addTreeToAllParentChildCounts(std::unordered_map<std::pair<size_t,size_t>,double>& parent_child_counts, Tree& tree, double &weight);
        void                                addTreeToAllRootSplitCounts(std::unordered_map<size_t,double>& root_split_counts, Tree& tree, double &weight);
        void                                countAllSubsplits(Tree& tree, std::unordered_map<std::pair<size_t,size_t>,double>& parent_child_counts, std::unordered_map<size_t,double>& root_split_counts, std::unordered_map<size_t,double>& q, bool doSA);
        void                                countAllSubsplits(const std::vector<std::vector<size_t> >& compact_tree, std::unordered_map<std::pair<size_t,size_t>,double>& parent_child_counts, std::unordered_map<size_t,double>& root_split_counts, std::unordered_map<size_t,double>& q, bool doSA);
        void                                enumerateAllSubsplits(std::vector<Tree> &trees);
        std::vector<std::pair<size_t,size_t> > getCasesFromCompactTree(const std::vector<std::vector<size_t> >& compact_tree, size_t row);
        void                                incrementParentChildCount(std::unordered_map<std::pair<size_t,size_t>,double> &parent_child_counts, std::pair<size_t,size_t> &this_parent_child, double &weight);
        void                                incrementRootSplitCount(std::unordered_map<size_t,double>& root_split_counts, size_t this_root_split, double &weight);
        bool                                isValid(void) const;
        bool                                isValidCPD(size_t parent_index) const;
        bool                                isValidRootDistribution(void) const;
        void                                regularizeCounts(std::unordered_map<std::pair<size_t,size_t>,double>& parent_child_counts, std::unordered_map<size_t,double>& root_split_counts, std::unordered_map<std::pair<size_t,size_t>,double>& pseudo_parent_child_counts, std::unordered_map<size_t,double>& pseudo_root_split_counts, double alpha);
        std::vector<std::pair<size_t,size_t> > subsplitCasesToIndexCases(std::vector<std::pair<Subsplit,Subsplit> > &subsplit_cases) const;
        void                                treesToCompactTrees(std::vector<Tree> &trees, std::vector<std::vector<std::vector<size_t> > >& compact_trees);
        // // Misc.
        std::vector<double>                  computeUnconditionalSubsplitProbabilities(void) const;
        std::vector<std::pair<size_t,double> >               getCPDForParent(void) const;
        std::vector<std::vector<std::pair<size_t,double> > > getCPDs(void) const;
        size_t                               getNumChildrenForParent(size_t parent) const;
        size_t                               getNumRootSplits(void) const;
        std::vector<std::pair<size_t,double> >               getRootSplits(void) const;
        double                               KL(SBNParameters &Q_x) const; // Takes the current object to be P(x) and the passed argument to be Q(x) and computes KL(P || Q)
        void                                 recursivelyComputeUnconditionalSubsplitProbabilities(std::vector<double> &subsplit_probs, size_t parent_index, double p) const;
        size_t                               size(void) const;
        // std::vector<std::pair<Split,double> > computeCladeProbabilities(void) const;
        // std::vector<std::pair<Split,double> > computeSplitProbabilities(void) const;

      private:
        // members
        // TODO: pull out root splits from all_subsplits to enable storing root splits and q as simple vectors of probabilities
        std::unordered_map<Subsplit,size_t>            all_subsplits; // Internally we work with subsplits as indices, this takes a subsplit and gives us its index
        std::string                                    branch_length_approximation_method;
        std::unordered_map<RbBitSet,std::vector<double> >        edge_length_distribution_parameters; // In an unconstrained SBN, we learn branch lengths as functions of the splits they represent
        size_t                                         num_taxa; // The number of taxa in the tree the SBN describes
        std::vector<std::pair<size_t,double> >         root_splits; // The root splits (indices) in the tree and their probabilities
        std::vector<std::vector<std::pair<size_t,double> > >     subsplit_cpds; // For each subsplit (index), its children (indices) and their (conditional) probabilities
        std::vector<Taxon>                             taxa; // The taxa in the tree the SBN describes
        bool                                           time_calibrated; // Is this a time-calibrated SBN?

    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const SBNParameters& x);                                           //!< Overloaded output operator


}

namespace std {
  template<>
  struct hash<std::pair<RevBayesCore::Subsplit,RevBayesCore::Subsplit> >
  {
    size_t operator()(const std::pair<RevBayesCore::Subsplit,RevBayesCore::Subsplit>& p) const
    {
      return p.first.getHash() ^ p.second.getHash();
    }
  };

  template<>
  struct equal_to<std::pair<RevBayesCore::Subsplit,RevBayesCore::Subsplit> >
  {
    size_t operator()(const std::pair<RevBayesCore::Subsplit,RevBayesCore::Subsplit>& lhs, const std::pair<RevBayesCore::Subsplit,RevBayesCore::Subsplit>& rhs) const
    {
      return lhs == rhs;
    }
  };

  template<>
  struct hash<std::pair<size_t,size_t> >
  {
    size_t operator()(const std::pair<size_t,size_t>& p) const
    {
      return std::hash<size_t>{}(p.first) ^ std::hash<size_t>{}(p.second);
    }
  };

  template<>
  struct equal_to<std::pair<size_t,size_t> >
  {
    size_t operator()(const std::pair<size_t,size_t>& lhs, const std::pair<size_t,size_t>& rhs) const
    {
      return lhs == rhs;
    }
  };

}

#endif
