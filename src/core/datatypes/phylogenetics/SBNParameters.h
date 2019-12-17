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
        const std::string &                                                 getBranchLengthApproximationMethod(void) const; // Tell someone how this SBN is approximating branch lengths
        std::map<RbBitSet,std::vector<double> >&                            getEdgeLengthDistributionParameters(void); // For each split in an unrooted tree, the parameters of the edge length distribution
        std::map<RbBitSet,std::vector<double> >&                            getNodeTimeDistributionParameters(void); // For each clade in a rooted tree, the parameters of the edge length distribution
        const size_t                                                        getNumTaxa(void) const; // The number of taxa in the tree the SBN describes
        std::vector<std::pair<Subsplit,double> >&                           getRootSplits(void); // The root splits in the tree and their probabilities
        std::map<Subsplit,std::vector<std::pair<Subsplit,double> > >&       getSubsplitCPDs(void); // For each subsplit, its children and their probabilities
        std::vector<Taxon>&                                                 getTaxa(void); // The taxa in the tree the SBN describes
        const std::vector<Taxon>&                                           getTaxa(void) const; // The taxa in the tree the SBN describes

        // Helper functions for SBN distributions
        double                              computeLnProbabilityRootedTopology( const Tree &tree ) const;
        std::vector<std::pair<Subsplit,double> >  computeLnProbabilityTopologyAndRooting( const Tree &tree ) const;
        double                              computeLnProbabilityUnrootedTopology( const Tree &tree ) const;
        double                              computeRootSplitProbability( const Subsplit &root_split ) const;
        double                              computeSubsplitTransitionProbability( const Subsplit &parent, const Subsplit &child ) const;
        Subsplit                            drawRootSplit( void ) const;
        Subsplit                            drawSubsplitForY( const Subsplit &s ) const;
        Subsplit                            drawSubsplitForZ( const Subsplit &s ) const;
        bool                                isValid(void) const;
        double                              KL(SBNParameters &Q_x) const; // Takes the current object to be P(x) and the passed argument to be Q(x) and computes KL(P || Q)

        // Functions for learning SBNs
        void                                countAllSubsplits(Tree& tree, std::map<std::pair<Subsplit,Subsplit>,double>& parent_child_counts, std::map<Subsplit,double>& root_split_counts, std::map<Subsplit,double>& q, bool doSA);
        void                                fitBranchLengthDistributions(std::vector<Tree> &trees);
        void                                fitNodeTimeDistributions(std::vector<Tree> &trees);
        bool                                isValidCPD(std::vector<std::pair<Subsplit,double> >& cpd, Subsplit& parent) const;
        bool                                isValidRootDistribution(void) const;
        void                                learnTimeCalibratedSBN( std::vector<Tree>& trees );
        void                                learnUnconstrainedSBNSA( std::vector<Tree> &trees );
        void                                learnUnconstrainedSBNEM( std::vector<Tree> &trees, double &alpha );
        void                                makeCPDs(std::map<std::pair<Subsplit,Subsplit>,double>& parent_child_counts);
        void                                makeRootSplits(std::map<Subsplit,double>& root_split_counts);
        void                                normalizeCPDForSubsplit(std::vector<std::pair<Subsplit,double> >& cpd, Subsplit& parent);
        void                                regularizeCounts(std::map<std::pair<Subsplit,Subsplit>,double>& parent_child_counts, std::map<Subsplit,double>& root_split_counts, std::map<std::pair<Subsplit,Subsplit>,double>& pseudo_parent_child_counts, std::map<Subsplit,double>& pseudo_root_split_counts, double alpha);

        // Helper functions for learning SBNs
        void                                addTreeToAllParentChildCounts(std::map<std::pair<Subsplit,Subsplit>,double>& parent_child_counts, Tree& tree, double &weight);
        void                                addTreeToAllRootSplitCounts(std::map<Subsplit,double>& root_split_counts, Tree& tree, double &weight);
        void                                incrementParentChildCount(std::map<std::pair<Subsplit,Subsplit>,double> &parent_child_counts, std::pair<Subsplit,Subsplit> &this_parent_child, double &weight);
        void                                incrementRootSplitCount(std::map<Subsplit,double>& root_split_counts, Subsplit &this_root_split, double &weight);

        // // Misc.
        std::map<Subsplit,double>            computeUnconditionalSubsplitProbabilities(void) const;
        void                                 recursivelyComputeUnconditionalSubsplitProbabilities(std::map<Subsplit,double> &subsplit_probs, Subsplit &parent, double p) const;
        std::map<Subsplit,double>            computeUnconditionalSubsplitProbabilities1(void);
        void                                 recursivelyComputeUnconditionalSubsplitProbabilities1(std::map<Subsplit,double> &subsplit_probs, Subsplit &parent, double p);
        // std::vector<std::pair<Split,double> > computeCladeProbabilities(void) const;
        // std::vector<std::pair<Split,double> > computeSplitProbabilities(void) const;

      private:
        // members
        std::string                                    branch_length_approximation_method;
        std::map<RbBitSet,std::vector<double> >        edge_length_distribution_parameters; // In an unconstrained SBN, we learn branch lengths as functions of the splits they represent
        size_t                                         num_taxa; // The number of taxa in the tree the SBN describes
        std::vector<std::pair<Subsplit,double> >       root_splits; // The root splits in the tree and their probabilities
        std::map<Subsplit,std::vector<std::pair<Subsplit,double> > >        subsplit_cpds; // For each subsplit, its children and their probabilities
        std::vector<Taxon>                             taxa; // The taxa in the tree the SBN describes
        bool                                           time_calibrated; // Is this a time-calibrated SBN?
    };

    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const SBNParameters& x);                                           //!< Overloaded output operator


}

#endif
