/**
 * @file
 * This file contains helper functions for manipulating trees in RevBayes.
 *
 * @brief Namespace containing helper functions for trees
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-07-05 16:47:08 +0200 (Thu, 05 Jul 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id: RlStringUtilities.h 1651 2012-07-05 14:47:08Z hoehna $
 */

#ifndef TreeUtilities_H
#define TreeUtilities_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "AverageDistanceMatrix.h"
#include "DistanceMatrix.h"
#include "MatrixReal.h"
#include "RbVector.h"
#include "Tree.h"
#include "TopologyNode.h"
#include <string>
#include <cstdint>
#include <vector>

#include <boost/unordered_set.hpp>


namespace RevBayesCore {

    namespace TreeUtilities {

        // these function are for public use
        std::vector<double>     calculateEDR(const Tree& t);                                                                            //!< get distribution of weighted times between speciation/coalescent events in a tree
        double                  calculateMPD(const Tree& t, const AbstractHomologousDiscreteCharacterData& c, size_t site_index, size_t state_index, bool zscore, bool branch_lengths, size_t num_randomizations); //!< calculate the Mean Pairwise Distance
        double                  calculateMNTD(const Tree& t, const AbstractHomologousDiscreteCharacterData& c, size_t site_index, size_t state_index, bool zscore, bool branch_lengths, size_t num_randomizations); //!< calculate the Mean Nearest Taxon Distance
        void                    climbUpTheTree(const TopologyNode& node, boost::unordered_set< const TopologyNode* >& node_root_path) ; //!< find path from given node to root
        double                  computeRobinsonFouldDistance(const std::vector<RbBitSet>& a, const std::vector<RbBitSet>& b, bool symmetric);//!< Robinson-Foulds distance
        double                  computeRobinsonFouldDistance(const Tree& a, const Tree& b, bool symmetric);                             //!< Robinson-Foulds distance
        Tree*                   convertTree(const Tree& t, bool resetIndex=true);                                                       //!< convert tree to time tree
        double                  getAgeOfMRCA(const Tree& t, const std::string& first, const std::string& second);                       //!< calculate age of MRCA based on tip names
        void                    getAges(const TopologyNode& n, std::vector<double>& ages, bool internals_only=true);                    //!< fill vector with node ages
        AverageDistanceMatrix   getAverageDistanceMatrix(const RbVector<DistanceMatrix>& matvect, const RbVector<double>* weights);     //!< calculate the (possibly weighted) average of multiple distance matrices
        int                     getCollessMetric(const TopologyNode&, int& size);                                                       //!< calculate the Colless metric
        DistanceMatrix*         getDistanceMatrix(const Tree& tree);                                                                    //!< get matrix of all distances
        int                     getFitchScore(const Tree& t, const AbstractHomologousDiscreteCharacterData &c);                         //!< calculate the parsimony score
        double                  getGammaStatistic(const Tree& t);                                                                       //!< calculate the Gamma statistic
        std::vector<double>     getInverseES(const Tree& t);                                                                            //!< calculate the inverse Equal Splits measure
        double                  getMeanInverseES(const Tree &t, const AbstractHomologousDiscreteCharacterData &c, size_t state_index);  //!< calculate the mean inverse Equal Splits measure
        size_t                  getMrcaIndex(const TopologyNode& l, const TopologyNode& r);                                             //!< get index of MRCA
        int                     getNodalDistance(const TopologyNode& l, const TopologyNode& r);                                         //!< get nodal distance between two nodes
        DistanceMatrix*         getNodalDistanceMatrix(const Tree& tree);                                                               //!< get matrix of nodal distances between all tips
        double                  getOldestTipAge(const TopologyNode& n);                                                                 //!< get the age of the oldest tip below specified node
        std::vector<double>     getPSSP(const Tree& t, const AbstractHomologousDiscreteCharacterData& c, size_t state_index);           //!< calculate the Parsimoniously Same State Paths
        void                    getTaxaInSubtree(TopologyNode& n, std::vector<TopologyNode*>& taxa );                                   //!< get taxa below specified node
        bool                    isConnectedNNI(const Tree& a, const Tree& b);                                                           //!< Check if the two trees are connected by a single NNI move
        void                    makeUltrametric(Tree& t);                                                                               //!< make the tree ultrametric by extending terminal branches
        void                    offsetTree(TopologyNode& n, double factor);                                                             //!< offset node and its children by a factor
        void                    rescaleSubtree(TopologyNode& n, double factor, bool v=false);                                           //!< rescale tree ages below a node by a factor, except tips
        void                    rescaleTree(TopologyNode& n, double factor);                                                            //!< rescale tree ages below a node by a factor
        void                    setAges(TopologyNode& n, const std::vector<double>& ages);                                              //!< set ages of a node and children from a vector
        void                    setAgesRecursively(TopologyNode& n, double age);                                                        //!< set age of a node and rescale its children
        Tree*                   startingTreeInitializer(Tree& treeToChange, std::vector<Taxon>& taxaToCopy, std::int64_t agePrecision);         //!< make sure starting tree satisfies age constraints, and assign min/max ages to its tips
        
        // internal helper functions
        void                    constructTimeTreeRecursively(TopologyNode& tn, const TopologyNode &n, std::vector<TopologyNode*> &nodes, std::vector<double> &ages, double depth); //!< helper function for time tree conversion
        void                    processDistsInSubtree(const TopologyNode& node, MatrixReal& matrix, std::vector< std::pair<std::string, double> >& distsToNodeFather, const std::map< std::string, int >& namesToId); //!< helper function for distance matrix calculation
        double                  getAgeOfMRCARecursive(const TopologyNode& node, boost::unordered_set <const TopologyNode* >& pathFromOtherNodeToRoot) ; //!< helper for MRCA age
        std::set<size_t>        recursivelyGetPSSP(const TopologyNode &node, const AbstractHomologousDiscreteCharacterData &c, std::vector<double> &branch_lengths, size_t state_index); //!< helper for PSSP calculation
        std::set<size_t>        recursivelyComputeFitch(const TopologyNode &node, const AbstractHomologousDiscreteCharacterData &c, size_t site, int &score); //!< helper for parsimony score

    }

}


#endif
