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
#include <vector>

#include <boost/unordered_set.hpp>


namespace RevBayesCore {

    namespace TreeUtilities {

        // these function are for public use
        double                  computeRobinsonFouldDistance(const Tree &a, const Tree &b);
        Tree*                   convertTree(const Tree &t, bool resetIndex=true);
//        AdmixtureTree*  convertToAdmixtureTree(const Tree &t, std::vector<std::string> names); // , const std::vector<std::string> names);
        void                    getAges(Tree *t, TopologyNode *n, std::vector<double>& ages, bool internalsOnly=true);
        AverageDistanceMatrix   getAverageDistanceMatrix(const RbVector<DistanceMatrix>& matvect);
        DistanceMatrix*         getDistanceMatrix(const Tree& tree);
        size_t                  getMrcaIndex(const TopologyNode *l, const TopologyNode *r);
        int                     getNodalDistance(const TopologyNode *l, const TopologyNode *r);
        DistanceMatrix*         getNodalDistanceMatrix(const Tree& tree);
        void                    getOldestTip(Tree* t, TopologyNode *n, double& oldest);
        void                    getTaxaInSubtree(TopologyNode *n, std::vector<TopologyNode*> &taxa );
        void                    offsetTree(Tree *t, TopologyNode *n, double factor);
        void                    makeUltrametric(Tree *t); // make the tree ultrametric by extending terminal branches
        void                    rescaleSubtree(Tree *t, TopologyNode *n, double factor, bool v=false);
        void                    rescaleTree(Tree *t, TopologyNode *n, double factor);
        void                    setAges(Tree *t, TopologyNode *n, std::vector<double>& ages);
        void                    setAgesRecursively(Tree *t, TopologyNode *n, double age);
        void                    setBranchLength(Tree *t, size_t index, double value);

        // internal helper functions
        void                    constructTimeTreeRecursively(TopologyNode *tn, const TopologyNode &n, std::vector<TopologyNode*> &nodes, std::vector<double> &ages, double depth);
//        void            constructAdmixtureTreeRecursively(AdmixtureNode *tn, const TopologyNode &n, std::vector<AdmixtureNode*> &nodes, std::vector<double> &ages);


		void                    processDistsInSubtree(const TopologyNode& node, MatrixReal& matrix, std::vector< std::pair<std::string, double> >& distsToNodeFather, const std::map< std::string, int >& namesToId);
        void                    climbUpTheTree(const TopologyNode& node, boost::unordered_set <const TopologyNode* >& pathFromNodeToRoot) ;
        double                  getAgeOfMRCARecursive(const TopologyNode& node, boost::unordered_set <const TopologyNode* >& pathFromOtherNodeToRoot) ;
        double                  getAgeOfMRCA(const Tree &t, std::string first, std::string second) ;
        int                     getCollessMetric(const TopologyNode&, int& size);
        double                  getGammaStatistic(const Tree &t);
        int                     getFitchScore(const Tree &t, const AbstractHomologousDiscreteCharacterData &c);
        std::set<size_t>        recursivelyComputeFitch(const TopologyNode &node, const AbstractHomologousDiscreteCharacterData &c, size_t site, int &score);
        double                  getMeanInverseES(const Tree &t, const AbstractHomologousDiscreteCharacterData &c, size_t state_index);
        std::vector<double>     getInverseES(const Tree &t);
        std::vector<double>     getPSSP(const Tree &t, const AbstractHomologousDiscreteCharacterData &c, size_t state_index);
        std::set<size_t>        recursivelyGetPSSP(const TopologyNode &node, const AbstractHomologousDiscreteCharacterData &c, std::vector<double> &branch_lengths, size_t state_index);
        double                  calculateMPD(const Tree &t, const AbstractHomologousDiscreteCharacterData &c, size_t site_index, size_t state_index, bool zscore, bool branch_lengths, size_t num_randomizations);
        double                  calculateMNTD(const Tree &t, const AbstractHomologousDiscreteCharacterData &c, size_t site_index, size_t state_index, bool zscore, bool branch_lengths, size_t num_randomizations);
        std::vector<double>     calculateEDR(Tree &t);


    }

}


#endif
