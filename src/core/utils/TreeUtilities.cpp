#include <math.h>
#include <cstdlib>
#include <boost/functional/hash/extensions.hpp>
#include <algorithm>
#include <ostream>
#include <iostream>
#include <limits>
#include <set>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "MatrixBoolean.h"
#include "MatrixReal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "AbstractDiscreteTaxonData.h"
#include "AverageDistanceMatrix.h"
#include "Cloneable.h"
#include "DiscreteCharacterState.h"
#include "DistanceMatrix.h"
#include "RbConstIterator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "boost/unordered_set.hpp"

using namespace RevBayesCore;


/**
 * Calculate distribution of times between speciation/coalescent events in a tree,
 * where the times are normalized by the number of lineages observed at that event.
 * The moments of this distribution can be informative statistics for distinguishing different birth-death processes.
 * @param t input tree
 * @return vector of weighted times
 */
std::vector<double> RevBayesCore::TreeUtilities::calculateEDR(const Tree& t)
{
    const TopologyNode& node = t.getRoot();
    std::vector<double> ages = std::vector<double>(t.getNumberOfNodes(), 0.0);
    TreeUtilities::getAges(node, ages);
    std::sort(ages.begin(), ages.end());

    std::vector<double> edr;
    double time = 0.0;
    size_t n = t.getNumberOfTips();
    for (size_t i = 0; i < ages.size(); ++i)
    {
        if (ages[i] > 0)
        {
            edr.push_back(n * (ages[i] - time));
            time = ages[i];
            --n;
        }
    }
    return edr;
}


/**
 * Calculate the Mean Pairwise Distance (MPD; Webb 2000; Webb et al 2002) for taxa with a certain observed character state.
 * The z-score of the MPD (optionally calculated using randomizations) is equivalent to the Net Relatedness Index (NRI; Webb 2000).
 * @param t input tree
 * @param c input character state alignment
 * @param site_index index of the character site
 * @param state_index index of the specific state
 * @param zscore whether to calculate the z-score
 * @param branch_lengths whether to use branch lengths in the calculation
 * @param num_randomizations number of randomizations used to calculate the z-score
 * @return observed MPD or z-score of the MPD
 */
double RevBayesCore::TreeUtilities::calculateMPD(const Tree& t, const AbstractHomologousDiscreteCharacterData& c, size_t site_index, size_t state_index,
                                                 bool zscore, bool branch_lengths, size_t num_randomizations)
{

    // use either pairwise branch length or nodal distances
    DistanceMatrix* distances;
    if (branch_lengths == true)
    {
        distances = getDistanceMatrix(t);
    }
    else
    {
        distances = getNodalDistanceMatrix(t);
    }

    // get the mean pairwise distances of all taxa in the observed state
    double total_dist = 0.0;
    double num_dist = 0.0;
    size_t num_taxa_in_state = 0;
    std::vector<std::string> tip_names = t.getTipNames();
    for (size_t i = 0; i < tip_names.size(); ++i)
    {
        std::string name1 = distances->getTaxa()[i].getName();
        size_t state1 = c.getTaxonData(name1).getCharacter(site_index).getStateIndex();
        if (state1 == state_index)
        {
            num_taxa_in_state += 1;
            for (size_t j = i + 1; j < tip_names.size(); ++j)
            {
                std::string name2 = distances->getTaxa()[j].getName();
                size_t state2 = c.getTaxonData(name2).getCharacter(site_index).getStateIndex();
                if (state2 == state_index)
                {
                    total_dist += distances->getMatrix()[i][j];
                    num_dist += 1.0;
                }
            }
        }
    }
    double obs_mean_distance = 0.0;
    if (num_dist != 0)
    {
        obs_mean_distance = total_dist/num_dist;
    }

    if (zscore == false)
    {
        return obs_mean_distance;
    }

    // randomizations to calculate z-score (NRI)
    std::vector<double> random_mean_distances;
    for (size_t k = 0; k < num_randomizations; ++k)
    {
        std::vector<std::string> new_tip_names;

        // pick random taxa
        while (new_tip_names.size() < num_taxa_in_state)
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            size_t u = rng->uniform01() * tip_names.size();
            std::string new_name = tip_names[u];
            if (std::find(new_tip_names.begin(), new_tip_names.end(), new_name) == new_tip_names.end())
            {
                new_tip_names.push_back(new_name);
            }
        }

        // calculate mean pairwise distances for the random taxa
        double total_dist = 0.0;
        double num_dist = 0.0;
        for (size_t i = 0; i < tip_names.size(); ++i)
        {
            std::string name1 = distances->getTaxa()[i].getName();
            if (std::find(new_tip_names.begin(), new_tip_names.end(), name1) != new_tip_names.end())
            {
                for (size_t j = i + 1; j < tip_names.size(); ++j)
                {
                    std::string name2 = distances->getTaxa()[j].getName();
                    if (std::find(new_tip_names.begin(), new_tip_names.end(), name2) != new_tip_names.end())
                    {
                        total_dist += distances->getMatrix()[i][j];
                        num_dist += 1.0;
                    }
                }
            }
        }

        double temp_mean_distance = 0.0;
        if (num_dist != 0)
        {
            temp_mean_distance = total_dist/num_dist;
        }
        random_mean_distances.push_back(temp_mean_distance);
    }

    // zscore = (mpd_obs - mpd_rand_mean) / mpd_rand_sd
    double random_mean = 0.0;
    for (size_t i = 0; i < random_mean_distances.size(); ++i)
    {
        random_mean += random_mean_distances[i];
    }
    if (random_mean_distances.size() != 0)
    {
        random_mean /= random_mean_distances.size();
    }
    double random_stdv = 0.0;
    for (size_t i = 0; i < random_mean_distances.size(); ++i)
    {
        random_stdv += pow(random_mean_distances[i] - random_mean, 2);
    }
    if (random_mean_distances.size() != 0)
    {
        random_stdv /= random_mean_distances.size();
    }
    random_stdv = sqrt(random_stdv);
    if (random_stdv != 0.0)
    {
        return (obs_mean_distance - random_mean) / random_stdv;
    }
    return 0.0;
}


/**
 * Calculate the Mean Nearest Taxon Distance (MNTD; Webb 2000; Webb et al 2002) for taxa with a certain observed character state.
 * The z-score of the MNTD (optionally calculated using randomizations) is equivalent to the Nearest Taxa Index (NTI; Webb 2000).
 * @param t input tree
 * @param c input character state alignment
 * @param site_index index of the character site
 * @param state_index index of the specific state
 * @param zscore whether to calculate the z-score
 * @param branch_lengths whether to use branch lengths in the calculation
 * @param num_randomizations number of randomizations used to calculate the z-score
 * @return observed MNTD or z-score of the MNTD
 */
double RevBayesCore::TreeUtilities::calculateMNTD(const Tree& t, const AbstractHomologousDiscreteCharacterData &c, size_t site_index, size_t state_index,
                                                  bool zscore, bool branch_lengths, size_t num_randomizations)
{

    // use either pairwise branch length or nodal distances
    DistanceMatrix* distances;
    if (branch_lengths == true)
    {
        distances = getDistanceMatrix(t);
    }
    else
    {
        distances = getNodalDistanceMatrix(t);
    }

    // get the distances of the closest relative of all taxa in the observed state
    double total_dist = 0.0;
    double num_dist = 0.0;
    size_t num_taxa_in_state = 0;
    std::vector<std::string> tip_names = t.getTipNames();
    for (size_t i = 0; i < tip_names.size(); ++i)
    {
        std::string name1 = distances->getTaxa()[i].getName();
        size_t state1 = c.getTaxonData(name1).getCharacter(site_index).getStateIndex();
        if (state1 == state_index)
        {
            num_taxa_in_state += 1;
            double min_dist = 0.0;
            for (size_t j = 0; j < tip_names.size(); ++j)
            {
                if (j != i)
                {
                    std::string name2 = distances->getTaxa()[j].getName();
                    size_t state2 = c.getTaxonData(name2).getCharacter(site_index).getStateIndex();
                    if (state2 == state_index)
                    {
                        double this_dist = distances->getMatrix()[i][j];
                        if (this_dist < min_dist || min_dist == 0.0)
                        {
                            min_dist = this_dist;
                        }
                    }
                }
            }
            total_dist += min_dist;
            num_dist += 1.0;
        }
    }
    double obs_mean_distance = 0.0;
    if (num_dist != 0)
    {
        obs_mean_distance = total_dist/num_dist;
    }

    if (zscore == false)
    {
        return obs_mean_distance;
    }

    // randomizations to calculate z-score (NTI)
    std::vector<double> random_mean_distances;
    for (size_t k = 0; k < num_randomizations; ++k)
    {
        double total_dist = 0.0;
        double num_dist = 0.0;
        std::vector<std::string> new_tip_names;

        // pick random taxa
        while (new_tip_names.size() < num_taxa_in_state)
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            size_t u = rng->uniform01() * tip_names.size();
            std::string new_name = tip_names[u];
            if (std::find(new_tip_names.begin(), new_tip_names.end(), new_name) == new_tip_names.end())
            {
                new_tip_names.push_back(new_name);
            }
        }

        // calculate mean nearest taxon distance for the random taxa
        for (size_t i = 0; i < tip_names.size(); ++i)
        {
            std::string name1 = distances->getTaxa()[i].getName();
            if (std::find(new_tip_names.begin(), new_tip_names.end(), name1) != new_tip_names.end())
            {
                double min_dist = 0.0;
                for (size_t j = 0; j < tip_names.size(); ++j)
                {
                    if (j != i)
                    {
                        std::string name2 = distances->getTaxa()[j].getName();
                        if (std::find(new_tip_names.begin(), new_tip_names.end(), name2) != new_tip_names.end())
                        {
                            double this_dist = distances->getMatrix()[i][j];
                            if (this_dist < min_dist || min_dist == 0.0)
                            {
                                min_dist = this_dist;
                            }
                        }
                    }
                }
                total_dist += min_dist;
                num_dist += 1.0;
            }
        }

        double temp_mean_distance = 0.0;
        if (num_dist != 0)
        {
            temp_mean_distance = total_dist/num_dist;
        }
        random_mean_distances.push_back(temp_mean_distance);
    }

    double random_mean = 0.0;
    for (size_t i = 0; i < random_mean_distances.size(); ++i)
    {
        random_mean += random_mean_distances[i];
    }
    if (random_mean_distances.size() != 0)
    {
        random_mean /= random_mean_distances.size();
    }
    double random_stdv = 0.0;
    for (size_t i = 0; i < random_mean_distances.size(); ++i)
    {
        random_stdv += pow(random_mean_distances[i] - random_mean, 2);
    }
    if (random_mean_distances.size() != 0)
    {
        random_stdv /= random_mean_distances.size();
    }
    random_stdv = sqrt(random_stdv);
    if (random_stdv != 0.0)
    {
        return (obs_mean_distance - random_mean) / random_stdv;
    }
    return 0.0;
}


/**
 * Find path from a given node to the root, recursively
 * @param[in] node current node
 * @param[out] pathFromNodeToRoot path found so far
 */
void RevBayesCore::TreeUtilities::climbUpTheTree(const TopologyNode& node, boost::unordered_set <const TopologyNode* >& node_root_path)

{

    if ( node.isRoot() == false )
    {
        node_root_path.insert(&node);
        climbUpTheTree(node.getParent(), node_root_path);
    }
    
}

/** Calculate Robinson-Foulds distance between two trees
 * @param a,b trees between which to calculate the distance
 * @return RF distance
 */
double RevBayesCore::TreeUtilities::computeRobinsonFouldDistance(const RevBayesCore::Tree& a, const RevBayesCore::Tree& b, bool symmetric)
{
    
    std::vector<RbBitSet>* bipartitions_a = a.getNodesAsBitset();
    std::vector<RbBitSet>* bipartitions_b = b.getNodesAsBitset();
    
    double RF_distance = TreeUtilities::computeRobinsonFouldDistance(*bipartitions_a, *bipartitions_b, symmetric);
    
    delete bipartitions_a;
    delete bipartitions_b;
    
    return RF_distance;
}

/** Calculate Robinson-Foulds distance between two trees
 * @param a,b trees between which to calculate the distance
 * @return RF distance
 */
double RevBayesCore::TreeUtilities::computeRobinsonFouldDistance(const std::vector<RevBayesCore::RbBitSet>& bipartitions_a, const std::vector<RevBayesCore::RbBitSet>& bipartitions_b, bool symmetric)
{

    bool found = false;
    double distance = 0.0;
    for (size_t i = 0; i< bipartitions_a.size(); ++i)
    {
        found = false;
        for (size_t j = 0; j < bipartitions_b.size(); ++j)
        {
            if (bipartitions_a[i] == bipartitions_b[j])
            {
                found = true;
                break;
            }
        }
        if (found == false)
        {
            distance += 1.0;
        }
    }
    
    if ( symmetric == true )
    {
        distance *= 2;
    }
    else
    {

        for (size_t i = 0; i< bipartitions_b.size(); ++i)
        {
            found = false;
            for (size_t j = 0; j < bipartitions_a.size(); ++j)
            {
                if (bipartitions_b[i] == bipartitions_a[j])
                {
                    found = true;
                    break;
                }
            }

            if (found == false)
            {
                distance += 1.0;
            }
        }
        
    }

    return distance;
}

/**
 * Helper function to recusively build a time tree
 * @param tn current time tree node
 * @param n current node
 * @param nodes current vector of nodes of the time tree
 * @param ages current vector of ages of the time tree
 * @param depth depth of the current node
 */
void RevBayesCore::TreeUtilities::constructTimeTreeRecursively(TopologyNode& tn, const TopologyNode& n, std::vector<TopologyNode*> &nodes, std::vector<double> &ages, double depth)
{

    // set the name
    tn.setName( n.getName() );

    // copy the index
    tn.setIndex( n.getIndex() );

    // copy the branch "comments"
    const std::vector<std::string> &branchComments = n.getBranchParameters();
    for (size_t i = 0; i < branchComments.size(); ++i)
    {
        std::string tmp = branchComments[i];
        if ( tmp[0] == '&')
        {
            tmp = tmp.substr(1,tmp.size());
        }
        std::vector<std::string> pair;
        StringUtilities::stringSplit(tmp, "=", pair);
        tn.addBranchParameter(pair[0], pair[1]);
    }
    // copy the node "comments"
    const std::vector<std::string> &nodeComments = n.getNodeParameters();
    for (size_t i = 0; i < nodeComments.size(); ++i)
    {
        std::string tmp = nodeComments[i];
        if ( tmp[0] == '&')
        {
            tmp = tmp.substr(1,tmp.size());
        }
        std::vector<std::string> pair;
        StringUtilities::stringSplit(tmp, "=", pair);
        tn.addNodeParameter(pair[0], pair[1]);
    }

    // set the node flags
    tn.setSampledAncestor( n.isSampledAncestorTip() );

    // remember the node
    nodes.push_back( &tn );

    // set the age
    double a = depth - n.getBranchLength();
    if ( a < 1E-5 )
    {
        a = 0.0;
    }
    ages.push_back( a );

    // create children
    for (size_t i = 0; i < n.getNumberOfChildren(); ++i)
    {
        const TopologyNode& child = n.getChild( i );
        TopologyNode* new_child = new TopologyNode();

        // set parent child relationship
        new_child->setParent( &tn );
        tn.addChild( new_child );

        // start recursive call
        constructTimeTreeRecursively(*new_child, child, nodes, ages, a);
    }

    // Mark tip nodes on a 0-length branch as sampled ancestors
    if ( tn.getNumberOfChildren() == 2)
    {
        // The time-tree doesn't have any branch lengths or ages yet, so we have to look
        // at the branch-length tree that we are copying from.

        auto& bl_c1 = n.getChild(0);
        auto& bl_c2 = n.getChild(1);

        if (bl_c1.isTip() and bl_c1.getBranchLength() == 0 and bl_c2.getBranchLength() > 0)
        {
            tn.getChild(0).setSampledAncestor( true );
        }

        if (bl_c2.isTip() and bl_c2.getBranchLength() == 0 and bl_c1.getBranchLength() > 0)
        {
            tn.getChild(1).setSampledAncestor( true );
        }
    }
}

/**
 * Convert tree to time tree (with clock)
 * @param t input tree
 * @param resetIndex whether to reorder node indexes
 * @return time tree
 */
RevBayesCore::Tree* RevBayesCore::TreeUtilities::convertTree(const Tree& t, bool reset_index)
{
    // create time tree object (topology + times)
    Tree *tt = new Tree();

    // clock trees should always be rooted
    tt->setRooted( true );

    // get the root of the original tree
    const TopologyNode& bln = t.getRoot();

    TopologyNode* root = new TopologyNode();

    // copy the root index
    root->setIndex( bln.getIndex() );

    std::vector<double> ages;
    std::vector<TopologyNode*> nodes;

    double max_depth = bln.getMaxDepth() + bln.getBranchLength();

    // recursive creation of the tree
    constructTimeTreeRecursively(*root, bln, nodes, ages, max_depth);

    // add the root which creates the topology
    tt->setRoot( root, reset_index );

    // set the ages
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->setAge( ages[i] );
    }

    // copy the root edge
    root->setBranchLength( bln.getBranchLength() );

    return tt;
}


/** Get ages of all nodes below a specific node
 * @param[in] t input tree
 * @param[in] n specified node
 * @param[out] ages return vector to fill with ages
 * @param[in] internalsOnly whether to get ages for only internal nodes and not tips
 */
void RevBayesCore::TreeUtilities::getAges(const TopologyNode& n, std::vector<double>& ages, bool internals_only)
{
    // we only get internal node ages if internalsOnly==true
    if ( n.isTip() == false )
    {
        // get the age of the node
        ages[n.getIndex()] = n.getAge();

        // get all children ages
        const std::vector<TopologyNode*>& children = n.getChildren();
        for (size_t i = 0; i < children.size(); i++)
        {
            getAges(  *children[i], ages);
        }
        
    }
    else if ( internals_only == false )
    {
        // get the age of the tip if internalsOnly==false
        ages[n.getIndex()] = n.getAge();
    }

}


/**
 * Helper function to get the MRCA age, recursively
 * @param node current node
 * @param pathFromOtherNodeToRoot path between other node and root of the tree
 * @return MRCA age
 */
double RevBayesCore::TreeUtilities::getAgeOfMRCARecursive(const TopologyNode& node, boost::unordered_set <const TopologyNode* >& node_root_path)
{

    if ( node.isRoot() || node_root_path.find(&node) != node_root_path.end() )
    {
        return node.getAge();
    }
    else
    {
        return getAgeOfMRCARecursive( node.getParent(), node_root_path );
    }
    
}


/**
 * Get age of the MRCA of two tips given by names
 * @param t input tree
 * @param first,second tip names
 * @return MRCA age
 */
double RevBayesCore::TreeUtilities::getAgeOfMRCA(const Tree& t, const std::string& first, const std::string& second)
{

    const TopologyNode& node1 = t.getTipNodeWithName(first) ;

    const TopologyNode& node2 = t.getTipNodeWithName(second) ;

    if (! (node2.equals( node1 ) ) )
    {
        boost::unordered_set <const TopologyNode* > pathFromNode1ToRoot;
        climbUpTheTree(node1, pathFromNode1ToRoot);

        double age = getAgeOfMRCARecursive(node2, pathFromNode1ToRoot);
        return age;
    }
    else
    {
        return node1.getAge();
    }

}


/**
 * Calculate average distance matrix from a vector of distance matrices
 * @param matvect vector of distance matrices
 * @param weights vector of weights, set to NULL by default. If no weights are provided, a vector of equal weights of the same length as matvect is created internally.
 * @return average distance matrix
 *
 */
RevBayesCore::AverageDistanceMatrix RevBayesCore::TreeUtilities::getAverageDistanceMatrix(const RbVector<RevBayesCore::DistanceMatrix>& matvect, const RbVector<double>* weights = NULL)
{
    // test if a vector of weights was provided, and if not, set it to a vector of 1s
    if (weights == NULL)
    {
        weights = new RbVector<double>(matvect.size(), 1.0);
    }
    
    // note that the check that matvect.size() == weights.size() is performed upstream within AvgDistanceMatrixFunction()
    
    // gather all taxa across all source matrices into a single vector
    std::vector<RevBayesCore::Taxon> allTaxa;

    for (RbConstIterator<DistanceMatrix> it = matvect.begin(); it != matvect.end(); ++it)
    {
        allTaxa.insert(allTaxa.end(), it->getTaxa().begin(), it->getTaxa().end());
    }

    // convert the vector of taxa into a vector of strings for easier sorting
    std::vector<std::string> allNames( allTaxa.size() );

    for (size_t i = 0; i < allTaxa.size(); i++)
    {
        allNames[i] = allTaxa[i].getName();
    }

    // get rid of duplicates by converting from vector to unordered_set
    boost::unordered_set<std::string> uniqueNames;

    for (size_t j = 0; j < allNames.size(); j++)
    {
        uniqueNames.insert(allNames[j]);
    }

    // repopulate the original vector with unique values only
    allNames.assign( uniqueNames.begin(), uniqueNames.end() );
    
    // sort alphabetically
    std::sort( allNames.begin(), allNames.end() );

    // initialize the sum and divisor matrices using the size-based constructor:
    RevBayesCore::MatrixReal sumMatrix = MatrixReal( allNames.size() );
    RevBayesCore::MatrixReal divisorMatrix = MatrixReal( allNames.size() );

    // initialize the corresponding Boolean matrix of the right dimensions, filled with 'false'
    RevBayesCore::MatrixBoolean mask = MatrixBoolean( allNames.size() );

    for (size_t mat = 0; mat != matvect.size(); ++mat)
    {
        std::vector<Taxon> taxa = matvect[mat].getTaxa();
        std::vector<size_t> matInd( taxa.size() );
        
        for (size_t i = 0; i != taxa.size(); i++)
        {
            matInd[i] = std::distance( allNames.begin(), std::find( allNames.begin(), allNames.end(), taxa[i].getName() ) );
        }
        
        for (size_t j = 0; j != taxa.size(); j++)
        {
            for (size_t k = j + 1; k != taxa.size(); ++k)
            {
                sumMatrix[ matInd[j] ][ matInd[k] ] += matvect[mat].getMatrix()[j][k] * (*weights)[mat];
                sumMatrix[ matInd[k] ][ matInd[j] ] += matvect[mat].getMatrix()[j][k] * (*weights)[mat]; // by symmetry
                divisorMatrix[ matInd[j] ][ matInd[k] ] += (*weights)[mat];
                divisorMatrix[ matInd[k] ][ matInd[j] ] += (*weights)[mat];                              // by symmetry
                mask[ matInd[j] ][ matInd[k] ] = true;
                mask[ matInd[k] ][ matInd[j] ] = true;                                                   // by symmetry
            }
        }
    }

    // divide the sum matrix by the divisor matrix
    RevBayesCore::MatrixReal averageMatrix = MatrixReal( allNames.size() );

    for (size_t i = 0; i != averageMatrix.getNumberOfRows(); i++)
    {
        for (size_t j = 0; j != averageMatrix.getNumberOfColumns(); j++)
        {
            if (divisorMatrix[i][j] > 0.0)
            {
                averageMatrix[i][j] = sumMatrix[i][j]/divisorMatrix[i][j];
            }
        }
    }

    // convert from a vector of strings back to a vector of taxa:
    std::vector<Taxon> uniqueTaxa( allNames.size() );

    for(size_t k = 0; k != allNames.size(); k++)
    {
        uniqueTaxa[k] = Taxon( allNames[k] );
    }

    // initialize a distance matrix based on the average matrix and taxon vector obtained above:
    DistanceMatrix dm = DistanceMatrix( averageMatrix, uniqueTaxa );

    // combine with the Boolean mask into an average distance matrix:
    AverageDistanceMatrix adm = AverageDistanceMatrix( dm, mask );

    return adm;
}



/**
 * Calculate the Colless metric for a specific node
 * @param node current node
 * @param size current size
 * @return metric for the node
 */
int RevBayesCore::TreeUtilities::getCollessMetric(const TopologyNode& node, int& size)
{
    if ( node.isTip() )
    {
        size = (node.getAge() == 0.0);
        return 0.0;
    }

    const TopologyNode& left  = node.getChild(0);
    const TopologyNode& right = node.getChild(1);

    int left_size  = 0;
    int right_size = 0;

    double left_metric  = getCollessMetric(left, left_size);
    double right_metric = getCollessMetric(right, right_size);

    size = left_size + right_size;

    int metric = std::abs( left_size - right_size);

    if ( left_size == 0 || right_size == 0 )
    {
        metric = 0;
    }

    return left_metric + right_metric + metric;
}


/**
 * Calculate distance matrix based on branch lengths between all tips of a tree
 * @param tree input tree
 * @return distance matrix
 */
RevBayesCore::DistanceMatrix* RevBayesCore::TreeUtilities::getDistanceMatrix(const Tree& tree)
{

    RevBayesCore::MatrixReal* matrix = new MatrixReal( tree.getNumberOfTips() );

    std::vector<Taxon> names = tree.getTaxa( ) ;

    std::map< std::string, int > namesToId;

    for (size_t i = 0; i < names.size(); ++i)
    {
        namesToId[ names[i].getName() ] = int(i);
    }

    std::vector< std::pair<std::string, double> > distsToRoot;

    processDistsInSubtree( tree.getRoot() , *matrix, distsToRoot, namesToId);

    DistanceMatrix* distMat = new DistanceMatrix(*matrix, names);

    // free memory
    delete matrix;

    return distMat;
}


/**
 * Calculate the parsimony score of a tree and alignment based on the algorithm from Fitch (1970) "Distinguishing Homologous from Analogous Proteins".
 * @param t input tree
 * @param c input character state alignment
 * @return parsimony score
 */
int RevBayesCore::TreeUtilities::getFitchScore(const Tree& t, const AbstractHomologousDiscreteCharacterData& c)
{
    int score = 0;
    for (size_t i = 0; i < c.getNumberOfCharacters(); i++)
    {
        recursivelyComputeFitch(t.getRoot(), c, i, score);
    }
    
    return score;
}


/**
 * Calculate the Gamma statistic from Pybus & Harvey (2000) equation 1
 * @param t input tree
 * @return gamma statistic
 */
double RevBayesCore::TreeUtilities::getGammaStatistic(const Tree& tree)
{
    std::vector<TopologyNode*> nodes = tree.getNodes();

    std::vector<double> ages;
    for (size_t i = 0; i < nodes.size(); i++)
    {
        ages.push_back(nodes[i]->getAge());
    }

    // calculate internode distances
    std::sort(ages.begin(), ages.end());
    std::vector<double> distances;
    for (size_t i = (ages.size() - 1); i > 0; i--)
    {
        distances.push_back(ages[i] - ages[i - 1]);
        if (ages[i - 1] == 0)
        {
            break;
        }
    }

    double n = tree.getNumberOfTips();
    if (n < 3)
    {
        //return NaN;
        return std::numeric_limits<double>::quiet_NaN();
    }

    double T = 0;
    for (int j = 2; j <= n; j++)
    {
        T = T + (j * distances[j - 2]);
    }

    double a = 1 / ( n - 2 );
    double b = 0;
    for (int i = 2; i <= (n - 1); i++)
    {
        double temp = 0;
        for (int k = 2; k <= i; k++)
        {
            temp = temp + (k * distances[k - 2]);
        }
        b = b + temp;
    }
    double num = (a * b) - (T / 2);

    double den = T * sqrt( (1 / (12 * (n - 2))) );


    return num/den;
}


/**
 * Calculate the inverse equal splits metric for all tips.
 * @param t input tree
 * @return vector of inverse ES for tips
 */
std::vector<double> RevBayesCore::TreeUtilities::getInverseES(const Tree& tree)
{
    if ( tree.isRooted() == false )
    {
        throw RbException("Inverse ES can only be calculated on rooted trees.");
    }

    std::vector<double> inverse_es;

    // calculate equal splits (ES) measure for each tip
    for (size_t i = 0; i < tree.getNumberOfTips(); i++)
    {
        double tip_es = 0;
        size_t node_index = tree.getTipNode(i).getIndex();

        // traverse from tip to root
        double depth = 1;
        while (true)
        {
            if (tree.getNode(node_index).isRoot() == true)
            {
                break;
            }
            tip_es += tree.getNode(node_index).getBranchLength() * (1 / pow(2, depth - 1));
            node_index = tree.getNode(node_index).getParent().getIndex();
            depth++;
        }
        if (tip_es != 0)
        {
            inverse_es.push_back(1/tip_es);
        }
    }

    return inverse_es;
}


/**
 * Calculate the mean inverse equal splits metric for tips in a single state as described in:
 * Rabosky and Goldberg (2017) "FiSSE: A simple nonparametric test for the effects of a binary character on lineage diversiﬁcation rates"
 * This metric is typically calculated for a single character at a time, but here it is extended over multiple characters.
 * @param t input tree
 * @param c input character state alignment
 * @param state_index index of the specific state
 * @return mean inverse ES
 */
double RevBayesCore::TreeUtilities::getMeanInverseES(const Tree& tree, const AbstractHomologousDiscreteCharacterData& c, size_t state_index)
{
    if ( tree.isRooted() == false )
    {
        throw RbException("Mean inverse ES can only be calculated on rooted trees.");
    }

    std::vector<double> summed_inverse_es = std::vector<double>(c.getNumberOfCharacters(), 0);
    std::vector<double> num_tips_in_state = std::vector<double>(c.getNumberOfCharacters(), 0);
    std::vector<std::string> tip_names = tree.getTipNames();

    // calculate equal splits (ES) measure for each tip as necessary
    for (size_t i = 0; i < tip_names.size(); i++)
    {
        bool calculated_for_tip = false;
        double tip_es = 0;
        size_t node_index = tree.getTipNodeWithName( tip_names[i] ).getIndex();

        for (size_t j = 0; j < c.getNumberOfCharacters(); j++)
        {
            size_t state = c.getTaxonData(tip_names[i]).getCharacter(j).getStateIndex();
            if (state == state_index)
            {
                num_tips_in_state[j] += 1;
                if (calculated_for_tip == true)
                {
                    if (tip_es != 0)
                    {
                        summed_inverse_es[j] += 1/tip_es;
                    }
                }
                else
                {
                    // traverse from tip to root
                    double depth = 1;
                    while (true)
                    {
                        if (tree.getNode(node_index).isRoot() == true)
                        {
                            break;
                        }
                        tip_es += tree.getNode(node_index).getBranchLength() * (1 / pow(2, depth - 1));
                        node_index = tree.getNode(node_index).getParent().getIndex();
                        depth++;
                    }
                    if (tip_es != 0)
                    {
                        summed_inverse_es[j] += 1/tip_es;
                    }
                    calculated_for_tip = true;
                }
            }
        }
    }

    // calculate mean inverse ES for the character state
    double mean_inverse_es = 0;
    for (size_t i = 0; i < c.getNumberOfCharacters(); i++)
    {
        if (num_tips_in_state[i] != 0)
        {
            mean_inverse_es += (1/num_tips_in_state[i]) * summed_inverse_es[i];
        }
    }
    return mean_inverse_es;
}


/**
 * Get node index of the MRCA of two nodes
 * @param left,right input nodes
 * @return index of MRCA
 */
size_t RevBayesCore::TreeUtilities::getMrcaIndex(const TopologyNode& left, const TopologyNode& right)
{

    if ( &left == &right )  //same
    {
        return left.getIndex();
    }
    else if ( left.getAge() < right.getAge() )
    {
        return RevBayesCore::TreeUtilities::getMrcaIndex( left.getParent(), right );
    }
    else
    {
        return RevBayesCore::TreeUtilities::getMrcaIndex( left, right.getParent() );
    }

}


/**
 * Get nodal distance between two nodes. Nodal distance is 1 for each node on the shortest path between the two.
 * @param left,right nodes between which to calculate the distance
 * @return nodal distance
 */
int RevBayesCore::TreeUtilities::getNodalDistance(const TopologyNode& left, const TopologyNode& right)
{
    if ( &left == &right || &left.getParent() == &right || &left == &right.getParent() )
    {
        return 0;
    }
    else if ( left.getAge() < right.getAge() )
    {
        return 1 + RevBayesCore::TreeUtilities::getNodalDistance( left.getParent(), right );
    }
    else
    {
        return 1 + RevBayesCore::TreeUtilities::getNodalDistance( left, right.getParent() );
    }
}

/**
 * Get matrix of nodal distances between all tips of a tree
 * @param tree input tree
 * @return distance matrix indexed by tip indexes
 */
RevBayesCore::DistanceMatrix* RevBayesCore::TreeUtilities::getNodalDistanceMatrix(const Tree& tree)
{
    RevBayesCore::MatrixReal matrix = MatrixReal( tree.getNumberOfTips() );

    std::vector<Taxon> names = tree.getTaxa( ) ;
    for (size_t i = 0; i < names.size(); i++)
    {
        for (size_t j = i + 1; j < names.size(); j++)
        {
            matrix[i][j] = matrix[j][i] = TreeUtilities::getNodalDistance(tree.getTipNode(i), tree.getTipNode(j));
        }
    }

    return new DistanceMatrix(matrix, names);
}

/**
 * Get the age of the oldest tip below a specific node of a binary tree
 * @param[in] t binary tree
 * @param[in] n top node
 * @param[out] oldest current oldest age found
 *
 * @todo for some reason, this function only works for binary trees
 */
double RevBayesCore::TreeUtilities::getOldestTipAge(const TopologyNode& n)
{

    if ( n.isTip() == false )
    {

        // assertion that we have binary trees
#ifdef ASSERTIONS_TREE
        if ( n->getNumberOfChildren() != 2 )
        {
            throw RbException("Oldest tip is only implemented for binary trees!");
        }
#endif

        double left_age  = getOldestTipAge( n.getChild(0) );
        double right_age = getOldestTipAge( n.getChild(1) );
        
        return (left_age > right_age ? left_age : right_age);
    }
    else
    {
        return n.getAge();
    }
}


/**
 * Calculate the Parsimoniously Same State Paths (PSSP). This is the set of branch lengths
 * from clades parsimoniously reconstructed to have the same state. Given (((A,B),C),(D,E)), if A, B,
 * D, and E are in state 0, then PSSP(0) will contain the four branch lengths in (A,B) and (D,E). Uses
 * Fitch's (1970) algorithm for parsimony ancestral state reconstruction.
 * @param t input tree
 * @param c input character state alignment
 * @param state_index index of the specific state
 * @return vector of PSSP
 */
std::vector<double> RevBayesCore::TreeUtilities::getPSSP(const Tree& tree, const AbstractHomologousDiscreteCharacterData& c, size_t state_index)
{
    std::vector<double> branch_lengths;
    if ( c.getNumberOfCharacters() != 1 )
    {
        throw RbException("getPSSP is only implemented for character alignments with a single site.");
    }
    recursivelyGetPSSP(tree.getRoot(), c, branch_lengths, state_index);
    
    return branch_lengths;
}

/**
 * Get all tips below specified node, recursively
 * @param n current node
 * @param[out] taxa vector of tips found so far
 */
void RevBayesCore::TreeUtilities::getTaxaInSubtree(TopologyNode& n, std::vector<TopologyNode*> &taxa )
{

    if ( n.isTip() )
    {
        taxa.push_back( &n );
    }
    else
    {

        // recursively add children to the list of nodes in this subtree
        for (size_t i = 0; i < n.getNumberOfChildren(); ++i)
        {
            TopologyNode& child = n.getChild( i );

            getTaxaInSubtree( child, taxa );
        }
    }

}

/**
 * Check if the two trees are connected by a single NNI move
 */
bool RevBayesCore::TreeUtilities::isConnectedNNI(const Tree& a, const Tree& b)
{
    
    size_t num_nodes = a.getNumberOfNodes();
    
//    // now exchange the two nodes
//    parent.removeChild( node_B );
//    node->removeChild( node_A );
//    parent.addChild( node_A );
//    node->addChild( node_B );
//    node_A->setParent( &parent );
//    node_B->setParent( node );
    
    return false;
}



/**
 * Make tree ultrametric by extending terminal branches to the present
 * @param t tree to be modified
 */
void RevBayesCore::TreeUtilities::makeUltrametric(Tree& tree)
{

    // make sure that the tree is currently using ages
    tree.getRoot().setUseAges(true, true);
    
    double max = 0.0;
    std::vector<double> ages ;
    for (size_t i = 0; i < tree.getNumberOfTips(); ++i)
    {
        TopologyNode* node = &(tree.getTipNode( i ) );
        double age = node->getBranchLength();
        node = &(node->getParent());

        while ( node->isRoot() == false )
        {
            age += node->getBranchLength();
            node = &(node->getParent());
        }
        if (age > max)
        {
            max = age;
        }
        ages.push_back(age);

    }

    // We extend terminal branches
    for (size_t i = 0; i < tree.getNumberOfTips(); ++i)
    {
        tree.getTipNode( i ).setBranchLength(tree.getTipNode( i ).getBranchLength() + max - ages[i]);
    }
    
    // finally, make sure that all the internal nodes have the ages properly set
    tree.getRoot().recomputeAge(true);


}


void RevBayesCore::TreeUtilities::minBLTimeScaling(Tree& treeToScale, const std::vector<Taxon>& taxa, const double minBrLen)
{
    // Check that the user-supplied tree contains the same number of tips as the vector of taxa
    size_t tip_num = treeToScale.getNumberOfTips();
    size_t tax_num = taxa.size();
    
    if (tip_num != tax_num)
    {
        throw RbException("Number of tips in the initial tree does not match the number of taxa.");
    }
    
    // Check that the tip labels of the user-supplied tree match those of the vector of taxa
    std::vector<std::string> tip_names;
    for (size_t i = 0; i < tip_num; ++i)
    {
        const TopologyNode& n = treeToScale.getTipNode( i );
        tip_names.push_back( n.getTaxon().getName() );
    }

    std::vector<std::string> taxon_names;
    for (size_t i = 0; i < tax_num; ++i)
    {
        taxon_names.push_back( taxa[i].getName() );
    }
    
    std::sort(tip_names.begin(), tip_names.end());
    std::sort(taxon_names.begin(), taxon_names.end());
    if (tip_names != taxon_names)
    {
        throw RbException("Tip names of the initial tree do not match the taxon names.");
    }
    
    // Alter the tip age values of p in place
    for (size_t i = 0; i < tax_num; ++i)
    {
        std::string tip_name = taxa[i].getName();
        treeToScale.setTaxonObject( tip_name, taxa[i] );
    }
    
    // Grab the first tip; which one we start with should not matter
    TopologyNode& node = treeToScale.getTipNode( 0 );
    node.setParentAge( minBrLen );
}


/**
 * Offset the age of a node and its children by a factor
 * @param t tree to be modified
 * @param n top node to offset
 * @param factor offset value
 */
void RevBayesCore::TreeUtilities::offsetTree(TopologyNode& node, double factor)
{
    // rescale the time of the node
    double new_age = node.getAge() + factor;
    node.setAge( new_age );

    // offset all children
    const std::vector<TopologyNode*>& children = node.getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        offsetTree( *children[i], factor);
    }

}


/**
 * Helper function for calculating distance matrix between all tips of a tree, recursively
 * @param[in] node current node
 * @param[out] matrix current distance matrix
 * @param[out] distsToNodeFather vector to store distance between leaves and node parent
 * @param[in] namesToId map of node names to their index
 */
void RevBayesCore::TreeUtilities::processDistsInSubtree(const TopologyNode& node, MatrixReal& matrix, std::vector< std::pair<std::string, double> >& distsToNodeFather, const std::map< std::string, int >& namesToId)
{
    distsToNodeFather.clear();

    // node-is-leaf case
    if ( node.isTip() == true )
    {
        distsToNodeFather.push_back(make_pair ( node.getName(), node.getBranchLength() ) );
        return;
    }

    // For all leaves in node's subtree, get leaf-to-node distances.
    // Leaves are ordered according to the children of node.
    std::map< size_t, std::vector< std::pair<std::string, double> > > leavesDists;
    for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
    {
        const TopologyNode* child = &( node.getChild(i) );
        processDistsInSubtree(*child, matrix, leavesDists[child->getIndex()], namesToId); // recursivity
    }

    // Write leaf-leaf distances to the distance matrix.
    // Only pairs in which the two leaves belong to different
    // children are considered.
    for (size_t son1_loc = 0; son1_loc < node.getNumberOfChildren(); ++son1_loc)
    {
        for (size_t son2_loc = 0; son2_loc < son1_loc; ++son2_loc)
        {
            const TopologyNode* son1 = &(node.getChild(son1_loc) );
            const TopologyNode* son2 = &(node.getChild(son2_loc) );

            for (std::vector< std::pair<std::string, double> >::iterator son1_leaf = leavesDists[son1->getIndex()].begin();
                 son1_leaf != leavesDists[son1->getIndex()].end();
                 ++son1_leaf)
            {
                for (std::vector< std::pair<std::string, double> >::iterator son2_leaf = leavesDists[son2->getIndex()].begin();
                     son2_leaf != leavesDists[son2->getIndex()].end();
                     ++son2_leaf)
                {
                    int son1_leaf_id = namesToId.at( son1_leaf->first );
                    int son2_leaf_id = namesToId.at( son2_leaf->first );
                    matrix[son1_leaf_id] [son2_leaf_id] =
                    matrix[son2_leaf_id] [son1_leaf_id] =
                    ( son1_leaf->second + son2_leaf->second );
                }
            }
        }
    }

    // node-is-root case
    if (node.isRoot())
    {
        // node-is-root-and-leaf case
        if (node.isTip() )
        {
            std::string root_name = node.getName();
            for (std::vector< std::pair<std::string, double> >::iterator other_leaf = leavesDists[node.getChild(0).getIndex()].begin();
                 other_leaf != leavesDists[node.getChild(0).getIndex()].end();
                 ++other_leaf)
            {
                matrix [ namesToId.at(root_name) ] [ namesToId.at(other_leaf->first) ]= matrix[ namesToId.at(other_leaf->first) ] [ namesToId.at(root_name) ] = other_leaf->second;
            }
        }

        return;
    }


    // Get distances from node's parent to considered leaves
    distsToNodeFather.clear();
    double nodeToFather = node.getBranchLength();

    for (std::map<size_t, std::vector<std::pair<std::string, double> > >::iterator son = leavesDists.begin(); son != leavesDists.end(); ++son)
    {
        for (std::vector< std::pair<std::string, double> >::iterator leaf = (son->second).begin(); leaf != (son->second).end(); ++leaf)
        {
            distsToNodeFather.push_back(make_pair(leaf->first, (leaf->second + nodeToFather)));
        }
    }

}


/**
 * Helper function for the parsimony score calculation
 * @param node current node
 * @param c input character state alignment
 * @param site current site index
 * @param score current parsimony score
 * @return set of parsimonious characters at the node
 */
std::set<size_t> TreeUtilities::recursivelyComputeFitch(const TopologyNode& node, const AbstractHomologousDiscreteCharacterData &c, size_t site, int &score)
{
    if ( node.isTip() == true )
    {
        std::set<size_t> tip_set;
        const std::string& name = node.getName();
        size_t state = c.getTaxonData(name).getCharacter(site).getStateIndex();
        tip_set.insert( state );
        return tip_set;
    }
    else
    {
        if ( node.getNumberOfChildren() != 2 )
        {
            throw RbException("Fitch score calculation is only implemented for binary trees.");
        }
        std::set<size_t> l = recursivelyComputeFitch(node.getChild(0), c, site, score);
        std::set<size_t> r = recursivelyComputeFitch(node.getChild(1), c, site, score);

        std::set<size_t> intersect;
        set_intersection(l.begin(), l.end(), r.begin(), r.end(), std::inserter(intersect, intersect.begin()));

        if ( intersect.size() == 0 )
        {
            score++;
            std::set<size_t> union_set;
            set_union(l.begin(), l.end(), r.begin(), r.end(), std::inserter(union_set, union_set.begin()));
            return union_set;
        }
        return intersect;
    }
}


/** Helper function for PSSP calculation
 * @param current node
 * @param c input character state alignment
 * @param branch_lengths PSSP calculated so far
 * @param state_index index of the specific state
 * @return set of branches in state given by state_index
 */
std::set<size_t> TreeUtilities::recursivelyGetPSSP(const TopologyNode& node, const AbstractHomologousDiscreteCharacterData& c, std::vector<double> &branch_lengths, size_t state_index)
{
    if ( node.isTip() == true )
    {
        std::set<size_t> tip_set;
        const std::string& name = node.getName();
        size_t state = c.getTaxonData(name).getCharacter(0).getStateIndex();
        tip_set.insert( state );
        return tip_set;
    }
    else
    {
        if ( node.getNumberOfChildren() != 2 )
        {
            throw RbException("getPSSP is only implemented for binary trees.");
        }
        std::set<size_t> l = recursivelyGetPSSP(node.getChild(0), c, branch_lengths, state_index);
        std::set<size_t> r = recursivelyGetPSSP(node.getChild(1), c, branch_lengths, state_index);

        std::set<size_t> intersect;
        set_intersection(l.begin(), l.end(), r.begin(), r.end(), std::inserter(intersect, intersect.begin()));

        if (intersect.size() == 0)
        {
            std::set<size_t> union_set;
            set_union(l.begin(), l.end(), r.begin(), r.end(), std::inserter(union_set, union_set.begin()));
            return union_set;
        }
        if (intersect.size() == 1)
        {
            if (intersect.find(state_index) != intersect.end())
            {
                branch_lengths.push_back(node.getChild(0).getBranchLength());
                branch_lengths.push_back(node.getChild(1).getBranchLength());
            }
        }
        return intersect;
    }
}


/**
* Rescale subtree by scaling the age of an internal node and all its children except tips by a factor
* @param t tree to be modified
* @param n top node to rescale
* @param factor factor used for rescaling
* @param verbose whether to explicitly check for the tree being binary
*
* @note the difference with rescaleTree is that rescaleSubtree will do nothing on tips
*/
void TreeUtilities::rescaleSubtree(TopologyNode& node, double factor, bool verbose)
{
    // we only rescale internal nodes which have no SA as children
    if ( !node.isTip() && !node.isSampledAncestorParent())
    {
        // rescale the age of the node
        double new_age = node.getAge() * factor;
        node.setAge(new_age);

        // assertion that we have binary trees
        if ( verbose == true )
        {
            if ( node.getNumberOfChildren() != 2 )
            {
                throw RbException("Subtree scaling is only implemented for binary trees!");
            }
        }

        // rescale both children
        rescaleSubtree( node.getChild(0), factor, verbose);
        rescaleSubtree( node.getChild(1), factor, verbose);
    }

}


/**
* Rescale tree by scaling the age of a node and all its children by a factor
* @param t tree to be modified
* @param n top node to rescale
* @param factor factor used for rescaling
*/
void RevBayesCore::TreeUtilities::rescaleTree(TopologyNode& node, double factor)
{
    // recursive call for internal nodes
    if ( node.isTip() == false )
    {
        // we only rescale non-tip nodes
        // rescale the time of the node
        double new_age = node.getAge() * factor;
        node.setAge( new_age);


        // assertion that we have binary trees
#ifdef ASSERTIONS_TREE
        if ( n->getNumberOfChildren() != 2 )
        {
            throw RbException("Tree scaling is only implemented for binary trees!");
        }
#endif

        // rescale both children
        rescaleTree( node.getChild(0), factor);
        rescaleTree( node.getChild(1), factor);
    }

}


/**
 * Set ages of a node and its children to values given by a vector
 * @param t tree to be modified
 * @param n top node
 * @param ages new age values
 */
void RevBayesCore::TreeUtilities::setAges(TopologyNode& node, const std::vector<double>& ages)
{
    // we only rescale internal nodes
    if ( node.isTip() == false )
    {
        // rescale the age of the node
        node.setAge( ages[node.getIndex()] );

        // rescale both children
        const std::vector<TopologyNode*>& children = node.getChildren();
        for (size_t i = 0; i < children.size(); i++)
        {
            setAges( *children[i], ages);
        }

    }

}

/**
 * Set age of node then rescale all children to keep the branch lengths constant
 * @param t tree to be modified
 * @param n node whose age is set
 * @param age new node age
 */
void RevBayesCore::TreeUtilities::setAgesRecursively(TopologyNode& node, double age)
{
    // first, we set the age of this node
    node.setAge( age );

    // we only rescale internal nodes
    if ( node.isTip() == false )
    {

        // rescale both children
        const std::vector<TopologyNode*>& children = node.getChildren();
        for (size_t i = 0; i < children.size(); ++i)
        {
            setAgesRecursively( *(children[i]), age-children[i]->getBranchLength());
        }

    }

}



/**
 * Make sure that the starting tree used to initialize serial-sampling birth-death analyses satisfies specified age constraints, and assign min/max ages to its tips
 * @param treeToChange user-supplied starting tree to be modified
 * @param taxaToCopy vector of Taxon objects corresponding to the tips of the tree
 * @param agePrecision how many decimal places to use when checking for compatibility between the tip ages from treeToChange and taxaToCopy
 */
Tree* RevBayesCore::TreeUtilities::startingTreeInitializer(Tree& treeToChange, std::vector<Taxon>& taxaToCopy, long agePrecision)
{
    // Check that the user-supplied tree contains the same number of tips as the vector of taxa
    size_t tip_num = treeToChange.getNumberOfTips();
    size_t tax_num = taxaToCopy.size();
    
    if (tip_num != tax_num)
    {
        throw RbException("Number of tips in the initial tree does not match the number of taxa.");
    }
    
    // Check that the tip labels of the user-supplied tree match those of the vector of taxa
    std::vector<std::string> tip_names;
    for (size_t i = 0; i < tip_num; ++i)
    {
        const TopologyNode& n = treeToChange.getTipNode( i );
        tip_names.push_back( n.getTaxon().getName() );
    }

    std::vector<std::string> taxon_names;
    for (size_t i = 0; i < tax_num; ++i)
    {
        taxon_names.push_back( taxaToCopy[i].getName() );
    }
    
    std::sort(tip_names.begin(), tip_names.end());
    std::sort(taxon_names.begin(), taxon_names.end());
    if (tip_names != taxon_names)
    {
        throw RbException("Tip names of the initial tree do not match the taxon names.");
    }
    
    double factor = pow(10.0, agePrecision);
    
    // Check that the ages of the taxa in the user-supplied tree fall within their age uncertainty ranges
    for (size_t i = 0; i < tax_num; ++i)
    {
        // Grab the minimum and maximum ages of the i-th taxon from the taxaToCopy vector
        double min_age = std::round(taxaToCopy[i].getMinAge() * factor) / factor;
        double max_age = std::round(taxaToCopy[i].getMaxAge() * factor) / factor;
        // Get the age of the tree tip of the same name
        double tip_age = std::round(treeToChange.getTipNodeWithName( taxaToCopy[i].getName() ).getAge() * factor) / factor;

        if (tip_age < min_age || tip_age > max_age)
        {
            std::ostringstream s;
            s << "The age of tip '" << taxaToCopy[i].getSpeciesName() << "' in the initial tree is outside of specified bounds.";
            throw RbException( s.str() );
        }
    }
    
    // Alter the tip age values of treeToChange in place
    for (size_t i = 0; i < tax_num; ++i)
    {
        std::string tip_name = taxaToCopy[i].getName();
        treeToChange.setTaxonObject( tip_name, taxaToCopy[i] );
    }
    
    RevBayesCore::Tree *p = &treeToChange;
    return p;
}
