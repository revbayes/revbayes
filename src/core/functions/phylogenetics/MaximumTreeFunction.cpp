#include "MaximumTreeFunction.h"

#include <limits>
#include <set>

#include "RbException.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeUtilities.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"
#include "boost/dynamic_bitset.hpp" // IWYU pragma: keep

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

MaximumTreeFunction::MaximumTreeFunction( const TypedDagNode< RbVector<Tree> > *ts ) : TypedFunction<Tree>( new Tree() ),
    trees( ts )
{
    // add the lambda parameter as a parent
    addParameter( ts );

    update();
}


MaximumTreeFunction::~MaximumTreeFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



MaximumTreeFunction* MaximumTreeFunction::clone( void ) const
{
    return new MaximumTreeFunction( *this );
}


void MaximumTreeFunction::update( void )
{
    getMinDepthMatrix (  );
    value = getSpeciesTreeFromMinDepths()->clone();
}



void MaximumTreeFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

        trees = static_cast<const TypedDagNode< RbVector<Tree> >* >( newP );

}

/**---------------------------------------------------------------------
 |   Code taken from mrBayes 3.2.2.
 |   GetMinDepthMatrix: converted from GetMinDists.
 |
 |   This algorithm scans the gene trees and calculates the minimum depth
 |   (height) separating species across gene trees. The complexity of the
 |   original algorithm was O(mn^3), where m is the number of gene trees and
 |   n is the number of taxa in each gene tree. I think this algorithm has
 |   complexity that is better on average, but the difference is small.
 |
 |   I have rewritten the algorithm also to show alternative techniques that
 |   could be used in this and other BEST algorithms.
 |
 ----------------------------------------------------------------------*/

void MaximumTreeFunction::getMinDepthMatrix (  )
{
    depthMatrix.clear();

    // First we get the list of all species present in the gene trees
    speciesNamesV = trees->getValue()[0].getSpeciesNames();
    double maxDepth       = trees->getValue()[0].getNode( trees->getValue()[0].getRoot().getIndex() ).getAge();
    std::set<std::string> speciesNames ;
    for (size_t i = 0; i < speciesNamesV.size(); ++i )
    {
        speciesNames.insert ( speciesNamesV[i] );
    }
    if (speciesNames.size() <= 1  ) {
        throw RbException( "Not enough species in the first gene tree. Are you sure you provided species information?" );
    }

    std::set<std::string>::iterator it;
    double age;
    for (size_t i = 1; i < trees->getValue().size(); ++i )
    {
        std::vector<std::string> s = trees->getValue()[i].getSpeciesNames();
        for (size_t j = 0; j < s.size(); ++j )
        {
            if ( speciesNames.find(s[j])==speciesNames.end() )
            {
                speciesNames.insert(s[j]);
            }
             age = trees->getValue()[0].getNode( trees->getValue()[0].getRoot().getIndex() ).getAge();
            if ( age > maxDepth )
            {
                maxDepth = age;
            }

        }
    }

    numSpecies = speciesNames.size();

    speciesToIndex.clear();
    speciesNamesV.clear();
    //Now we have all species names.
    boost::dynamic_bitset<> speciesAbsent = boost::dynamic_bitset<> ( numSpecies ) ;
    //All leaves are present in completeBitVector:
    size_t index = 0;

    for (it = speciesNames.begin(); it != speciesNames.end() ; ++it)
    {
        speciesAbsent[index] = 0;
        speciesNamesV.push_back( *it );
        std::pair<std::string, size_t> pair = std::pair<std::string, size_t> (*it, index);
        speciesToIndex.insert( pair );
        index = index + 1;

    }

    // Set initial max depth for upper triangular matrix
    size_t numUpperTriang = (numSpecies * (numSpecies - 1)) / 2;
    for (size_t i=0; i<numUpperTriang; i++)
        depthMatrix.push_back( maxDepth );



    speciesPairSets.clear();
    for (size_t i=0; i<numSpecies-1; i++) {
        for (size_t j=i+1; j<numSpecies; j++) {
            speciesPairSets.push_back( std::pair< std::string, std::string > ( speciesNamesV[i], speciesNamesV[j] ) ) ;
        }
    }


    // Now we are ready to cycle over gene trees
    std::vector < boost::dynamic_bitset<> > allBitsets;
    std::vector < TopologyNode* > nodes;
    std::map < double, boost::dynamic_bitset<> > nodeAgesToSpeciesBitset;
    for (size_t w=0; w < trees->getValue().size(); w++) {
        allBitsets.clear();
        nodeAgesToSpeciesBitset.clear();
        // Set species sets for interior nodes and get their age. O(n)
        nodes = trees->getValue()[w].getNodes();
        for (size_t i=0; i< nodes.size(); i++) {
            if (nodes[i]->isTip()) {
                boost::dynamic_bitset<> nodeBitset =  speciesAbsent;
                nodeBitset[ speciesToIndex[nodes[i]->getSpeciesName()] ] = 1;
                allBitsets.push_back(nodeBitset);
            }
            else { //Interior node, we do the union of the children bitsets
                boost::dynamic_bitset<> nodeBitset = allBitsets[ nodes[i]->getChild(0).getIndex() ];
                for (size_t j=1; j< nodes[i]->getNumberOfChildren() ; j++) {
                    nodeBitset |= allBitsets[ nodes[i]->getChild(j).getIndex() ];
                }
                allBitsets.push_back(nodeBitset);
                double age = nodes[i]->getAge();
                while ( nodeAgesToSpeciesBitset.find(age) != nodeAgesToSpeciesBitset.end() ) {
                    age += std::numeric_limits<double>::min();
                }
                nodeAgesToSpeciesBitset[ age ] = nodeBitset;
            }
        }

        // Now we have bitsets and ages for all nodes in the gene tree
        // nodeAgesToSpeciesBitset is ordered according to node age, which is useful for the next loop.

        // Finally find the minimum for each cell in the upper triangular matrix
        // This is the time critical step with complexity O(n^3) in the simplest
        // algorithm version. This algorithm should do a little better in most cases.
        std::map < double, boost::dynamic_bitset<> >::iterator jt;
        for (size_t i=0; i<numUpperTriang; ++i) {
            jt = nodeAgesToSpeciesBitset.begin();
            while ( jt->first < depthMatrix[i] && jt != nodeAgesToSpeciesBitset.end() ) {
                //If the subtree contains both species...
                if ( jt->second[ speciesToIndex[speciesPairSets[i].first] ] && jt->second[ speciesToIndex[speciesPairSets[i].second] ] ) {
                    depthMatrix[i] = jt->first;
                    break;
                }
                //if not...
                jt++;
            }

         }

    }   // Next gene tree

    return ;
}




/**---------------------------------------------------------------------
 |   Taken from MrBayes 3.2.2
 |   GetSpeciesTreeFromMinDepths: converted from GetConstraints, Startsptree,
 |   and MaximumTree.
 |
 |   This is a clustering algorithm based on minimum depths for species pairs.
 |   It reduces an n choose 2 upper triangular min depth matrix to an array
 |   of n-1 node depths, which fit onto a tree.
 |
 |   @param      speciesTree     The species tree to be filled  (out)
 |   @param      depthMatrix     The min depth matrix, upper triangular array (in)
 |   @returns    Returns NO_ERROR if success, ERROR if negative brlens occur
 ----------------------------------------------------------------------*/
Tree* MaximumTreeFunction::getSpeciesTreeFromMinDepths (  )
{

    //We need to go through depthMatrix from smallest to largest.
    boost::dynamic_bitset<> speciesAbsent = boost::dynamic_bitset<> ( numSpecies ) ;
    for (size_t i = 0; i < numSpecies ; ++i)
    {
        speciesAbsent[i] = 0;
    }
    //std::map< double, std::pair< std::string, std::string >  > depthToPairs ;
    std::map< double, boost::dynamic_bitset< > > depthToPairs ;

    for (size_t i = 0; i< depthMatrix.size(); ++i) {
        boost::dynamic_bitset< > currentPair = speciesAbsent;
        currentPair[ speciesToIndex[speciesPairSets[i].first] ] = 1;
        currentPair[ speciesToIndex[speciesPairSets[i].second] ] = 1;
        depthToPairs[ depthMatrix[i] ] = currentPair;
    }

    //Now we do the loop
    std::map< double, boost::dynamic_bitset< >  >::iterator it;
    std::vector< std::string > speciesDone ;
    boost::dynamic_bitset< > bitSpeciesDone = speciesAbsent;
    for (size_t i = 0; i< numSpecies; ++i) {
        if (depthToPairs.begin()->second[i]) {
            speciesDone.push_back (speciesNamesV[i]);
            bitSpeciesDone[i] = 1;
        }
    }

    double currentBl = depthToPairs.begin()->first;
    it = depthToPairs.begin();
    depthToPairs.erase( it );
    std::map < TopologyNode*, double > nodeToBl;
    std::map < size_t, TopologyNode* > idToNode; //To fill
    //Build the initial cherry tree
    size_t spId = speciesToIndex[speciesDone[0]];
    TopologyNode* node1 = new TopologyNode (speciesDone[0], spId);
    idToNode[spId] = node1;
    spId = speciesToIndex[speciesDone[1]];
    TopologyNode* node2 = new TopologyNode (speciesDone[1], spId);
    idToNode[spId] = node2;
    TopologyNode* Root = new TopologyNode ("");
    node1->setParent ( Root );
    node2->setParent ( Root );
    Root->addChild(node1);
    Root->addChild(node2);
    nodeToBl[node1] = currentBl;
    nodeToBl[node2] = currentBl;
    std::map< double, boost::dynamic_bitset< >  >::iterator currentBestPair;
    double currentMinDepth;
    double totalDepth = currentBl;

    while ( speciesDone.size() < numSpecies )
    {
        currentMinDepth = std::numeric_limits<double>::max();
        bool found = false;
        for ( it = depthToPairs.begin(); it != depthToPairs.end(); ++it )
        {
            //Look for the smallest pair of species for which one of the two species is already included in the current tree.
            if ( bitSpeciesDone.intersects(it->second) && !it->second.is_subset_of( bitSpeciesDone )  && it->first < currentMinDepth) {
                currentMinDepth = it->first;
                currentBestPair = it;
                found=true;
            }

            if ( it->first > currentMinDepth )
            {
                break;
            }
        }

        //We have found the best pair to introduce
        //We get the species that will be added to the current species tree
        //And the other species already in the current species tree as well
        boost::dynamic_bitset< > spToAdd = currentBestPair->second ;
        spToAdd -= bitSpeciesDone;
        boost::dynamic_bitset< > spAlreadyThere = currentBestPair->second ;
        spAlreadyThere -= spToAdd;
        size_t spToAddId = spToAdd.find_first();
        size_t spAlreadyThereId = spAlreadyThere.find_first();
        currentBl = currentBestPair->first;
        TopologyNode* node = new TopologyNode ( speciesNamesV[spToAddId], spToAddId );
        idToNode[spToAddId] = node;
        TopologyNode* nodeR = new TopologyNode ( "" );
        node->setParent( nodeR );
        nodeR->addChild(node);
        nodeToBl[node] = currentBl;
        //Then we need to go through the current tree to know where to branch the new tip on the path from leaf
        //spAlreadyThereId to the current root, or above the current root.
        TopologyNode* futureSibling = idToNode[spAlreadyThereId];
        TopologyNode* futureParent = &( futureSibling->getParent() );
        //Now climb up the tree
        double parentDepth = nodeToBl[futureSibling];
        while (parentDepth < currentBl && !futureParent->isRoot()) {
            futureSibling = futureParent;
            futureParent = &(futureParent->getParent());
            parentDepth += nodeToBl[futureSibling];
        }
        if (parentDepth >= currentBl) {
            double dist = nodeToBl[futureSibling] ;
            double newBl = currentBl - (parentDepth - dist);
            futureParent->removeChild(futureSibling);
            futureSibling->setParent( nodeR );
            nodeR->addChild(futureSibling);
            nodeToBl[futureSibling] = newBl;
            nodeToBl[nodeR] = dist - newBl;
            nodeR->setParent(futureParent);
            futureParent->addChild(nodeR);
        }
        else {
            Root->setParent( nodeR );
            nodeR->addChild(Root);
            nodeToBl[Root] = currentBl - totalDepth;
            totalDepth = currentBl;
            Root = nodeR;
        }

        bitSpeciesDone[spToAddId] = 1;
        speciesDone.push_back( speciesNamesV[spToAddId] );
        depthToPairs.erase( currentBestPair );
    }
    
    //In principle we have just built a proper ultrametric tree.
    //Now we make a tree of it.
    Tree* tree = new Tree();
    tree->setRoot ( Root, true );

    // set the branch lengths
    std::map < TopologyNode*, double >::iterator jt;


    for (jt = nodeToBl.begin() ; jt != nodeToBl.end(); ++jt)
    {
        tree->getNode( jt->first->getIndex() ).setBranchLength( jt->second );
    }

    // set the root age and branch length before calling convertTree
    TopologyNode& root = tree->getRoot();
    root.setAge(totalDepth);
    root.setBranchLength(0.0);

    Tree* ttree = TreeUtilities::convertTree ( *tree ) ;

    // return the tree, the caller is responsible for destruction
    return ttree;




//
//
//    int         i, j, numUpperTriang, nLongsNeeded, index, nextNodeIndex;
//    Depth       *minDepth;
//    PolyTree    *polyTree;
//    PolyNode    *p, *q, *r, *u, *qPrev, *rPrev;
//
//    nLongsNeeded    = ((numSpecies - 1) / nBitsInALong) + 1;
//    numUpperTriang  = numSpecies*(numSpecies - 1) / 2;
//    minDepth        = (Depth *) SafeCalloc (numUpperTriang, sizeof(Depth));
//
//	// Convert depthMatrix to an array of Depth structs
//    index = 0;
//    for (i=0; i<numSpecies; i++) {
//        for (j=i+1; j<numSpecies; j++) {
//            minDepth[index].depth   = depthMatrix[index];
//            minDepth[index].pairSet = speciesPairSets[index];
//            index++;
//        }
//	}
//
//    // Sort the array of distance structs (O(log n^2))
//    qsort((void *)(minDepth), (size_t)(numUpperTriang), sizeof(Depth), CompareDepths);
//
//    // The algorithm below reduces the upper triangular matrix (n choose 2) to an n-1
//    // array in O(n^2log(n)) time. We build the tree at the same time, since we can
//    // find included pairs in the tree in log(n) time. We use a polytomous tree for this.
//
//    // Allocate space for polytomous tree and set up partitions
//    polyTree = AllocatePolyTree(numSpecies);
//    AllocatePolyTreePartitions(polyTree);
//
//    // Build initial tree (a bush)
//    polyTree->isRooted = YES;
//    polyTree->isClock = YES;
//    polyTree->root = &polyTree->nodes[2*numSpecies-2];
//    for (i=0; i<numSpecies; i++) {
//        p = &polyTree->nodes[i];
//        p->index = i;
//        p->depth = 0.0;
//        p->left = NULL;
//        if (i<numSpecies-1)
//            p->sib = &polyTree->nodes[i+1];
//        else
//            p->sib = NULL;
//        p->anc = polyTree->root;
//    }
//    p = polyTree->root;
//    p->index = 2*numSpecies - 2;
//    p->left = &polyTree->nodes[0];
//    p->sib = NULL;
//    p->anc = NULL;
//    p->depth = -1.0;
//    polyTree->nNodes = numSpecies + 1;
//    polyTree->nIntNodes = 1;
//    GetPolyDownPass(polyTree);
//    ResetPolyTreePartitions(polyTree);      /* set bitsets (partitions) for initial tree */
//
//    // Resolve bush using sorted depth structs
//    nextNodeIndex = numSpecies;
//    for (i=0; i<numUpperTriang; i++) {
//
//        // Find tip corresponding to first taxon in pair
//        p = &polyTree->nodes[FirstTaxonInPartition(minDepth[i].pairSet, nLongsNeeded)];
//
//        // Descend tree until we find a node within which the pair set is nested
//        do {
//            p = p->anc;
//        } while (!IsPartNested(minDepth[i].pairSet, p->partition, nLongsNeeded));
//
//        if (p->left->sib->sib != NULL) {
//
//            // This node is still a polytomy
//
//            // Find left and right descendants of new node
//            qPrev = NULL;
//            for (q=p->left; IsSectionEmpty(q->partition, minDepth[i].pairSet, nLongsNeeded); q=q->sib)
//                qPrev = q;
//            rPrev = q;
//            for (r=q->sib;  IsSectionEmpty(r->partition, minDepth[i].pairSet, nLongsNeeded); r=r->sib)
//                rPrev = r;
//
//            // Introduce the new node
//            u = &polyTree->nodes[nextNodeIndex];
//            u->index = nextNodeIndex;
//            nextNodeIndex++;
//            polyTree->nIntNodes++;
//            polyTree->nNodes++;
//            u->left = q;
//            u->anc = p;
//            if (p->left == q)
//                p->left = u;
//            else
//                qPrev->sib = u;
//            // former upstream sibling to r should point to r->sib
//            if (rPrev == q)
//                u->sib = r->sib;
//            else
//                rPrev->sib = r->sib;
//            if (q->sib == r)
//                u->sib = r->sib;
//            else
//                u->sib = q->sib;
//            u->depth = minDepth[i].depth;   // because minDepth structs are sorted, we know this is the min depth
//            assert (u->depth > 0.0);
//
//            // Create new taxon set with bitfield operations
//            for (j=0; j<nLongsNeeded; j++)
//                u->partition[j] = q->partition[j] | r->partition[j];
//
//            // Patch the tree together with the new node added
//            q->sib  = r;
//            r->sib = NULL;
//            q->anc = u;
//            r->anc = u;
//        }
//        else if (p == polyTree->root && p->depth < 0.0) {
//
//            // This is the first time we hit the root of the tree && it is resolved
//            p->depth = minDepth[i].depth;
//            assert (p->depth > 0.0);
//
//        }
//        // other cases should not be added to tree
//    }
//
//    // Make sure we have a complete species tree
//    assert (polyTree->nIntNodes == numSpecies - 1);
//
//    // Set traversal sequences
//    GetPolyDownPass(polyTree);
//
//    // Set branch lengths from node depths (not done automatically for us)
//    // Make sure all branch lengths are nonnegative (we can have 0.0 brlens, they
//    // should not be problematic in a species tree; they occur when there are
//    // ties in the min depth matrix that have not been modified by the move)
//    for (i=0; i<polyTree->nNodes; i++) {
//        p = polyTree->allDownPass[i];
//        if (p->anc == NULL)
//            p->length = 0.0;
//        else
//            p->length = p->anc->depth - p->depth;
//        if (p->length < 0.0 ) {
//            FreePolyTree(polyTree);
//            free (minDepth);
//            return (ERROR);
//        }
//    }
//
//    // Copy to species tree from polytomous tree
//    CopyToSpeciesTreeFromPolyTree (speciesTree, polyTree);
//
//    // Free locally allocated variables
//    FreePolyTree(polyTree);
//    free (minDepth);
//
//    return(NO_ERROR);
}
