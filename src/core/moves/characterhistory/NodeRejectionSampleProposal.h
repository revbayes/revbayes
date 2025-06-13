#ifndef NodeRejectionSampleProposal_H
#define NodeRejectionSampleProposal_H

#include "AbstractRateMatrix.h"
#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "DeterministicNode.h"
#include "HomologousDiscreteCharacterData.h"
#include "DistributionBinomial.h"
#include "DistributionPoisson.h"
#include "PathRejectionSampleProposal.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateGeneratorSequence.h"
#include "RbException.h"
#include "Simplex.h"
#include "StochasticNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TreeChangeEventMessage.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace RevBayesCore {

    /**
     * The scaling operator.
     *
     * A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5)
     * and scales the current vale by a scaling factor
     * sf = exp( lambda * u )
     * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Michael Landis)
     * @since 2009-09-08, version 1.0
     *
     */

    template<class charType>
    class NodeRejectionSampleProposal : public Proposal {

    public:
        NodeRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l=1.0, double r=0.234);
        // PL comments: should it be *rf be NULL? //!<  constructor
        NodeRejectionSampleProposal( const NodeRejectionSampleProposal& p );                                                        //!<  constructor
        virtual                                                    ~NodeRejectionSampleProposal(void);                              //!< Destructor
        
        NodeRejectionSampleProposal&                                operator=(const NodeRejectionSampleProposal& p);

        // Basic utility functions
        void                                                        assignNode(TopologyNode* nd);
        void                                                        assignSiteIndexSet(const std::set<size_t>& s);
        NodeRejectionSampleProposal*                                clone(void) const;                                              //!< Clone object
        void                                                        cleanProposal(void);
//        virtual double                                              computeLnProposal();
        double                                                      doProposal(void);                                               //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                                      getProposalTuningParameter(void) const;
        double                                                      getRootBranchLength(void);                                     //!< get the length of the root branch
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;                   //!< Print the parameter summary
        void                                                        prepareProposal(void);                                          //!< Prepare the proposal
        std::set<size_t>                                            chooseCharactersToSample(double p);
        void                                                        setSampledCharacters(const std::set<size_t>& s);
        void                                                        sampleNodeCharacters(void);                                     //!< Sample the characters at the node
        void                                                        sampleRootCharacters(void);                                     //!< Sample the characters at the root
        void                                                        setProposalTuningParameter(double tp);
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);         //!< Set the rate generator.
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d); //!< Set the rate generator.
        // void                                                        setRootFrequencies(const TypedDagNode< Simplex > *rf); //!< Set the root frequencies
        void                                                        tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                        undoProposal(void);                                             //!< Reject the proposal

    protected:

        void                                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on

        // parameters
        StochasticNode<AbstractHomologousDiscreteCharacterData>*    ctmc;
        const TypedDagNode<RateGenerator>*                          q_map_site;
        const TypedDagNode<RateGeneratorSequence>*                  q_map_sequence;

        // dimensions
        size_t                                                      numCharacters;
        size_t                                                      numStates;

        // proposal
        std::vector<size_t>                                         storedNodeState;
        std::vector<size_t>                                         storedSubrootState;

        TopologyNode*                                               node;
        double                                                      storedLnProb;
        double                                                      proposedLnProb;

        PathRejectionSampleProposal<charType>*                      nodeProposal;
        PathRejectionSampleProposal<charType>*                      leftProposal;
        PathRejectionSampleProposal<charType>*                      rightProposal;

        TransitionProbabilityMatrix                                 nodeTpMatrix;
        TransitionProbabilityMatrix                                 leftTpMatrix;
        TransitionProbabilityMatrix                                 rightTpMatrix;

        double                                                      lambda;
        std::set<size_t>                                            sampledCharacters;
        std::set<size_t>                                            allCharacters;
        
        // const TypedDagNode< Simplex >*                              rootFrequencies;
    };

}



/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::NodeRejectionSampleProposal<charType>::NodeRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l, double r) : Proposal(r),
    ctmc(n),
    q_map_site( NULL ),
    q_map_sequence( NULL ),
    numCharacters(n->getValue().getNumberOfCharacters()),
    numStates(2),
    node( NULL ),
    nodeTpMatrix(2),
    leftTpMatrix(2),
    rightTpMatrix(2),
    lambda(l)
    // rootFrequencies(rf)
{

    addNode( ctmc );

    nodeProposal  = new PathRejectionSampleProposal<charType>(n, l, r);
    leftProposal  = new PathRejectionSampleProposal<charType>(n, l, r);
    rightProposal = new PathRejectionSampleProposal<charType>(n, l, r);
    
    for (size_t i = 0; i < numCharacters; i++)
    {
        allCharacters.insert(i);
    }

}


/**
 * Copy constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::NodeRejectionSampleProposal<charType>::NodeRejectionSampleProposal( const NodeRejectionSampleProposal& p ) : Proposal(p),
    ctmc( p.ctmc ),
    q_map_site( p.q_map_site ),
    q_map_sequence( p.q_map_sequence ),
    numCharacters( p.numCharacters ),
    numStates( p.numStates ),
    node( p.node ),
    nodeTpMatrix( p.nodeTpMatrix ),
    leftTpMatrix( p.leftTpMatrix ),
    rightTpMatrix( p.rightTpMatrix ),
    lambda( p.lambda )
    // rootFrequencies( p.rootFrequencies )
{

    addNode( ctmc );
    // addNode( q_map_site );
    // addNode( q_map_sequence );

    nodeProposal  = p.nodeProposal->clone();
    leftProposal  = p.leftProposal->clone();
    rightProposal = p.rightProposal->clone();
    
    sampledCharacters   = p.sampledCharacters;
    allCharacters       = p.allCharacters;

    storedNodeState     = p.storedNodeState;
    storedSubrootState  = p.storedSubrootState;;

    storedLnProb        = p.storedLnProb;
    proposedLnProb      = p.proposedLnProb;

}


/**
  * Destructor
  *
  */
 template<class charType>
 RevBayesCore::NodeRejectionSampleProposal<charType>::~NodeRejectionSampleProposal( void )
 {
         
     delete nodeProposal;
     delete leftProposal;
     delete rightProposal;

 }


/**
 * Copy constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::NodeRejectionSampleProposal<charType>& RevBayesCore::NodeRejectionSampleProposal<charType>::operator=( const NodeRejectionSampleProposal& p )
{

    if (this != &p)
    {
        delete nodeProposal;
        delete leftProposal;
        delete rightProposal;
        
        removeNode( ctmc );
        removeNode( q_map_site );
        removeNode( q_map_sequence );
        
        ctmc                = p.ctmc;
        q_map_site          = p.q_map_site;
        q_map_sequence      = p.q_map_sequence;
        numCharacters       = p.numCharacters;
        numStates           = p.numStates;
        node                = p.node;
        nodeTpMatrix        = p.nodeTpMatrix;
        leftTpMatrix        = p.leftTpMatrix;
        rightTpMatrix       = p.rightTpMatrix;
        lambda              = p.lambda;
        // rootFrequencies     = p.rootFrequencies;

        addNode( ctmc );
        addNode( q_map_site );
        addNode( q_map_sequence );
        
        nodeProposal  = p.nodeProposal->clone();
        leftProposal  = p.leftProposal->clone();
        rightProposal = p.rightProposal->clone();
    
        
        sampledCharacters   = p.sampledCharacters;
        allCharacters       = p.allCharacters;

        storedNodeState     = p.storedNodeState;
        storedSubrootState  = p.storedSubrootState;;

        storedLnProb        = p.storedLnProb;
        proposedLnProb      = p.proposedLnProb;
    }

    
    return *this;
}


template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::assignNode(TopologyNode* nd)
{
    node = nd;
}


template<class charType>
std::set<size_t> RevBayesCore::NodeRejectionSampleProposal<charType>::chooseCharactersToSample(double p)
{
    
    std::set<size_t> s;
    s.insert(GLOBAL_RNG->uniform01() * numCharacters);
    for (size_t i = 0; i < numCharacters; i++)
    {
        if (GLOBAL_RNG->uniform01() < p)
        {
            s.insert(i);
        }
    }
    return s;
    
}



template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::cleanProposal( void )
{
    nodeProposal->cleanProposal();
    rightProposal->cleanProposal();
    leftProposal->cleanProposal();

    storedNodeState.clear();
    storedSubrootState.clear();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::NodeRejectionSampleProposal<charType>* RevBayesCore::NodeRejectionSampleProposal<charType>::clone( void ) const
{
    return new NodeRejectionSampleProposal( *this );
}

/**
 * Perform the Proposal.
 *
 * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
 * and scales the current vale by a scaling factor
 * sf = exp( lambda * u )
 * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
 *
 * \return The hastings ratio.
 */
template<class charType>
double RevBayesCore::NodeRejectionSampleProposal<charType>::doProposal( void )
{
    proposedLnProb = 0.0;
    
    double proposedLnProbRatio = 0.0;
    
    // backward proposal probability
    //    proposedLnProbRatio += computeLnProposal();
    
    // TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    
    if ( node->isRoot() )
    {
        std::cout << "Root node." << std::endl;
        sampleRootCharacters();
        
        //  string rtb  = rootBranch;
        double root_branch_length = getRootBranchLength();
        if ( root_branch_length > 0 )
        {
            // update parent incident path
            proposedLnProbRatio += nodeProposal->doProposal();
        }
        // update 2x child incident paths
        proposedLnProbRatio += leftProposal->doProposal();
        proposedLnProbRatio += rightProposal->doProposal();
    }
    else
    {
        std::cout << "Not root node." << std::endl;
        sampleNodeCharacters();
        // update 3x child incident paths
        proposedLnProbRatio += nodeProposal->doProposal();
        proposedLnProbRatio += leftProposal->doProposal();
        proposedLnProbRatio += rightProposal->doProposal();
    }


    // forward proposal probability
//    proposedLnProbRatio -= computeLnProposal();

    return proposedLnProbRatio;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::NodeRejectionSampleProposal<charType>::getProposalName( void ) const
{
    static std::string name = "NodeRejectionSampleProposal";

    return name;
}


template<class charType>
double RevBayesCore::NodeRejectionSampleProposal<charType>::getProposalTuningParameter( void ) const
{
    return lambda;
}


template<class charType>
double RevBayesCore::NodeRejectionSampleProposal<charType>::getRootBranchLength( void )
{
    // PL comments: this method is here because getRootBranchLength() in TreeHistoryCtmc.h gives a fake root branch length if the true one is 0
    // PL comments: can be remove if that is fixed
    
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    const Tree& tree = c->getTree();
    std::vector<TopologyNode*> nds = tree.getNodes();
    node = nds[tree.getRoot().getIndex()];
    
    double origin = c->getRootBranchLength();
    double root_age = node->getAge();
    
    double root_branch_length = ( origin - root_age == 0 ) ? 0 : origin;
}


/**
 *
 */
template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::prepareProposal( void )
{

    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }

    storedLnProb = 0.0;
    proposedLnProb = 0.0;

    // sample a node from the tree that is not a tip node
    const Tree& tree = p->getTree();
    std::vector<TopologyNode*> nds = tree.getNodes();
    node = NULL;
    do
    {
        size_t idx = GLOBAL_RNG->uniform01() * nds.size();
        node = nds[idx];
    } while ( node->isTip() == true );

    // prepare the path proposals

    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    double root_branch_length = getRootBranchLength();
    if ( not ( node->isRoot() && root_branch_length == 0 ) ) {
        nodeProposal->assignNode(node);
        nodeProposal->prepareProposal();
    }

    leftProposal->assignNode(&node->getChild(0));
    leftProposal->prepareProposal();

    rightProposal->assignNode(&node->getChild(1));
    rightProposal->prepareProposal();


    // store node state values
    storedNodeState.clear();
    storedSubrootState.clear();

    size_t num_sites = p->getNumberOfSites();
    storedNodeState.resize(num_sites,0);
    const std::vector<CharacterEvent*>& nodeState = p->getHistory(*node).getChildCharacters();
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(nodeState[site_index])->getState();
        storedNodeState[site_index] = s;
    }
    
    if (node->isRoot()) {
        TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
        storedSubrootState.resize(num_sites,0);
        const std::vector<CharacterEvent*>& subrootState = p->getHistory(*node).getParentCharacters();
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = static_cast<CharacterEventDiscrete*>(subrootState[site_index])->getState();
            storedSubrootState[site_index] = s;
        }
    }

    
    // sample characters to be updated and pass to proposals
    sampledCharacters = chooseCharactersToSample(lambda);
    nodeProposal->setSampledCharacters(sampledCharacters);
    leftProposal->setSampledCharacters(sampledCharacters);
    rightProposal->setSampledCharacters(sampledCharacters);

}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }
}

template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::sampleNodeCharacters( void )
{
    // PL comments: edit this to sample root node
    if ( node->isTip()  )
    {
        return;
    }
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }

    CharacterHistoryDiscrete& histories = c->getHistories();

    // get local node and branch information
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);

    size_t node_index = node->getIndex();
    size_t left_index = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    double node_age   = node->getAge();
    double left_age   = left_child.getAge();
    double right_age  = right_child.getAge();

    double node_rate  = c->getBranchRate( node_index );
    double left_rate  = c->getBranchRate( left_index );
    double right_rate = c->getBranchRate( right_index );

    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& leftChildState  = histories[ left_index ].getChildCharacters();
    const std::vector<CharacterEvent*>& rightChildState = histories[ right_index ].getChildCharacters();

    // states to update
    std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index ].getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState  = histories[ left_index ].getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[ right_index ].getParentCharacters();

    // we may also update this if it is the root state (not const)
    std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()].getParentCharacters();
    
    // get parent age
    double parent_age = 0.0;
    parent_age = node->getParent().getAge();
    
    // get transition probs
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
    rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, nodeTpMatrix);
    rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  leftTpMatrix);
    rm.calculateTransitionProbabilities(node_age, right_age, right_rate, rightTpMatrix);

    
    std::set<size_t>::iterator it_s;

//  for (size_t site_index = 0; site_index < num_sites; ++site_index)
    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
    {
        size_t site_index = *it_s;
        size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
        size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
        size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();

        double u = GLOBAL_RNG->uniform01();

        std::vector<double> state_probs(numStates,0.0);
        double sum = 0.0;
        for ( size_t i=0; i<numStates; ++i )
        {
            double p = nodeTpMatrix[ancS][i] * leftTpMatrix[i][desS1] * rightTpMatrix[i][desS2];
            sum += p;
            state_probs[i] = p;
        }

        unsigned int s = 0;
        for ( size_t i=0; i<numStates; ++i )
        {
            u -= (state_probs[i]/sum);
            if ( u <= 0.0 )
            {
                break;
            }
            ++s;
        }

        static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
    }

}


template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::sampleRootCharacters( void )
{
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // const Simplex rf = c->getRootFrequencies();
    const Simplex& rf = c->getRootFrequencies();
    // TypedDagNode< Simplex >*
    std::cout << "Root frequencies: " << rf << std::endl;
    
    double root_branch_length = getRootBranchLength();
    if ( root_branch_length == 0)
    {
    //     if ( rf.size() == 0 )
    //     {
    //         const RateMatrix *rm = ( q_map_sequence != NULL ? dynamic_cast<const RateMatrix *>( &q_map_sequence->getValue() ) : // dynamic_cast<const RateMatrix *>( &q_map_site->getValue() ) );
    //         if ( rm != NULL )
    //         {
    //             rf = rm->getStationaryFrequencies();
    //             std::cout << "Root frequencies at stationarity calculated: " << rf << std::endl;
    //         }
    //         else
    //         {
    //             throw RbException("You either need to use a rate-matrix or specify root frequencies.");
    //         }
    //     }
        
        CharacterHistoryDiscrete& histories = c->getHistories();
        
        // get local node and branch information
        TopologyNode &left_child  = node->getChild(0);
        TopologyNode &right_child = node->getChild(1);

        size_t node_index = node->getIndex();
        size_t left_index = left_child.getIndex();
        size_t right_index = right_child.getIndex();
        
        double node_age   = node->getAge();
        double left_age   = left_child.getAge();
        double right_age  = right_child.getAge();

        double node_rate  = 0.0;
        double left_rate  = c->getBranchRate( left_index );
        double right_rate = c->getBranchRate( right_index );

        // states for conditional sampling probs
        const std::vector<CharacterEvent*>& leftChildState  = histories[ left_index ].getChildCharacters();
        const std::vector<CharacterEvent*>& rightChildState = histories[ right_index ].getChildCharacters();

        // states to update
        std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index ].getChildCharacters();
        std::vector<CharacterEvent*>& leftParentState  = histories[ left_index ].getParentCharacters();
        std::vector<CharacterEvent*>& rightParentState = histories[ right_index ].getParentCharacters();

        
        // no need when no root branch
        // std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()].getParentCharacters();
        
        // get transition probs
        const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
        rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  leftTpMatrix);
        rm.calculateTransitionProbabilities(node_age, right_age, right_rate, rightTpMatrix);
        
        std::set<size_t>::iterator it_s;

        for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
        {
            size_t site_index = *it_s;
            
            double u = GLOBAL_RNG->uniform01();
            unsigned int s = 0;
            for ( size_t i=0; i<numStates; ++i )
            {
                u -= rf[i];
                if ( u <= 0.0 )
                {
                    break;
                }
                ++s;
            }
            //            std::cout << s;
            static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
            static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
        }
    }
    else
    {
        // if ( rtb == "INPUT" )
        // {
        //     double node_rate  = c->getBranchRate( node_index );
        // }
        // double left_rate  = c->getBranchRate( left_index );
        // double right_rate = c->getBranchRate( right_index );
    //
        // // states for conditional sampling probs
        // const std::vector<CharacterEvent*>& leftChildState  = histories[ left_index ].getChildCharacters();
        // const std::vector<CharacterEvent*>& rightChildState = histories[ right_index ].getChildCharacters();
    //
        // // states to update
        // std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index ].getChildCharacters();
        // std::vector<CharacterEvent*>& leftParentState  = histories[ left_index ].getParentCharacters();
        // std::vector<CharacterEvent*>& rightParentState = histories[ right_index ].getParentCharacters();
    //
        //
        // // we may also update this if it is the root state (not const)
        // std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()].getParentCharacters();
        //
        //
        // // get origin age and root branch length (if present)
        // double parent_age         = 0.0;
        // if ( root_branch_length > 0 )
        // {
        //     parent_age = node->getAge() + c->getRootBranchLength();
        // }
    //
        // // get transition probs
        // const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
        // rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  leftTpMatrix);
        // rm.calculateTransitionProbabilities(node_age, right_age, right_rate, rightTpMatrix);
        // if ( rtb == "INPUT")
        // {
        //     rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, nodeTpMatrix);
        // }
        //
        //
        // std::set<size_t>::iterator it_s;
        // // sample root/origin state
        // if ( rtb == "INPUT")
        // {
        //     for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
        //     {
        //         size_t site_index = *it_s;
        //
        //         double u = GLOBAL_RNG->uniform01();
        //         unsigned int s = 0;
        //         for ( size_t i=0; i<numStates; ++i )
        //         {
        //             u -= rf[i];
        //             if ( u <= 0.0 )
        //             {
        //                 break;
        //             }
        //             ++s;
        //         }
        //         //            std::cout << s;
        //         static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
        //     }
        //
        //     //  for (size_t site_index = 0; site_index < num_sites; ++site_index)
        //     for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
        //     {
        //         size_t site_index = *it_s;
        //
        //         size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
        //         size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
        //         size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
    //
        //
        //         double u = GLOBAL_RNG->uniform01();
    //
        //         std::vector<double> state_probs(numStates,0.0);
        //         double sum = 0.0;
        //         for ( size_t i=0; i<numStates; ++i )
        //         {
        //             double p = nodeTpMatrix[ancS][i] * leftTpMatrix[i][desS1] * rightTpMatrix[i][desS2];
        //             sum += p;
        //             state_probs[i] = p;
        //         }
    //
        //         unsigned int s = 0;
        //         for ( size_t i=0; i<numStates; ++i )
        //         {
        //             u -= (state_probs[i]/sum);
        //             if ( u <= 0.0 )
        //             {
        //                 break;
        //             }
        //             ++s;
        //         }
        //         static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
        //         static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
        //         static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
        //     }
        // }

    }
}


template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{

    // removeNode( q_map_site );
    // removeNode( q_map_sequence );
    
    q_map_site = d;
    numStates = q_map_site->getValue().getNumberOfStates();

    nodeProposal->setRateGenerator( q_map_site );
    leftProposal->setRateGenerator( q_map_site );
    rightProposal->setRateGenerator( q_map_site );

    nodeTpMatrix = TransitionProbabilityMatrix(numStates);
    leftTpMatrix = TransitionProbabilityMatrix(numStates);
    rightTpMatrix = TransitionProbabilityMatrix(numStates);
    
    // addNode( q_map_site );
    // q_map_sequence = NULL;

}


template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{

    // removeNode( q_map_site );
    // removeNode( q_map_sequence );
    
    q_map_sequence = d;
    numStates = q_map_sequence->getValue().getNumberOfStates();

    nodeProposal->setRateGenerator( q_map_sequence );
    leftProposal->setRateGenerator( q_map_sequence );
    rightProposal->setRateGenerator( q_map_sequence );

    nodeTpMatrix = TransitionProbabilityMatrix(numStates);
    leftTpMatrix = TransitionProbabilityMatrix(numStates);
    rightTpMatrix = TransitionProbabilityMatrix(numStates);

    // addNode( q_map_sequence );
    // q_map_site = NULL;
    
}


// template<class charType>
// void RevBayesCore::NodeRejectionSampleProposal<charType>::setRootFrequencies(const TypedDagNode< Simplex > *rf)
// {
//     rootFrequencies = rf;
// }


template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::setSampledCharacters(const std::set<size_t>& s)
{
    sampledCharacters = s;
}

/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    if (oldN == ctmc)
    {
        ctmc = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN) ;
    }
    else if (oldN == q_map_site)
    {
        q_map_site = static_cast<DeterministicNode<RateGenerator>* >(newN);
    }
    else if (oldN == q_map_sequence)
    {
        q_map_sequence = static_cast<DeterministicNode<RateGeneratorSequence>* >(newN);
    }

    nodeProposal->swapNode(oldN, newN);
    leftProposal->swapNode(oldN, newN);
    rightProposal->swapNode(oldN, newN);
}


template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::tune( double rate )
{
    double p = this->targetAcceptanceRate;
    if ( rate > p )
    {
        lambda /= (1.0 + ((rate-p)/(1.0 - p)));
    }
    else
    {
        lambda *= (2.0 - rate/p);
        if (lambda > 1.0)
            lambda = 1.0;
    }
}



/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the ctmc/DAG-node to its original value.
 */
template<class charType>
void RevBayesCore::NodeRejectionSampleProposal<charType>::undoProposal( void )
{
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = p->getNumberOfSites();
    
    CharacterHistoryDiscrete& histories = p->getHistories();
    
    // restore node state
    std::vector<CharacterEvent*>& nodeChildState = histories[node->getIndex()].getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState = histories[node->getChild(0).getIndex() ].getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[node->getChild(1).getIndex()].getParentCharacters();
    
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = storedNodeState[site_index];
        static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
    }
    
    // restore subroot state if needed
    if (node->isRoot()) {
        //        std::cout << "restore subrootState : ";
        std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()].getParentCharacters();
        
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = storedSubrootState[site_index];
            //            std::cout << s;
            static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
        }
        //        std::cout << "\n";
    }
    
    
    // restore path state
    nodeProposal->undoProposal();
    rightProposal->undoProposal();
    leftProposal->undoProposal();
    
    // clear sampled character set
    sampledCharacters.clear();
}

#endif /* defined(__rb_mlandis__NodeRejectionSampleProposal__) */
