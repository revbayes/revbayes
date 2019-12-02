//#ifndef NodeUniformizationSampleProposal_H
//#define NodeUniformizationSampleProposal_H
//
//#include "BranchHistory.h"
//#include "DeterministicNode.h"
//#include "HomologousDiscreteCharacterData.h"
//#include "DistributionBinomial.h"
//#include "DistributionPoisson.h"
//#include "PathUniformizationSampleProposal.h"
//#include "Proposal.h"
//#include "RandomNumberFactory.h"
//#include "RandomNumberGenerator.h"
//#include "RateGeneratorSequence.h"
//#include "RbException.h"
//#include "StochasticNode.h"
////#include "TransitionProbability.h"
//#include "TopologyNode.h"
//#include "TypedDagNode.h"
//
//#include <cmath>
//#include <iostream>
//#include <set>
//#include <string>
//
//namespace RevBayesCore {
//    
//    /**
//     * The scaling operator.
//     *
//     * A scaling proposal draws a random uniform number u ~ unif (-0.5,0.5)
//     * and scales the current vale by a scaling factor
//     * sf = exp( lambda * u )
//     * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
//     *
//     * @copyright Copyright 2009-
//     * @author The RevBayes Development Core Team (Michael Landis)
//     * @since 2009-09-08, version 1.0
//     *
//     */
//    
//    template<class charType>
//    class NodeUniformizationSampleProposal : public Proposal {
//        
//    public:
//        NodeUniformizationSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, TypedDagNode<Tree>* t, double l, TopologyNode* nd=NULL );                                                                //!<  constructor
//        NodeUniformizationSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, TypedDagNode<Tree>* t, PathUniformizationSampleProposal<charType>* p, double l, TopologyNode* nd=NULL );                                                                //!<  constructor
//        
//        // Basic utility functions
//        void                                    assignNode(TopologyNode* nd);
//        void                                    assignSiteIndexSet(const std::set<size_t>& s);
//        NodeUniformizationSampleProposal*       clone(void) const;                                                                  //!< Clone object
//        void                                    cleanProposal(void);
//        double                                  doProposal(void);                                                                   //!< Perform proposal
//        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
//        void                                    printParameterSummary(std::ostream &o) const;                                       //!< Print the parameter summary
//        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
//        void                                    sampleNodeCharacters(const std::set<size_t>& indexSet);   //!< Sample the characters at the node
//        //        void                            sampleNodeCharacters2(const TopologyNode& node, const std::set<size_t>& indexSet);   //!< Sample the characters at the node
//        double                                  sampleRootCharacters(const std::set<size_t>& indexSet);
//        void                                    setRateGenerator(const TypedDagNode<RateGenerator> *d);                                   //!< Set the rate generator.
//        void                                    setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d);                                   //!< Set the rate generator.
//        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
//        void                                    undoProposal(void);                                                                 //!< Reject the proposal
//        
//    protected:
//        
//        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
//        
//        
//        // parameters
//        StochasticNode<AbstractHomologousDiscreteCharacterData>*  ctmc;
//        TypedDagNode<Tree>*                     tau;
//        const TypedDagNode<RateGenerator>*            q_map_site;
//        const TypedDagNode<RateGeneratorSequence>*    q_map_sequence;
//        
//        // dimensions
//        size_t                                  numNodes;
//        size_t                                  numCharacters;
//        size_t                                  numStates;
//        
//        // proposal
//        std::vector<size_t>                     storedNodeState;
//        std::vector<size_t>                     storedRootState;
//        TopologyNode*                           node;
//        std::set<size_t>                        siteIndexSet;
//        double                                  storedLnProb;
//        double                                  proposedLnProb;
//        
//        PathUniformizationSampleProposal<charType>* nodeProposal;
//        PathUniformizationSampleProposal<charType>* leftProposal;
//        PathUniformizationSampleProposal<charType>* rightProposal;
//        
//        TransitionProbabilityMatrix nodeTpMatrix;
//        TransitionProbabilityMatrix leftTpMatrix;
//        TransitionProbabilityMatrix rightTpMatrix;
//        
//        double                                  lambda;
//        
//        // flags
//        bool                                    fixNodeIndex;
//        bool                                    sampleNodeIndex;
//        bool                                    sampleSiteIndexSet;
//        
//    };
//    
//}
//
//
//
///**
// * Constructor
// *
// * Here we simply allocate and initialize the Proposal object.
// */
//template<class charType>
//RevBayesCore::NodeUniformizationSampleProposal<charType>::NodeUniformizationSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, TypedDagNode<Tree> *t, double l, TopologyNode* nd) : Proposal(),
//    ctmc(n),
//    tau(t),
//    q_map_site(NULL),
//    q_map_sequence(NULL),
//    numNodes(t->getValue().getNumberOfNodes()),
//    numCharacters(n->getValue().getNumberOfCharacters()),
//    numStates(2),
//    node(nd),
//    nodeTpMatrix(numStates),
//    leftTpMatrix(numStates),
//    rightTpMatrix(numStates),
//    lambda(l),
//    sampleNodeIndex(true),
//    sampleSiteIndexSet(true)
//{
//    
//    //    std::cout << numStates << "\n";
//    addNode(ctmc);
//    addNode(tau);
//    
//    nodeProposal = new PathUniformizationSampleProposal<charType>(n,t,l,nd);
//    leftProposal = new PathUniformizationSampleProposal<charType>(n,t,l,nd);
//    rightProposal = new PathUniformizationSampleProposal<charType>(n,t,l,nd);
//    
//    fixNodeIndex = (node != NULL);
//}
//
//
//template<class charType>
//RevBayesCore::NodeUniformizationSampleProposal<charType>::NodeUniformizationSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, TypedDagNode<Tree> *t, PathUniformizationSampleProposal<charType>* p, double l, TopologyNode* nd) : Proposal(),
//    ctmc(n),
//    tau(t),
//    q_map_site(NULL),
//    q_map_sequence(NULL),
//    numNodes(t->getValue().getNumberOfNodes()),
//    numCharacters(n->getValue().getNumberOfCharacters()),
//    numStates(2),
//    node(nd),
//    nodeTpMatrix(numStates),
//    leftTpMatrix(numStates),
//    rightTpMatrix(numStates),
//    lambda(l),
//    sampleNodeIndex(true),
//    sampleSiteIndexSet(true)
//{
//    
//    //    std::cout << numStates << "\n";
//    addNode( ctmc );
//    addNode( tau );
//    
//    nodeProposal = new PathUniformizationSampleProposal<charType>(*p);
//    leftProposal = new PathUniformizationSampleProposal<charType>(*p);
//    rightProposal = new PathUniformizationSampleProposal<charType>(*p);
//    
//    fixNodeIndex = (node != NULL);
//}
//
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::cleanProposal( void )
//{
//    //    AbstractTreeHistoryCtmc<charType>* p = static_cast< AbstractTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
//    //    const std::vector<BranchHistory*>& histories = p->getHistories();
//    
//    nodeProposal->cleanProposal();
//    rightProposal->cleanProposal();
//    leftProposal->cleanProposal();
//}
//
///**
// * The clone function is a convenience function to create proper copies of inherited objected.
// * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
// *
// * \return A new copy of the proposal.
// */
//template<class charType>
//RevBayesCore::NodeUniformizationSampleProposal<charType>* RevBayesCore::NodeUniformizationSampleProposal<charType>::clone( void ) const
//{
//    return new NodeUniformizationSampleProposal( *this );
//}
//
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::assignNode(TopologyNode* nd)
//{
//    node = nd;
//    sampleNodeIndex = false;
//}
//
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::assignSiteIndexSet(const std::set<size_t>& s)
//{
//    siteIndexSet = s;
//    sampleSiteIndexSet = false;
//}
//
//
///**
// * Get Proposals' name of object
// *
// * \return The Proposals' name.
// */
//template<class charType>
//const std::string& RevBayesCore::NodeUniformizationSampleProposal<charType>::getProposalName( void ) const
//{
//    static std::string name = "NodeUniformizationSampleProposal";
//    
//    return name;
//}
//
//
///**
// * Perform the Proposal.
// *
// * A scaling Proposal draws a random uniform number u ~ unif (-0.5,0.5)
// * and scales the current vale by a scaling factor
// * sf = exp( lambda * u )
// * where lambda is the tuning parameter of the Proposal to influence the size of the proposals.
// *
// * \return The hastings ratio.
// */
//template<class charType>
//double RevBayesCore::NodeUniformizationSampleProposal<charType>::doProposal( void )
//{
//    proposedLnProb = 0.0;
//    
//    double proposedLnProbRatio = 0.0;
//    
//    //    AbstractTreeHistoryCtmc<charType>* p = static_cast< AbstractTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
//    //    const std::vector<BranchHistory*>& histories = p->getHistories();
//    
//    // update node value
//    
//    //    if (node->isRoot())
//    //        sampleRootCharacters(siteIndexSet);
//    
//    sampleNodeCharacters(siteIndexSet);
//    
//    // update 3x incident paths
//    proposedLnProbRatio += nodeProposal->doProposal();
//    proposedLnProbRatio += leftProposal->doProposal();
//    proposedLnProbRatio += rightProposal->doProposal();
//    
//    return proposedLnProbRatio;
//}
//
//
///**
// *
// */
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::prepareProposal( void )
//{
//    AbstractTreeHistoryCtmc<charType>* p = static_cast< AbstractTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
//    
//    storedLnProb = 0.0;
//    proposedLnProb = 0.0;
//    
//    const Tree& tree = tau->getValue();
//    //    size_t numTips = tree.getNumberOfTips();
//    if (sampleNodeIndex && !fixNodeIndex)
//    {
//        node = NULL;
//        std::vector<TopologyNode*> nds = tree.getNodes();
//        while (node == NULL || node->isTip()) {
//            size_t idx = GLOBAL_RNG->uniform01() * nds.size(); //numTips + GLOBAL_RNG->uniform01() * (numNodes-numTips);
//            node = nds[idx];
//        };
//    }
//    sampleNodeIndex = true;
//    
//    if (sampleSiteIndexSet)
//    {
//        siteIndexSet.clear();
//        siteIndexSet.insert(GLOBAL_RNG->uniform01() * numCharacters); // at least one is inserted
//        for (size_t i = 0; i < numCharacters; i++)
//        {
//            if (GLOBAL_RNG->uniform01() < lambda)
//            {
//                siteIndexSet.insert(i);
//            }
//        }
//    }
//    sampleSiteIndexSet = true;
//    
//    // prepare the path proposals
//    nodeProposal->assignNode(node);
//    nodeProposal->assignSiteIndexSet(siteIndexSet);
//    nodeProposal->prepareProposal();
//    
//    leftProposal->assignNode(&node->getChild(0));
//    leftProposal->assignSiteIndexSet(siteIndexSet);
//    leftProposal->prepareProposal();
//    
//    rightProposal->assignNode(&node->getChild(1));
//    rightProposal->assignSiteIndexSet(siteIndexSet);
//    rightProposal->prepareProposal();
//    
//    // store node state values
//    storedNodeState.clear();
//    storedRootState.clear();
//    
//    storedNodeState.resize(numCharacters,0);
//    const std::vector<CharacterEvent*>& nodeState = p->getHistory(*node).getChildCharacters();
//    for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//    {
//        size_t s = nodeState[*it]->getState();
//        storedNodeState[*it] = s;
//    }
//    if (node->isRoot())
//    {
//        storedRootState.resize(numCharacters,0);
//        const std::vector<CharacterEvent*>& rootState = p->getHistory(*node).getParentCharacters();
//        for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//        {
//            size_t s = rootState[*it]->getState();
//            storedRootState[*it] = s;
//        }
//    }
//}
//
//
///**
// * Print the summary of the Proposal.
// *
// * The summary just contains the current value of the tuning parameter.
// * It is printed to the stream that it passed in.
// *
// * \param[in]     o     The stream to which we print the summary.
// */
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::printParameterSummary(std::ostream &o) const
//{
//    o << "lambda = " << lambda;
//}
//
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::sampleNodeCharacters(const std::set<size_t>& indexSet)
//{
//    
//    if (!node->isTip())
//    {
//        bool isRoot = node->isRoot();
//        AbstractTreeHistoryCtmc<charType>* p = static_cast< AbstractTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
//        std::vector<BranchHistory*> histories = p->getHistories();
//        
//        // get transition probs
//        const RateGenerator& rm = q_map_site->getValue();
//        
//        TopologyNode &left_child  = node->getChild(0);
//        TopologyNode &right_child = node->getChild(1);
//        
//        double node_age     = node->getAge();
//        double left_age     = left_child.getAge();
//        double right_age    = right_child.getAge();
//        
//        double node_rate  = 1.0;
//        double left_rate  = 1.0;
//        double right_rate = 1.0;
//
//        rm.calculateTransitionProbabilities(leftTpMatrix, node_age, left_age, left_rate);
//        rm.calculateTransitionProbabilities(rightTpMatrix, node_age, right_age, right_rate);
//        if ( isRoot == false )
//        {
//            double parent_age   = node->getParent().getAge();
//            rm.calculateTransitionProbabilities(nodeTpMatrix, parent_age, node_age, node_rate);
//        }
//        
//        // states for conditional sampling probs
//        const std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()]->getParentCharacters();
//        const std::vector<CharacterEvent*>& leftChildState  = histories[node->getChild(0).getIndex()]->getChildCharacters();
//        const std::vector<CharacterEvent*>& rightChildState = histories[node->getChild(1).getIndex()]->getChildCharacters();
//
//        std::vector<CharacterEvent*> nodeChildState   = histories[node->getIndex()]->getChildCharacters();
//        std::vector<CharacterEvent*> leftParentState  = histories[node->getChild(0).getIndex()]->getParentCharacters();
//        std::vector<CharacterEvent*> rightParentState = histories[node->getChild(1).getIndex()]->getParentCharacters();
//        
//        for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//        {
//                  
//            size_t desS1 = leftChildState[*it]->getState();
//            size_t desS2 = rightChildState[*it]->getState();
//            
//            std::vector<double> g(numStates, 0.0);
//            double gSum = 0.0;
//            for (size_t i = 0; i < numStates; i++)
//            {
//                g[i] = leftTpMatrix[i][desS1] * rightTpMatrix[i][desS2];
//                if (!isRoot)
//                {
//                    size_t ancS = nodeParentState[*it]->getState();
//                    g[i] *= nodeTpMatrix[ancS][i];
//                }
//                gSum += g[i];
//            }
//            
//            double u = GLOBAL_RNG->uniform01() * gSum;
//            unsigned int s = 0;
//            for (unsigned i = 0; i < numStates; i++)
//            {
//                u -= g[i];
//                if (u <= 0.0)
//                {
//                    s = i;
//                    break;
//                }
//            }
//        
//            leftParentState[*it]->setState(s);
//            rightParentState[*it]->setState(s);
//            if (!isRoot)
//                nodeChildState[*it]->setState(s);
//
//        }
//    }
//}
//
//template<class charType>
//double RevBayesCore::NodeUniformizationSampleProposal<charType>::sampleRootCharacters(const std::set<size_t>& indexSet)
//{
//    double lnP = 0.0;
//    if (!node->isRoot())
//        return 0.0;
//    
//    AbstractTreeHistoryCtmc<charType>* p = static_cast< AbstractTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
//    BranchHistory* bh = &p->getHistory(*node);
//    std::vector<CharacterEvent*> parentState = bh->getParentCharacters();
//    
//    double r0 = q_map_site->getValue().getRate(1,0);
//    double r1 = q_map_site->getValue().getRate(0,1);
//    unsigned n1_old = 0;
//    unsigned n1_new = 0;
//    
//    double p1 = r1 / (r0 + r1);
//    for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//    {
//        unsigned s = 0;
//        double u = GLOBAL_RNG->uniform01();
//        if (u < p1)
//        {
//            s = 1;
//            n1_new++;
//        }
//        if (parentState[*it]->getState() == 1)
//        {
//            n1_old++;
//        }
//        
//        parentState[*it]->setState(s); //new CharacterEvent(*it,s,0.0);
//    }
//    //bh->setParentCharacters(parentState);
//    
//    size_t n = siteIndexSet.size();
//    size_t n0_old = n - n1_old;
//    size_t n0_new = n - n1_new;
//    double p0 = 1.0 - p1;
//    
//    lnP = n1_old * log(p1) + n0_old * log(p0) - n1_new * log(p1) - n0_new * log(p0);
//    return 0.0;
//    return lnP;
//}
//
//
///**
// * Reject the Proposal.
// *
// * Since the Proposal stores the previous value and it is the only place
// * where complex undo operations are known/implement, we need to revert
// * the value of the ctmc/DAG-node to its original value.
// */
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::undoProposal( void )
//{
//    AbstractTreeHistoryCtmc<charType>* p = static_cast< AbstractTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
//    const std::vector<BranchHistory*>& histories = p->getHistories();
//    
//    // restore path state
//    nodeProposal->undoProposal();
//    rightProposal->undoProposal();
//    leftProposal->undoProposal();
//    
//    // restore node state
//    std::vector<CharacterEvent*> nodeChildState = histories[node->getIndex()]->getChildCharacters();
//    std::vector<CharacterEvent*> leftParentState = histories[node->getChild(0).getIndex() ]->getParentCharacters();
//    std::vector<CharacterEvent*> rightParentState = histories[node->getChild(1).getIndex()]->getParentCharacters();
//    
//    for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//    {
//        size_t s = storedNodeState[*it];
//        nodeChildState[*it]->setState(s);
//        leftParentState[*it]->setState(s);
//        rightParentState[*it]->setState(s);
//    }
//    
//    // restore root state
//    if (node->isRoot())
//    {
//        std::vector<CharacterEvent*> rootState = histories[node->getIndex()]->getParentCharacters();
//        for (std::set<size_t>::iterator it = siteIndexSet.begin(); it != siteIndexSet.end(); it++)
//        {
//            size_t s = storedRootState[*it];
//            rootState[*it]->setState(s);
//        }
//        
//    }
//}
//
//
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
//{
//    
//    q_map_site = d;
//    numStates = q_map_site->getValue().getNumberOfStates();
//    
//}
//
//
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
//{
//    
//    q_map_sequence = d;
//    numStates = q_map_sequence->getValue().getNumberOfStates();
//    
//}
//
//
///**
// * Swap the current ctmc for a new one.
// *
// * \param[in]     oldN     The old ctmc that needs to be replaced.
// * \param[in]     newN     The new ctmc.
// */
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
//{
//    
//    if (oldN == ctmc)
//    {
//        ctmc = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN) ;
//    }
//    else if (oldN == tau)
//    {
//        tau = static_cast<StochasticNode<Tree>* >(newN);
//    }
//    else if (oldN == q_map_sequence)
//    {
//        q_map_sequence = static_cast<TypedDagNode<RateGeneratorSequence>* >(newN);
//    }
//    else if (oldN == q_map_site)
//    {
//        q_map_site = static_cast<TypedDagNode<RateGenerator>* >(newN);
//    }
//    
//    nodeProposal->swapNode(oldN, newN);
//    leftProposal->swapNode(oldN, newN);
//    rightProposal->swapNode(oldN, newN);
//}
//
//
///**
// * Tune the Proposal to accept the desired acceptance ratio.
// */
//template<class charType>
//void RevBayesCore::NodeUniformizationSampleProposal<charType>::tune( double rate )
//{
//    ; // do nothing
//}
//
//#endif /* defined(__rb_mlandis__NodeUniformizationSampleProposal__) */
