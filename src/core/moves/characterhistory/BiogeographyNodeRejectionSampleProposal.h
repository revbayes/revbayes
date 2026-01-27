#ifndef BiogeographicNodeRejectionSampleProposal_H
#define BiogeographicNodeRejectionSampleProposal_H

#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "DeterministicNode.h"
#include "DistributionBinomial.h"
#include "DistributionPoisson.h"
#include "HomologousDiscreteCharacterData.h"
#include "PathRejectionSampleProposal.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateGeneratorSequence.h"
#include "RbException.h"
#include "RbMathCombinatorialFunctions.h"
#include "StochasticNode.h"
//#include "TransitionProbability.h"
#include "TreeChangeEventMessage.h"
#include "GeneralTreeHistoryCtmc.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

#include <algorithm>
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
    class BiogeographicNodeRejectionSampleProposal : public Proposal {
        
    public:
        BiogeographicNodeRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l=1.0, double r=0.234 );                                  //!<  constructor
        
        // Basic utility functions
        void                                                        assignNode(TopologyNode* nd);
        void                                                        assignSiteIndexSet(const std::set<size_t>& s);
        BiogeographicNodeRejectionSampleProposal*                   clone(void) const;                                              //!< Clone object
        void                                                        cleanProposal(void);
                double                                              doProposal(void);                                               //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                                      getProposalTuningParameter(void) const;
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;                   //!< Print the parameter summary
        void                                                        prepareProposal(void);                                          //!< Prepare the proposal
        std::set<size_t>                                            chooseCharactersToSample(double p);
        void                                                        setSampledCharacters(const std::set<size_t>& s);
        double                                                      sampleNodeCharacters(void);                                     //!< Sample the characters at the node
        double                                                      sampleCladogenesisCharacters(void);
        double                                                      sampleRootCharacters(void);                                     //!< Sample the characters at the root
        void                                                        setProposalTuningParameter(double tp);
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);         //!< Set the rate generator.
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d); //!< Set the rate generator.
        void                                                        tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                        undoProposal(void);                                             //!< Reject the proposal
        
    protected:
        virtual double                                              computeNodeLnProposal();
        virtual double                                              computeCladogenesisLnProposal();
        virtual double                                              computeAnagenesisConditionLnProposal();
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
        std::vector<size_t>                                         storedLeftState;
        std::vector<size_t>                                         storedRightState;
        std::vector<size_t>                                         storedSubrootState;
        std::string                                                 storedCladogeneticEvent;
        
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
        
    };
    
}



/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::BiogeographicNodeRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l, double r ) : Proposal(r),
ctmc(n),
q_map_site( NULL ),
q_map_sequence( NULL ),
numCharacters(n->getValue().getNumberOfCharacters()),
numStates(2),
node( NULL ),
nodeTpMatrix(2),
leftTpMatrix(2),
rightTpMatrix(2),
lambda(1.0), // for now, lambda always == 1.0
storedLnProb(0.0),
proposedLnProb(0.0),
storedCladogeneticEvent("")
{
    
    addNode( ctmc );

    nodeProposal = new PathRejectionSampleProposal<charType>(n, l, r);
    leftProposal = new PathRejectionSampleProposal<charType>(n, l, r);
    rightProposal = new PathRejectionSampleProposal<charType>(n, l, r);
    
    for (size_t i = 0; i < numCharacters; i++)
    {
        allCharacters.insert(i);
    }
    
}


template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::cleanProposal( void )
{
    nodeProposal->cleanProposal();
    rightProposal->cleanProposal();
    leftProposal->cleanProposal();
    
    storedNodeState.clear();
    storedLeftState.clear();
    storedRightState.clear();
    storedSubrootState.clear();
    storedCladogeneticEvent = "";
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>* RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::clone( void ) const
{
    return new BiogeographicNodeRejectionSampleProposal( *this );
}


template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::assignNode(TopologyNode* nd)
{
    node = nd;
}


template<class charType>
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::computeAnagenesisConditionLnProposal(void)
{
    double lnP = 0.0;
    if ( node->isTip() )
    {
        return 0.0;
    }
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // collect phylogenetic info
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    
    double parent_age = 0.0;
    if ( node->isRoot() )
    {
        parent_age = node->getAge() + c->getRootBranchLength();
    }
    else
    {
        parent_age = node->getParent().getAge();
    }
    double node_age   = node->getAge();
    double left_age   = left_child.getAge();
    double right_age  = right_child.getAge();
    
    double node_rate  = c->getBranchRate( node->getIndex() );
    double left_rate  = c->getBranchRate( left_child.getIndex() );
    double right_rate = c->getBranchRate( right_child.getIndex() );
    
    size_t node_index = node->getIndex();
    size_t left_index = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    // get transition probs
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
    rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, nodeTpMatrix);
    rm.calculateTransitionProbabilities(node_age, left_age,  left_rate, leftTpMatrix);
    rm.calculateTransitionProbabilities(node_age, right_age, right_rate, rightTpMatrix);
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index  ].getChildCharacters();
    const std::vector<CharacterEvent*>& nodeParentState  = histories[ node_index  ].getParentCharacters();
    const std::vector<CharacterEvent*>& leftChildState   = histories[ left_index  ].getChildCharacters();
    const std::vector<CharacterEvent*>& leftParentState  = histories[ left_index  ].getParentCharacters();
    const std::vector<CharacterEvent*>& rightChildState  = histories[ right_index ].getChildCharacters();
    const std::vector<CharacterEvent*>& rightParentState = histories[ right_index ].getParentCharacters();
    
    // get transition probabilities for all site-states
    std::set<size_t>::iterator it_s;
    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
    {
        size_t site_index = *it_s;
        
        size_t srcSA = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
        size_t desSA = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
        double probA = nodeTpMatrix[srcSA][desSA];
        
        size_t srcS1 = static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->getState();
        size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
        double prob1 = leftTpMatrix[srcS1][desS1];
        
        size_t srcS2 = static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->getState();
        size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
        double prob2 = rightTpMatrix[srcS2][desS2];
        
        
        lnP += log( probA * prob1 * prob2 );
    }
    
    // this is the denominator used to impose a condition on a joint pdf
    // hence, we want the negative of this log probability
    lnP = -lnP;
    
    return lnP;
}


template<class charType>
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::computeCladogenesisLnProposal(void)
{
    
    double p = 1.0;
    
    if (node->isTip()) {
        return 0.0;
    }
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    /*
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& leftParentState  = histories[node->getChild(0).getIndex()]->getParentCharacters();
    const std::vector<CharacterEvent*>& rightParentState = histories[node->getChild(1).getIndex()]->getParentCharacters();
    const std::vector<CharacterEvent*>& nodeChildState  = histories[node->getIndex()]->getChildCharacters();
    
    // determine the cladogenetic state based on un/shared areas
    size_t n_n_on = 0;
    size_t n_l_on = 0;
    size_t n_r_on = 0;
    size_t n_nlr_off = 0;
    size_t n_lr_mismatch = 0;
    size_t n_nlr_on = 0;
    size_t n_jump_mismatch = 0;
    for (size_t i = 0; i < nodeChildState.size(); i++)
    {
        size_t sn = static_cast<CharacterEventDiscrete*>(nodeChildState[i])->getState();
        size_t sl = static_cast<CharacterEventDiscrete*>(leftParentState[i])->getState();
        size_t sr = static_cast<CharacterEventDiscrete*>(rightParentState[i])->getState();
        
        n_n_on += sn;
        n_l_on += sl;
        n_r_on += sr;
        
        // copy sympatry
        if (sn==1&&sl==1&&sr==1)
        {
            n_nlr_on += 1;
        }
    
        // jump dispersal
        if ( (sn==0&&sr==1&&sl==0) || (sn==0&&sr==0&&sl==1) )
        {
            n_jump_mismatch += 1;
        }
        
        // allopatry
        if ( (sn==1&&sl==1&&sr==0) || (sn==1&&sl==0&&sr==1) )
        {
            n_lr_mismatch += 1;
        }
        
        if (sn==0&&sr==0&&sl==0) {
            n_nlr_off += 1;
        }
        
    }
    
    if (n_nlr_on > 1) {
        
        throw RbException("unknown cladogenetic state");
    }
    
    
    std::string clado_type = "";
    if (n_n_on == 0) {
        clado_type = "null_copy";
    }
    else if (n_nlr_on == 1 && n_nlr_off == (numCharacters-1))
    {
        clado_type = "sympatry_copy";
    }
    else if (n_n_on == n_lr_mismatch)
    {
        clado_type = "allopatry";
    }
    else if (n_jump_mismatch == 1)
    {
        clado_type = "jump_dispersal";
    }
    else if (n_nlr_on == 1 && n_lr_mismatch==(n_n_on-n_nlr_on))
    {
        clado_type = "sympatry_subset";
    }
    else
    {
        
        throw RbException("Unknown cladogenetic event type!");
    }
    
    std::string clado_type_2 = c->getCladogeneticEvent( node->getIndex() );
    if (clado_type != clado_type_2) {
        
        throw RbException("weird!");
    }
    */
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& nodeChildState  = histories[node->getIndex()].getChildCharacters();
    
    // determine the cladogenetic state based on un/shared areas
    size_t n_n_on = 0;
    for (size_t i = 0; i < nodeChildState.size(); i++)
    {
        size_t sn = static_cast<CharacterEventDiscrete*>(nodeChildState[i])->getState();
        n_n_on += sn;
    }
    std::string clado_type = c->getCladogeneticEvent( node->getIndex() );
    
    // the proposal prob
    std::map<std::string, double> clado_probs;
    clado_probs[ "sympatry_subset" ] = 0.5;
    clado_probs[ "allopatry" ]       = 0.5;
    clado_probs[ "jump_dispersal" ]  = 0.0;

    if ( clado_type == "null_copy") {
        p = 1.0;
    }
    else if ( clado_type == "sympatry_copy" )
    {
        p = 1.0;
    }
    else if ( clado_type == "allopatry" )
    {
        // Any combination of bits for one range (with the sister having its complement)
        // excluding the all-zero range and the all-one range (hence, -2)
        size_t n_events = std::pow(2, (double)n_n_on) - 2;
        p = clado_probs[ "allopatry" ] * (1.0 / n_events);
    }
    else if ( clado_type == "sympatry_subset" )
    {
        // Any single ancestral bit may be set across the two daughter ranges
        size_t n_events = 2 * n_n_on;
        p = clado_probs[ "sympatry_subset" ] * (1.0 / n_events);
    }
    else if ( clado_type == "jump_dispersal" )
    {
        // Any single non-ancestral bit may be set across the two daughter ranges
        size_t n_events = 2 * (numCharacters - n_n_on);
        p = clado_probs[ "jump_dispersal" ] * (1.0 / n_events);
    }
    else
    {
        throw RbException() << "cladogenetic event type \"" << clado_type << "\" not found!"; 
    }
    
    double lnP = log(p);
    
    return lnP;
}

template<class charType>
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::computeNodeLnProposal(void)
{
    
    // only apply move if this is an internal node
    if ( node->isTip() )
    {
        return 0.0;
    }
    
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    
    // save forward sampling prob
    double lnP = 0.0;
    
    // gather branch history variables
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    double parent_age = 0.0;
    if ( node->isRoot() )
    {
        parent_age = node->getAge() + c->getRootBranchLength();
    }
    else
    {
        parent_age = node->getParent().getAge();
    }
    double node_age    = node->getAge();
    double left_age    = left_child.getAge();
    double right_age   = right_child.getAge();
    
    double node_rate   = c->getBranchRate( node_index );
    double left_rate   = c->getBranchRate( left_index );
    double right_rate  = c->getBranchRate( right_index );
    
    // get transition probs
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
    rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, nodeTpMatrix);
    rm.calculateTransitionProbabilities(node_age, left_age,  left_rate, leftTpMatrix);
    rm.calculateTransitionProbabilities(node_age, right_age, right_rate, rightTpMatrix);
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& nodeParentState = histories[ node_index ].getParentCharacters();
    const std::vector<CharacterEvent*>& leftChildState  = histories[ left_index ].getChildCharacters();
    const std::vector<CharacterEvent*>& rightChildState = histories[ right_index ].getChildCharacters();
    
    // states to update
    const std::vector<CharacterEvent*>& nodeChildState  = histories[ node_index ].getChildCharacters();
    
    // compute sampling prob
    std::set<size_t>::iterator it_s;
    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
    {
        size_t site_index = *it_s;
        size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
        size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
        size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
        
        std::vector<double> state_probs(numStates,0.0);
        double sum = 0.0;
        for ( size_t i=0; i<numStates; ++i )
        {
            double p = nodeTpMatrix[ancS][i] * leftTpMatrix[i][desS1] * rightTpMatrix[i][desS2];
            sum += p;
            state_probs[i] = p;
        }
        
        // get the probability of the sampled state *before* cladogenesis
        size_t s = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
        lnP += log( state_probs[s]/sum );
    }

    return lnP;
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
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::doProposal( void )
{
    proposedLnProb = 0.0;
    
    //    // compute the reverse proposal probability;
    //    // note, we later extract the reverse proposal probabilities
    //    // from the three path samplers when calling doProposal()
//    std::cout << "REVERSE\n";
    storedLnProb = 0.0;
    double p_node_bwd = computeNodeLnProposal();
    double p_clado_bwd = computeCladogenesisLnProposal();
    double p_ana_cond_bwd = computeAnagenesisConditionLnProposal();
    storedLnProb = p_node_bwd + p_clado_bwd + p_ana_cond_bwd;
    
    // this is the hastings ratio
    double proposedLnProbRatio = 0.0;
    
    // update ancestral node state
    double p_node_fwd = sampleNodeCharacters();
    // update daughter states following cladogenesis
    double p_clado_fwd = sampleCladogenesisCharacters();
    // account for the conditional probability for anagenetic rej. sampling
    // NOTE: this term cancels with the node proposal density when
    //       daughter lineages identically inherit the ancestral range
    //       but not when cladogenetic events are present!
    double p_ana_cond_fwd = computeAnagenesisConditionLnProposal();
    
    proposedLnProb = p_node_fwd + p_clado_fwd + p_ana_cond_fwd;

    // compute the proposal density ratio for the node sampler
    proposedLnProbRatio = storedLnProb - proposedLnProb;
    
    // update 3x incident paths, and account for the joint probabilities
    // of sampling those path histories and ending in whatever states;
    // their contribution to the ratio also accounts for the reverse
    // proposal probabilities
    proposedLnProbRatio += nodeProposal->doProposal();
    proposedLnProbRatio += leftProposal->doProposal();
    proposedLnProbRatio += rightProposal->doProposal();
    
    // return hastings ratio
    return proposedLnProbRatio;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::getProposalName( void ) const
{
    static std::string name = "BiogeographicNodeRejectionSampleProposal";
    
    return name;
}


template<class charType>
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 *
 */
template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::prepareProposal( void )
{
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // sample an internal node
    const Tree& tree = c->getTree();
    std::vector<TopologyNode*> nds = tree.getNodes();
    node = NULL;
    do
    {
        size_t idx = GLOBAL_RNG->uniform01() * nds.size();
        node = nds[idx];
    } while ( node->isTip() );
    
    
    // sample characters to be updated and pass to proposals
    sampledCharacters = chooseCharactersToSample(lambda);
    nodeProposal->setSampledCharacters(sampledCharacters);
    leftProposal->setSampledCharacters(sampledCharacters);
    rightProposal->setSampledCharacters(sampledCharacters);
    
    // re-initialize proposal probabilities
    storedLnProb = 0.0;
    proposedLnProb = 0.0;
    
    // prepare the path proposals;
    // this computes the reverse proposal probs, among other things
    nodeProposal->assignNode(node);
    nodeProposal->prepareProposal();
    
    leftProposal->assignNode(&node->getChild(0));
    leftProposal->prepareProposal();
    
    rightProposal->assignNode(&node->getChild(1));
    rightProposal->prepareProposal();
    
    // store node,left,right cladogenetic state values
    storedNodeState.clear();
    storedLeftState.clear();
    storedRightState.clear();
    storedSubrootState.clear();
    
    const TopologyNode& left_node = node->getChild(0);
    const TopologyNode& right_node = node->getChild(1);
    const std::vector<CharacterEvent*>& nodeChildState   = c->getHistory(*node).getChildCharacters();
    const std::vector<CharacterEvent*>& leftParentState  = c->getHistory(left_node).getParentCharacters();
    const std::vector<CharacterEvent*>& rightParentState = c->getHistory(right_node).getParentCharacters();
    size_t num_sites = c->getNumberOfSites();
    storedNodeState.resize(num_sites,0);
    storedLeftState.resize(num_sites,0);
    storedRightState.resize(num_sites,0);
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t desSA = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
        size_t srcS1 = static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->getState();
        size_t srcS2 = static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->getState();
        storedNodeState[site_index] = desSA;
        storedLeftState[site_index] = srcS1;
        storedRightState[site_index] = srcS2;
    }
    
    // store cladogenetic event state
    storedCladogeneticEvent = c->getCladogeneticEvent( node->getIndex() );
    
    // store the subroot state if the root node is updated
    if (node->isRoot()) {
        storedSubrootState.resize(num_sites,0);
        const std::vector<CharacterEvent*>& subrootState = c->getHistory(*node).getParentCharacters();
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = static_cast<CharacterEventDiscrete*>(subrootState[site_index])->getState();
            storedSubrootState[site_index] = s;
        }
    }

    return;
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
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    //    o << "lambda = " << lambda;
}

template<class charType>
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::sampleNodeCharacters( void )
{
    
    // only apply move if this is an internal node
    if ( node->isTip() )
    {
        return 0.0;
    }
    
    // verify tree history object
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }

    // save forward sampling prob
    double lnP = 0.0;

    // gather branch history variables
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    double parent_age = 0.0;
    if ( node->isRoot() ) {
        parent_age = node->getAge() + c->getRootBranchLength();
    } else {
        parent_age = node->getParent().getAge();
    }
    double node_age    = node->getAge();
    double left_age    = left_child.getAge();
    double right_age   = right_child.getAge();
    
    double node_rate   = c->getBranchRate( node_index );
    double left_rate   = c->getBranchRate( left_index );
    double right_rate  = c->getBranchRate( right_index );
    
    // get transition probs
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
    
    rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, nodeTpMatrix);
    rm.calculateTransitionProbabilities(node_age, left_age,  left_rate,  leftTpMatrix);
    rm.calculateTransitionProbabilities(node_age, right_age, right_rate, rightTpMatrix);
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& leftChildState  = histories[ left_index ].getChildCharacters();
    const std::vector<CharacterEvent*>& rightChildState = histories[ right_index ].getChildCharacters();
    
    // states to update
    std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index ].getChildCharacters();
    std::vector<CharacterEvent*>& nodeParentState  = histories[ node_index ].getParentCharacters();
    std::vector<CharacterEvent*>& leftParentState  = histories[ left_index ].getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[ right_index ].getParentCharacters();
 
    // sample states
    std::set<size_t>::iterator it_s;
    if ( node->isRoot() )
    {
        for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
        {
            size_t site_index = *it_s;
            
            double u = GLOBAL_RNG->uniform01();
            std::vector<double> state_probs(numStates, 1.0/numStates);
            unsigned int s = 0;
            for ( size_t i=0; i<numStates; ++i )
            {
                u -= state_probs[i];
                if ( u <= 0.0 )
                {
                    break;
                }
                ++s;
            }
//            size_t old_s = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
            static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
        }
    }

    // sample node characters
    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
    {
        size_t site_index = *it_s;
        size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
        size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
        size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();
        
        std::vector<double> state_probs(numStates,0.0);
        double sum = 0.0;
        for ( size_t i=0; i<numStates; ++i )
        {
            double p = nodeTpMatrix[ancS][i] * leftTpMatrix[i][desS1] * rightTpMatrix[i][desS2];
            sum += p;
            state_probs[i] = p;
        }
        
        double u = GLOBAL_RNG->uniform01() * sum;
        unsigned int s = 0;
        for ( size_t i=0; i<numStates; ++i )
        {
            u -= state_probs[i];
            if ( u <= 0.0 )
            {
                //                        s = i;
                break;
            }
            ++s;
        }
        
        lnP += std::log( state_probs[s] / sum );
        
        static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(s);
        static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(s);
    }
    
    return lnP;
 
    
}



/**
 * Sample the cladogenetic states following speciation.
 *
 * Subdivides the ancestral range into right and left daughter ranges based on the cladogenetic model.
 */

template<class charType>
double RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::sampleCladogenesisCharacters(void)
{
    
    // verify tree history object
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // forward proposal probability
    double lnP = 0.0;
    
    // gather branch history variables
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);

    // states for conditional sampling probs
//    const std::vector<CharacterEvent*>& leftChildState  = histories[left_child.getIndex()]->getChildCharacters();
//    const std::vector<CharacterEvent*>& rightChildState = histories[right_child.getIndex()]->getChildCharacters();
    
    // states to update
    std::vector<CharacterEvent*>& nodeChildState   = histories[node->getIndex()].getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState  = histories[left_child.getIndex()].getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[right_child.getIndex()].getParentCharacters();
    
    // initialize by setting both daughter lineages to be identical to the parent state
    std::vector<size_t> idx_on;
    std::vector<size_t> idx_off;
    for (size_t i = 0; i < nodeChildState.size(); i++)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(nodeChildState[i])->getState();
        static_cast<CharacterEventDiscrete*>(leftParentState[i])->setState(s);
        static_cast<CharacterEventDiscrete*>(rightParentState[i])->setState(s);
        if (s == 1)
        {
            idx_on.push_back(i);
        }
        else
        {
            idx_off.push_back(i);
        }
        
    }
    
    // if the ancestral node has 0 or 1 areas, then that range is inherited
    // identically by both daughter lineages
    size_t n_anc_1 = idx_on.size();
    
    // sample a cladogenetic event type
    // NOTE: need to add a way to toggle what events are used;
    //       this will involve querying the clado prob mtx
    std::map<std::string, double> clado_probs;
    clado_probs[ "sympatry_subset" ] = 0.5;
    clado_probs[ "allopatry" ] = 0.5;
    clado_probs[ "jump_dispersal" ] = 0.0;

    std::string clado_type = "";
    double u_type = GLOBAL_RNG->uniform01();
    for (std::map<std::string, double>::iterator it = clado_probs.begin(); it != clado_probs.end(); it++)
    {
        u_type -= it->second;
        if (u_type <= 0.0) {
            clado_type = it->first;
            break;
        }
    }

    if (n_anc_1 == 0) {
        clado_type = "null_copy";
    } else if (n_anc_1 == 1) {
        clado_type = "sympatry_copy";
    }
    
    // sample cladogenetic states
    if (clado_type == "null_copy")
    {
        lnP = 0.0;
    }
    else if (clado_type == "sympatry_copy")
    {
        lnP = 0.0;
    }
    else if (clado_type == "sympatry_subset")
    {
        // which branch is the trunk and, implicitly, which is the branch?)
        bool left_branch_is_trunk = GLOBAL_RNG->uniform01() < 0.5;
        
        // which characters are inherited by the branch?
        std::vector<size_t> bud_char_idx;
        bud_char_idx.push_back( (size_t)(GLOBAL_RNG->uniform01() * idx_on.size() ) );
        
        // sample the cladogenetic state depending on bud/trunk and character indices
        if (left_branch_is_trunk) {
            
            // set all right-parent nodes to zero
            for (size_t i = 0; i < idx_on.size(); i++) {
                size_t j = idx_on[i];
                static_cast<CharacterEventDiscrete*>(rightParentState[j])->setState(0);
            }
            // set the sampled right-parent nodes to one
            for (std::vector<size_t>::iterator it = bud_char_idx.begin(); it != bud_char_idx.end(); it++)
            {
                size_t j = idx_on[*it];
                static_cast<CharacterEventDiscrete*>(rightParentState[j])->setState(1);
            }
        }
        else {
            // set all left-parent nodes to zero
            for (size_t i = 0; i < idx_on.size(); i++) {
                size_t j = idx_on[i];
                static_cast<CharacterEventDiscrete*>(leftParentState[j])->setState(0);
            }
            // set the sampled left-parent nodes to one
            for (std::vector<size_t>::iterator it = bud_char_idx.begin(); it != bud_char_idx.end(); it++)
            {
                size_t j = idx_on[*it];
                static_cast<CharacterEventDiscrete*>(leftParentState[j])->setState(1);
            }
        }
        
        size_t n_events = 2 * idx_on.size();
        lnP = std::log( clado_probs["sympatry_subset"] * (1.0 / n_events) );
    }
    else if (clado_type == "allopatry")
    {
        // sample how left/right characters will be divided
        deprecated::random_shuffle(idx_on.begin(), idx_on.end());
        int n_child1_on = (int)(GLOBAL_RNG->uniform01() * (idx_on.size()-1) ) + 1;
        
        // set the left/right child states
        for (size_t i = 0; i < idx_on.size(); i++) {
            
            size_t j = idx_on[i];

            if (n_child1_on > 0) {
                static_cast<CharacterEventDiscrete*>(leftParentState[j])->setState(0);
                static_cast<CharacterEventDiscrete*>(rightParentState[j])->setState(1);
            } else {
                static_cast<CharacterEventDiscrete*>(leftParentState[j])->setState(1);
                static_cast<CharacterEventDiscrete*>(rightParentState[j])->setState(0);
            }
            
            n_child1_on -= 1;
        }
        size_t n_events = std::pow((double)(2.0), (double)(idx_on.size())) - 2;
        lnP = std::log( clado_probs["allopatry"] * (1.0 / n_events) );
    }
    else if (clado_type == "jump_dispersal")
    {
        // which branch is the trunk and, implicitly, which is the branch?)
        bool left_branch_is_trunk = GLOBAL_RNG->uniform01() < 0.5;
        
        // select an area that is *not* part of the ancestral range
        size_t bud_char_idx = idx_off[ (size_t)GLOBAL_RNG->uniform01() * idx_off.size() ];
        
        // sample the cladogenetic state depending on bud/trunk and character indices
        if (left_branch_is_trunk) {
            for (size_t i = 0; i < idx_on.size(); i++) {
                size_t j = idx_on[i];
                static_cast<CharacterEventDiscrete*>(rightParentState[j])->setState(0);
            }
            static_cast<CharacterEventDiscrete*>(rightParentState[bud_char_idx])->setState(1);
            
        }
        else {
            for (size_t i = 0; i < idx_on.size(); i++) {
                size_t j = idx_on[i];
                static_cast<CharacterEventDiscrete*>(leftParentState[j])->setState(0);
            }
            static_cast<CharacterEventDiscrete*>(leftParentState[bud_char_idx])->setState(1);
            
        }
        size_t n_events =  2 * (numCharacters - idx_on.size());
        lnP = std::log( clado_probs["jump_dispersal"] * (1.0 / n_events) );
        
    }
    else
    {
        
        throw RbException(clado_type + " not a recognized cladogenetic event type!");
    }
    
    c->setCladogeneticEvent( clado_type, node->getIndex() );
    
    return lnP;
}


template<class charType>
std::set<size_t> RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::chooseCharactersToSample(double p)
{
    
    if (p >= 1.0)
    {
        return allCharacters;
    }
    
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
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::setSampledCharacters(const std::set<size_t>& s)
{
    sampledCharacters = s;
}


template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{
    
    q_map_site = d;
    numStates = q_map_site->getValue().getNumberOfStates();
    
    nodeProposal->setRateGenerator( q_map_site );
    leftProposal->setRateGenerator( q_map_site );
    rightProposal->setRateGenerator( q_map_site );
    
    nodeTpMatrix = TransitionProbabilityMatrix(numStates);
    leftTpMatrix = TransitionProbabilityMatrix(numStates);
    rightTpMatrix = TransitionProbabilityMatrix(numStates);
    
}


template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{
    
    q_map_sequence = d;
    numStates = q_map_sequence->getValue().getNumberOfStates();
    
    nodeProposal->setRateGenerator( q_map_sequence );
    leftProposal->setRateGenerator( q_map_sequence );
    rightProposal->setRateGenerator( q_map_sequence );
    
    nodeTpMatrix = TransitionProbabilityMatrix(numStates);
    leftTpMatrix = TransitionProbabilityMatrix(numStates);
    rightTpMatrix = TransitionProbabilityMatrix(numStates);
    
}


/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
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
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType>
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::tune( double rate )
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
void RevBayesCore::BiogeographicNodeRejectionSampleProposal<charType>::undoProposal( void )
{
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = c->getNumberOfSites();
    
    CharacterHistoryDiscrete& histories = c->getHistories();
    
    // restore node state
    std::vector<CharacterEvent*>& nodeChildState   = histories[node->getIndex()].getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState  = histories[node->getChild(0).getIndex() ].getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[node->getChild(1).getIndex()].getParentCharacters();
    
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t desSA = storedNodeState[site_index];
        size_t srcS1 = storedLeftState[site_index];
        size_t srcS2 = storedRightState[site_index];
        
        static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(desSA);
        static_cast<CharacterEventDiscrete*>(leftParentState[site_index])->setState(srcS1);
        static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->setState(srcS2);
    }
    
    // restore subroot state if needed
    if (node->isRoot()) {
//        std::cout << "restore subrootState : ";
        std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()].getParentCharacters();
        
        for (size_t site_index = 0; site_index < num_sites; ++site_index)
        {
            size_t s = storedSubrootState[site_index];
//                        std::cout << s;
            static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
        }
//        std::cout << "\n";
    }
    
    // restore cladogenetic event type
    c->setCladogeneticEvent( storedCladogeneticEvent, node->getIndex() );
    
    // restore path state
    nodeProposal->undoProposal();
    rightProposal->undoProposal();
    leftProposal->undoProposal();
    
    // clear sampled character set
    sampledCharacters.clear();
}

#endif /* defined(__rb_mlandis__BiogeographicNodeRejectionSampleProposal__) */
