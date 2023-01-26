#ifndef BiogeographicCladogeneticRejectionShiftProposal_H
#define BiogeographicCladogeneticRejectionShiftProposal_H

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
    class BiogeographicCladogeneticRejectionShiftProposal : public Proposal {
        
    public:
        BiogeographicCladogeneticRejectionShiftProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l=1.0, double r=0.234 );                                  //!<  constructor
        
        // Basic utility functions
        void                                                        assignNode(TopologyNode* nd);
        void                                                        assignSiteIndexSet(const std::set<size_t>& s);
        BiogeographicCladogeneticRejectionShiftProposal*            clone(void) const;                                              //!< Clone object
        void                                                        cleanProposal(void);
        double                                                      doProposal(void);                                               //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;                   //!< Print the parameter summary
        void                                                        prepareProposal(void);                                          //!< Prepare the proposal
        
        double getProposalTuningParameter(void) const;
        void setProposalTuningParameter(double tp);
        
        
        void                                                        setSampledCharacters(const std::set<size_t>& s);
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);         //!< Set the rate generator.
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d); //!< Set the rate generator.
        void                                                        tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                        undoProposal(void);                                             //!< Reject the proposal
        
    protected:
        void                                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        double                                                      assignCharactersToSample(double p);
        double                                                      computeAllopatryLnProposal(void);
        double                                                      computeAnagenesisConditionLnProposal();
        double                                                      computeJumpDispersalLnProposal(void);
        double                                                      computeNodeLnProposal(std::string clado_type);
        double                                                      computeSwapLnProposal(void);
        double                                                      computeSympatryLnProposal(void);
        void                                                        computeTransitionProbabilities(void);
        void                                                        initializeCladogeneticEventTypes(void);
        void                                                        printCladogeneticState(void);
        std::string                                                 sampleCladogenesisEventType(void);
        double                                                      sampleRootCharacters(void);                                     //!< Sample the characters at the root
        double                                                      sampleNodeCharacters(void);                                     //!< Sample the characters at the node
        double                                                      sampleAllopatryCharacters(void);
        double                                                      sampleSympatryCharacters(void);
        double                                                      sampleJumpDispersalCharacters(void);
        void                                                        setSwapCharacters(void);
        
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
        std::string                                                 proposedCladogeneticEvent;
        
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
//        std::set<size_t>                                            expandCharacters;
//        std::set<size_t>                                            contractCharacters;
        std::set<size_t>                                            sampledCharacters;
        std::vector<size_t>                                         onCharacters;
        std::vector<size_t>                                         offCharacters;
        std::set<size_t>                                            allCharacters;
        std::vector<size_t>                                         swapCharacters;
        std::vector<size_t>                                         expandCharacters;
        std::vector<size_t>                                         contractCharacters;
        std::vector<size_t>                                         leftCharacters;
        std::vector<size_t>                                         rightCharacters;

        std::map<std::string, double>                               cladogeneticSamplingProbabilities;
        
        std::string                                                 proposal_type;
        double                                                      expand_contract_lnProb;
        bool                                                        perform_swap;
        
        size_t                                                      n_left_fwd, n_left_bwd, n_right_fwd, n_right_bwd;
        
    };
    
}



/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::BiogeographicCladogeneticRejectionShiftProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l, double r ) : Proposal(r),
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
    
    initializeCladogeneticEventTypes();
    
}


template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::cleanProposal( void )
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
RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>* RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::clone( void ) const
{
    return new BiogeographicCladogeneticRejectionShiftProposal( *this );
}


template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::assignNode(TopologyNode* nd)
{
    node = nd;
}


template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::assignCharactersToSample(double p)
{
    
    // reset characters
    onCharacters.clear();
    offCharacters.clear();
    sampledCharacters.clear();
    swapCharacters.clear();
    expandCharacters.clear();
    contractCharacters.clear();
    proposal_type = "";
    
    // cladogenetic event type typically remains constant
    // unless edge case encountered
    proposedCladogeneticEvent = storedCladogeneticEvent;
    
    // get tree history
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    // collect on and off characters
    size_t node_index  = node->getIndex();
    const std::vector<CharacterEvent*>& nodeChildState = histories[ node_index  ]->getChildCharacters();
    for (size_t i = 0; i < nodeChildState.size(); i++) {
        size_t j = static_cast<CharacterEventDiscrete*>( nodeChildState[i] )->getState();
        if (j == 1) {
            onCharacters.push_back( i );
        }
        else {
            offCharacters.push_back( i );
        }
    }
    
    // number of on/off characters
    size_t n_on = onCharacters.size();
    size_t n_off = offCharacters.size();
    
    // get forwards/backwards proposal probabilities
    double p_fwd = 0.5;
    double p_bwd = 0.5;
    if (n_on == 1) {
        proposal_type = "expand";
        p_fwd = 1.0;
    }
    else if (n_off == 0) {
        proposal_type = "contract";
        p_fwd = 1.0;
    }
    else if (GLOBAL_RNG->uniform01() < 0.5) {
        proposal_type = "expand";
    }
    else {
        proposal_type = "contract";
    }
    
    // reverse move must be contract
    if ( proposal_type == "expand" && n_on == (numCharacters-1) ) {
        p_bwd = 1.0;
    }
    
    // reverse move must be expand
    if ( proposal_type == "contract" && n_on == 2 ) {
        p_bwd = 1.0;
    }
    
    // add one left- and one right-character to "swap"
    perform_swap = false;
    if (GLOBAL_RNG->uniform01() < 0.5) {
        perform_swap = true;
        setSwapCharacters();
        for (size_t i = 0; i < swapCharacters.size(); i++) {
            sampledCharacters.insert( swapCharacters[i] );
        }
    }
    
    // population expand/contract/sampled characters
    if (proposal_type == "expand")
    {
        
        // singleton range requires that we resample region
        if ( n_on == 1) {
            expandCharacters.push_back( onCharacters[0] );
            perform_swap = false;
        }

        // sample an off character to sample
        size_t u = (size_t)( GLOBAL_RNG->uniform01() * offCharacters.size() );
        expandCharacters.push_back( offCharacters[u] );

        // make sure expanded character paths are updated
        for (size_t i = 0; i < expandCharacters.size(); i++) {
            sampledCharacters.insert( expandCharacters[i] );
        }
    }
    else if (proposal_type == "contract")
    {
        
        // do not swap
        if ( n_on == 2 ) {
            perform_swap = false;
            contractCharacters.push_back( onCharacters[0] );
            contractCharacters.push_back( onCharacters[1] );
        }
        else {
            size_t u = -1;
            size_t idx = -1;
            bool is_duplicate = true;
            
            do {
                u = (size_t)( GLOBAL_RNG->uniform01() * onCharacters.size() );
                idx = onCharacters[u];
                is_duplicate = ( sampledCharacters.find(idx) != sampledCharacters.end() );
            } while ( is_duplicate );
            
            contractCharacters.push_back( idx );
        }
        
        // make sure contracted character paths are updated
        for (size_t i = 0; i < contractCharacters.size(); i++) {
            sampledCharacters.insert( contractCharacters[i] );
        }
        
    } else {
        throw RbException("Unknown proposal_type");
    }
    
    
    // more proposal bookkeeping for edge cases involving range size of 1
    bool expansion_from_1 = (n_on == 1 && proposal_type == "expand");
    bool contraction_from_2 = (n_on == 2 && proposal_type == "contract");
    if ( expansion_from_1 || contraction_from_2 ) {
        
        // if range expands from size 1
        if ( expansion_from_1 ) {
            
            std::string clado_type = "";
            double p_clado = 1.0;
            double u_type = GLOBAL_RNG->uniform01();
            std::map<std::string, double>::iterator it;
            for (it = cladogeneticSamplingProbabilities.begin(); it != cladogeneticSamplingProbabilities.end(); it++)
            {
                u_type -= it->second;
                if (u_type <= 0.0) {
                    p_clado = it->second;
                    clado_type = it->first;
                    break;
                }
            }
            
            // override
//            clado_type = "s";
//            p_clado = 1.0;
            
            // assign new clado state
            c->setCladogeneticEvent( clado_type, node->getIndex() );
            proposedCladogeneticEvent = clado_type;
            
            // collect forward prob
            p_fwd *= p_clado;
        }
        // if range contracts from size 2
        else if (contraction_from_2 ) {
            
            // get bwd expansion prob
            double p_clado = cladogeneticSamplingProbabilities[ storedCladogeneticEvent ];
//            p_clado = 1.0;
            
            // only one possible cladogenetic event type for all-1 statetwe
            proposedCladogeneticEvent = "s";
            
            // collect backward prob
            p_bwd *= p_clado;
        }
    }

    
    
    // account for probability of selecting random gain/lost area
    double p_fwd_area = 1.0;
    double p_bwd_area = 1.0;
    if (proposal_type == "expand") {
        p_fwd_area = (1.0 / (n_off));
        p_bwd_area = (1.0 / (n_on+1));
    } else if (proposal_type == "contract") {
        p_fwd_area = (1.0 / (n_on));
        p_bwd_area = (1.0 / (n_off+1));
    }
    // possibly other ratios to consider for L/R swap
    

//    
//    std::cout << proposal_type << "\n";
//    std::cout << storedCladogeneticEvent << " -> " << proposedCladogeneticEvent << "\n";
//    std::cout << "n_on = " << n_on << "  " << "n_off = " << n_off << "\n";
//    std::cout << "p_fwd = " << p_fwd << "\n";
//    std::cout << "p_bwd = " << p_bwd << "\n";
//    std::cout << "p_bwd / p_fwd = " << p_bwd / p_fwd << "\n";
//    std::cout << "p_fwd_area = " << p_fwd_area << "\n";
//    std::cout << "p_bwd_area = " << p_bwd_area << "\n";
//    std::cout << "p_bwd_area / p_fwd_area = " << p_bwd_area / p_fwd_area << "\n";
//    std::cout << "p_MH_ratio = " << (p_bwd*p_bwd_area) / (p_fwd*p_fwd_area) << "\n";
//
//    
    if ( expansion_from_1 ) {
        
        ;
    }
    if ( contraction_from_2 ) {
        
        ;
    }
    
    // return
    return log(  (p_bwd*p_bwd_area) / (p_fwd*p_fwd_area) );
}




template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeAnagenesisConditionLnProposal(void)
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
    std::vector<BranchHistory*> histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index = node->getIndex();
    size_t left_index = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    // states for conditional sampling probs
    const std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index  ]->getChildCharacters();
    const std::vector<CharacterEvent*>& nodeParentState  = histories[ node_index  ]->getParentCharacters();
    const std::vector<CharacterEvent*>& leftChildState   = histories[ left_index  ]->getChildCharacters();
    const std::vector<CharacterEvent*>& leftParentState  = histories[ left_index  ]->getParentCharacters();
    const std::vector<CharacterEvent*>& rightChildState  = histories[ right_index ]->getChildCharacters();
    const std::vector<CharacterEvent*>& rightParentState = histories[ right_index ]->getParentCharacters();
    
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
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeAllopatryLnProposal(void)
{
    
    return 0.0;
    double lnP = 0.0;
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    const std::vector<CharacterEvent*>& leftChildState      = histories[ left_index  ]->getChildCharacters();
    const std::vector<CharacterEvent*>& rightChildState     = histories[ right_index ]->getChildCharacters();
    const std::vector<CharacterEvent*>& leftParentState     = histories[ left_index  ]->getParentCharacters();
    const std::vector<CharacterEvent*>& rightParentState    = histories[ right_index ]->getParentCharacters();
    const std::vector<CharacterEvent*>& nodeChildState      = histories[ node_index  ]->getChildCharacters();
    const std::vector<CharacterEvent*>& nodeParentState     = histories[ node_index  ]->getParentCharacters();
    
    //    printCladogeneticState();
    //    std::cout << "\n";
    
    // Sample cladogenetic states for L/R pairs

    
    // if forward proposal_type is "contract", we want the reverse "expand" proposal probability
    if ( proposal_type == "contract" ) {
        
        // iterates over sites
        std::set<size_t>::iterator it_s;
        
        for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
        {
            size_t site_index = *it_s;
            
            // start states
            size_t A_i = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index] )->getState();
            size_t R_i = static_cast<CharacterEventDiscrete*>(rightParentState[site_index])->getState();
            size_t L_i = static_cast<CharacterEventDiscrete*>(leftParentState[site_index] )->getState();
            
            // end states
            size_t A_j = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index]  )->getState();
            size_t R_j = static_cast<CharacterEventDiscrete*>(rightChildState[site_index] )->getState();
            size_t L_j = static_cast<CharacterEventDiscrete*>(leftChildState[site_index]  )->getState();
            
            // get transition probs
            std::vector<double> state_probs(3, 0.0);
//            state_probs[0] = nodeTpMatrix[A_i][0] * leftTpMatrix[0][L_j] * rightTpMatrix[0][R_j];
            state_probs[0] = nodeTpMatrix[A_i][1] * leftTpMatrix[1][L_j] * rightTpMatrix[0][R_j];
            state_probs[1] = nodeTpMatrix[A_i][1] * leftTpMatrix[0][L_j] * rightTpMatrix[1][R_j];
            double sum = state_probs[0] + state_probs[1]; // + state_probs[2];
            
            // determine the area's cladogenetic state
            size_t s = 0;
//            if (A_j == 0 && L_i == 0 && R_i == 0) {
//                s = 0;
//            } else
            if (A_j == 1 && L_i == 1 && R_i == 0) {
                s = 0;
            } else if (A_j == 1 && L_i == 0 && R_i == 1) {
                s = 1;
            } else {
                std::cout << "node_index = " << node->getIndex() << "\n";
                std::cout << "node = " << node_index << "\n";
                throw RbException("unknown allopatric state!");
            }
            
            // record the state's probability
            lnP += log( state_probs[s]/sum );
        }
    }
    // if forward proposal_type is "expand", we want the reverse "contract" proposal probability
    else if (proposal_type == "expand") {
        
        // this probability is handled elsewhere
        lnP = 0.0;
    }
    
    
    
    return lnP;
    
}


template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeSwapLnProposal(void)
{
    double lnP = 0.0;
    double p_fwd = 0.0;
    double p_bwd = 0.0;
    
//    size_t n_left_fwd, n_right_fwd, n_left_bwd, n_right_bwd;
    
    p_fwd = (1.0 / n_left_fwd) * (1.0 / n_right_fwd);
    p_bwd = (1.0 / n_left_bwd) * (1.0 / n_right_bwd);
    
    return lnP;
}

template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeSympatryLnProposal(void)
{
    return 0.0;
    
    if (onCharacters.size() == 1) {
        return 0.0;
    }
    
    // MJL: Note to self. The proposal either adds to the left or right trunk w.p. 1.
    
    
    double lnP = 0.0;
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    // get some history states (will get more later)
    const std::vector<CharacterEvent*>& rightParentState = histories[ right_index ]->getParentCharacters();
    const std::vector<CharacterEvent*>& nodeParentState  = histories[ node_index  ]->getParentCharacters();
    //    const std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index  ]->getChildCharacters();
    
    // which branch is the trunk?
    std::string trunk_branch = "left";
    size_t n_right_on = 0;
    for (size_t i = 0; i < numCharacters; i++)
    {
        n_right_on += static_cast<CharacterEventDiscrete*>(rightParentState[i])->getState();
        if (n_right_on > 1) {
            trunk_branch = "right";
            break;
        }
    }
    
    // Retrieve history states depending on bud/trunk sampling
    size_t trunk_index = left_index;
    size_t bud_index = right_index;
    if (trunk_branch == "left") {
        ; // do nothing
    }
    else if (trunk_branch == "right") {
        trunk_index = right_index;
        bud_index = left_index;
    }
    else {
        throw RbException("trunk_branch state unknown");
    }
    
    //    std::cout << "compute cladogenetic ln prob\n";
    //    printCladogeneticState();
 
    // the forward move is proposal_type == "contract", so the reverse move is "expand"
    if (proposal_type == "contract") {
        // condition on for sampling
        const std::vector<CharacterEvent*>& trunkChildState   = histories[ trunk_index ]->getChildCharacters();
        const std::vector<CharacterEvent*>& budChildState     = histories[ bud_index   ]->getChildCharacters();
        const std::vector<CharacterEvent*>& trunkParentState  = histories[ trunk_index ]->getParentCharacters();
        const std::vector<CharacterEvent*>& budParentState    = histories[ bud_index   ]->getParentCharacters();
        
        const TransitionProbabilityMatrix& trunkTpMatrix = ( trunk_branch == "left" ? leftTpMatrix  : rightTpMatrix );
        const TransitionProbabilityMatrix& budTpMatrix   = ( trunk_branch == "left" ? rightTpMatrix : leftTpMatrix  );
        
        // Sample cladogenetic states for trunk
        std::set<size_t>::iterator it_s;
        std::vector<size_t> trunk_areas;
        for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
        {
            // get ancestral and trunk states
            size_t site_index = *it_s;
            size_t A_i = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
            size_t T_i = static_cast<CharacterEventDiscrete*>(trunkParentState[site_index])->getState();
//            size_t B_i = static_cast<CharacterEventDiscrete*>(budParentState[site_index])->getState();
            size_t T_j = static_cast<CharacterEventDiscrete*>(trunkChildState[site_index])->getState();
            size_t B_j = static_cast<CharacterEventDiscrete*>(budChildState[site_index])->getState();
            
            // get transition probs
            std::vector<double> state_probs(2, 0.0);
            state_probs[0] = nodeTpMatrix[A_i][0] * trunkTpMatrix[0][T_j] * budTpMatrix[0][B_j];
            state_probs[1] = nodeTpMatrix[A_i][1] * trunkTpMatrix[1][T_j] * budTpMatrix[1][B_j];
            double sum = state_probs[0] + state_probs[1];
            
            // get the actual sampled state
            size_t s = T_i;
            if (s == 1) {
                trunk_areas.push_back( site_index );
            }
            lnP += log( state_probs[s] / sum );
        }
        
        // Compute bud state sampling probs
        size_t n_trunk_areas = trunk_areas.size();
        std::vector<double> bud_probs_on( n_trunk_areas, 0.0 );
        double bud_sum_on = 0.0;
        size_t bud_area_on = 0;
        for (size_t i = 0; i < trunk_areas.size(); i++)
        {
            size_t site_index = trunk_areas[i];
            
            size_t B_i = static_cast<CharacterEventDiscrete*>(budParentState[site_index])->getState();
            size_t B_j = static_cast<CharacterEventDiscrete*>(budChildState[site_index])->getState();
            
            std::vector<double> state_probs(2, 0.0);
            state_probs[0] = budTpMatrix[0][B_j]; // * rightTpMatrix[R_i][R_j];
            state_probs[1] = budTpMatrix[1][B_j]; // * rightTpMatrix[R_i][R_j];
            
            bud_probs_on[i] = state_probs[1] / (state_probs[0] + state_probs[1]);
            bud_sum_on += bud_probs_on[i];
            
            if (B_i == 1) {
                bud_area_on = site_index;
            }
        }
        
        lnP += std::log( bud_probs_on[ bud_area_on ] / bud_sum_on );
    }
    // if forward proposal_type is "expand", we want the reverse "contract" proposal probability
    else if (proposal_type == "expand") {
        
        // this probability is handled elsewhere
        lnP = 0.0;
    }

    
    return lnP;
    
}


template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeJumpDispersalLnProposal(void)
{
    double lnP = 0.0;
    
    return lnP;
}


template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeNodeLnProposal(std::string clado_type)
{
    
    double lnP = 0.0;
    
    if (clado_type == "s") {
        lnP = computeSympatryLnProposal();
    } else if (clado_type == "a") {
        lnP = computeAllopatryLnProposal();
    } else if (clado_type == "j") {
        lnP = computeJumpDispersalLnProposal();
    } else {
        
        throw RbException("Unexpected proposedCladogeneticEvent type!");
    }
    
    return lnP;
}


template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::computeTransitionProbabilities(void)
{
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
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
    
    return;
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
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::doProposal( void )
{
//    std::cout << "BEFORE\n";
//    printCladogeneticState();
    
    double proposedLnProbRatio = 0.0;
    
    // compute expand-contract move probability
    // (carried over from prepareProposal()
    double p_expand_contract = expand_contract_lnProb; //sampleExpandContractProposal();
    
    // compute the backwards conditional transition probabilities
    double p_ana_cond_bwd = computeAnagenesisConditionLnProposal();
    
    // compute the root proposal probability ratio
    double p_root_ratio = sampleRootCharacters();
    
    // compute the cladogenetic proposal probability ratio
    double p_clado_ratio = sampleNodeCharacters();
    
    // compute the forwards conditional transition probabilities
    double p_ana_cond_fwd = computeAnagenesisConditionLnProposal();
    
    // count probability ratio for expand/contract move
    proposedLnProbRatio = p_expand_contract;
    
    // count probability ratio for sampled node states
    proposedLnProbRatio += p_clado_ratio;
    
    // update 3x incident paths, and account for the joint probabilities
    // of sampling those path histories and ending in whatever states;
    // their contribution to the ratio also accounts for the reverse
    // proposal probabilities
    proposedLnProbRatio += nodeProposal->doProposal();
    proposedLnProbRatio += leftProposal->doProposal();
    proposedLnProbRatio += rightProposal->doProposal();
    
    // include the conditional probabilities for rejection sampling
    double p_ana_cond_ratio = (p_ana_cond_bwd - p_ana_cond_fwd);
    
    proposedLnProbRatio += p_ana_cond_ratio;
    proposedLnProbRatio += p_root_ratio;
    
//    std::cout << "AFTER\n";
//    printCladogeneticState();
//    std::cout << "\n";
    
    // return hastings ratio
    return proposedLnProbRatio;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::getProposalName( void ) const
{
    static std::string name = "BiogeographicCladogeneticRejectionShiftProposal";
    
    return name;
}




template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::getProposalTuningParameter(void) const
{
    return lambda;
}



template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::initializeCladogeneticEventTypes(void)
{
    
    // verify tree history object
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    //    DECCladogeneticStateFunction
//    const AbstractCladogenicStateFunction* acf = dynamic_cast<const AbstractCladogenicStateFunction* >( &homogeneousCladogeneticProbabilityMatrix->getFunction() );
//    
//    const DECCladogeneticStateFunction* dcf = dynamic_cast<const DECCladogeneticStateFunction* >( acf );
    
//    std::vector<std::string> event_types = dcf->getEventTypes();
    
    const CladogeneticProbabilityMatrix& p = c->getCladogeneticProbabilityMatrix();
    std::vector<std::string> event_types = p.getEventTypes();
    
    for (size_t i = 0; i < event_types.size(); i++) {
        cladogeneticSamplingProbabilities[ event_types[i] ] = 1.0 / event_types.size();
    }
    
    return;
    //
    //    size_t node_index = nd->getIndex();
    //    size_t left_index = nd->getChild(0).getIndex();
    //    size_t right_index = nd->getChild(1).getIndex();
    //    std::string clado_state = cf->simulateDataAugmentedCladogeneticState( this->histories, node_index, left_index, right_index );
    
    
    // we need to learn this dynamically from the ctmc's clado matrix
//    cladogeneticSamplingProbabilities[ "s" ] = 0.5;
//    cladogeneticSamplingProbabilities[ "a" ] = 0.5;
//    cladogeneticSamplingProbabilities[ "j" ] = 0.0;
}



/**
 *
 */
template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::prepareProposal( void )
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
    
    // store cladogenetic event state
    storedCladogeneticEvent = c->getCladogeneticEvent( node->getIndex() );

    // sample characters to be updated and pass to proposals
    expand_contract_lnProb = assignCharactersToSample(lambda);
    
    // set characters for be updates for incident branches
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
    
    // compute the transition probability matrices
    computeTransitionProbabilities();
    
    return;
}


template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::printCladogeneticState(void)
{
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // save forward sampling prob
    //    double lnP = 0.0;
    //    double ln_prob_fwd = 0.0;
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    const std::vector<CharacterEvent*>& leftChildState      = histories[ left_index  ]->getChildCharacters();
    const std::vector<CharacterEvent*>& rightChildState     = histories[ right_index ]->getChildCharacters();
    const std::vector<CharacterEvent*>& leftParentState     = histories[ left_index  ]->getParentCharacters();
    const std::vector<CharacterEvent*>& rightParentState    = histories[ right_index ]->getParentCharacters();
    const std::vector<CharacterEvent*>& nodeChildState      = histories[ node_index ]->getChildCharacters();
    const std::vector<CharacterEvent*>& nodeParentState     = histories[ node_index ]->getParentCharacters();
    
    std::cout << "L0 : ";
    for (size_t i = 0; i < numCharacters; i++) { std::cout << static_cast<CharacterEventDiscrete*>( leftParentState[i] )->getState(); }
    std::cout << "\n";
    std::cout << "R0 : ";
    for (size_t i = 0; i < numCharacters; i++) { std::cout << static_cast<CharacterEventDiscrete*>( rightParentState[i] )->getState(); }
    std::cout << "\n";
    std::cout << "A1 : ";
    for (size_t i = 0; i < numCharacters; i++) { std::cout << static_cast<CharacterEventDiscrete*>( nodeChildState[i] )->getState(); }
    std::cout << "\n";
    
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
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    //    o << "lambda = " << lambda;
}


template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleAllopatryCharacters()
{
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // get phylo coordinates
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
//    const std::vector<CharacterEvent*>& leftChildState      = histories[ left_index  ]->getChildCharacters();
//    const std::vector<CharacterEvent*>& rightChildState     = histories[ right_index ]->getChildCharacters();
//    const std::vector<CharacterEvent*>& nodeParentState     = histories[ node_index ]->getParentCharacters();
    
    std::vector<CharacterEvent*>& leftParentState           = histories[ left_index  ]->getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState          = histories[ right_index ]->getParentCharacters();
    std::vector<CharacterEvent*>& nodeChildState            = histories[ node_index ]->getChildCharacters();
    
    // Sample cladogenetic states for trunk
    std::vector<size_t> trunk_areas;
    
//    std::cout << "BEFORE " << proposal_type << "\n";
//    printCladogeneticState();
    
    // forward proposal is "expand"
    if ( proposal_type == "expand" ) {
        
//        std::set<size_t>::iterator it_s;
        for (size_t i = 0; i < expandCharacters.size(); i++)
        {
            // get ancestral and trunk states
            size_t site_index = expandCharacters[i];
            
            if ( GLOBAL_RNG->uniform01() < 0.5 ) {
                static_cast<CharacterEventDiscrete*>( nodeChildState[site_index]   )->setState( 1 );
                static_cast<CharacterEventDiscrete*>( leftParentState[site_index]  )->setState( 1 );
                static_cast<CharacterEventDiscrete*>( rightParentState[site_index] )->setState( 0 );
            }
            else {
                static_cast<CharacterEventDiscrete*>( nodeChildState[site_index]   )->setState( 1 );
                static_cast<CharacterEventDiscrete*>( leftParentState[site_index]  )->setState( 0 );
                static_cast<CharacterEventDiscrete*>( rightParentState[site_index] )->setState( 1 );
            }
            
//            lnP += log( state_probs[s]/sum );
        }
    }
    // forward proposal is "contract"
    else if ( proposal_type == "contract" ) {
        
//        std::set<size_t>::iterator it_s;
        
//        for (it_s = contractCharacters.begin(); it_s != contractCharacters.end(); it_s++)
        for (size_t i = 0; i < contractCharacters.size(); i++)
        {
            // get ancestral and trunk states
//            size_t site_index = *it_s;
            size_t site_index = contractCharacters[i];

            static_cast<CharacterEventDiscrete*>( nodeChildState[site_index]   )->setState( 0 );
            static_cast<CharacterEventDiscrete*>( leftParentState[site_index]  )->setState( 0 );
            static_cast<CharacterEventDiscrete*>( rightParentState[site_index] )->setState( 0 );
        }

    }
    
    if (perform_swap)
    {
        size_t first_site_index = swapCharacters[0];
        size_t second_site_index = swapCharacters[1];
        
        size_t first_state = static_cast<CharacterEventDiscrete*>( leftParentState[first_site_index] )->getState();
        size_t not_first_state = ( first_state == 1 ? 0 : 1 );
        
        static_cast<CharacterEventDiscrete*>( leftParentState[first_site_index] )->setState( not_first_state );
        static_cast<CharacterEventDiscrete*>( rightParentState[first_site_index] )->setState( first_state );

        static_cast<CharacterEventDiscrete*>( leftParentState[second_site_index] )->setState( first_state );
        static_cast<CharacterEventDiscrete*>( rightParentState[second_site_index] )->setState( not_first_state );
    }
 
    return 0.0;
    
}

template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleSympatryCharacters()
{
    double lnP = 0.0;
    
    // get history object
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    // get phylo coordinates
    TopologyNode &left_child  = node->getChild(0);
    TopologyNode &right_child = node->getChild(1);
    
    size_t node_index  = node->getIndex();
    size_t left_index  = left_child.getIndex();
    size_t right_index = right_child.getIndex();
    
    // learn left/right states
    const std::vector<CharacterEvent*>& leftParentState   = histories[ left_index  ]->getParentCharacters();
    const std::vector<CharacterEvent*>& rightParentState  = histories[ right_index  ]->getParentCharacters();
    std::string trunk_branch = "left";
    size_t trunk_index = left_index;
    size_t bud_index = right_index;
    size_t n_left_on = 0;
    size_t n_right_on = 0;
    for (size_t i = 0; i < leftParentState.size(); i++) {
        n_left_on += static_cast< CharacterEventDiscrete* >( leftParentState[i] )->getState();
        n_right_on += static_cast< CharacterEventDiscrete* >( rightParentState[i] )->getState();
        
        if (n_left_on > 1) {
            break;
        }
        if (n_right_on > 1) {
            trunk_index = right_index;
            bud_index = left_index;
            trunk_branch = "right";
            break;
        }
    }
    
    // update
    std::vector<CharacterEvent*>& trunkParentState        = histories[ trunk_index  ]->getParentCharacters();
    std::vector<CharacterEvent*>& budParentState          = histories[ bud_index    ]->getParentCharacters();
    std::vector<CharacterEvent*>& nodeChildState          = histories[ node_index   ]->getChildCharacters();
   
//    std::cout << "BEFORE " << proposal_type << "\n";
//    printCladogeneticState();
    
//    bool do_swap = true;
    
    if ( proposal_type == "contract" ) {

        // get random set of characters to contract (remove)
//        std::vector<size_t> contractCharacters( sampledCharacters.begin(), sampledCharacters.end() );
        deprecated::random_shuffle( contractCharacters.begin(), contractCharacters.end());
        
        // remove areas in contractCharacters
        for (size_t i = 0; i < contractCharacters.size(); i++) {
            size_t site_index = contractCharacters[i];
            static_cast<CharacterEventDiscrete*>( nodeChildState[site_index] )->setState( 0 );
            static_cast<CharacterEventDiscrete*>( trunkParentState[site_index] )->setState( 0 );
            static_cast<CharacterEventDiscrete*>( budParentState[site_index] )->setState( 0 );
        }
        
        // if we are losing the last character
        if ( onCharacters.size() == 2 ) {
            
//            deprecated::random_shuffle( onCharacters.begin(), onCharacters.end());
            
            size_t site_index = contractCharacters[0];
            static_cast<CharacterEventDiscrete*>( nodeChildState[site_index] )->setState( 1 );
            static_cast<CharacterEventDiscrete*>( trunkParentState[site_index] )->setState( 1 );
            static_cast<CharacterEventDiscrete*>( budParentState[site_index] )->setState( 1 );
            
//            do_swap = false;
        }
//        else {
        
//        }
    }
    else if ( proposal_type == "expand" ) {
//         Sample cladogenetic states for trunk
        // get random set of characters to contract (remove)
//        std::vector<size_t> expandCharacters( sampledCharacters.begin(), sampledCharacters.end() );
        deprecated::random_shuffle( expandCharacters.begin(), expandCharacters.end());
        
        std::set<size_t>::iterator it_s;
        for (size_t i = 0; i < expandCharacters.size(); i++)
        {
            // get ancestral and trunk states
            size_t site_index = expandCharacters[i];
//            size_t A_i = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();
//            size_t T_j = static_cast<CharacterEventDiscrete*>(trunkChildState[site_index])->getState();
//            size_t B_j = static_cast<CharacterEventDiscrete*>(budChildState[site_index])->getState();
//            
            // record the new state and its probability
            static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState( 1 );
            static_cast<CharacterEventDiscrete*>(trunkParentState[site_index])->setState( 1 );
            static_cast<CharacterEventDiscrete*>(budParentState[site_index])->setState( 0 );
        }
        
        if ( onCharacters.size() == 1 ) {

            size_t site_index = expandCharacters[0];
            static_cast<CharacterEventDiscrete*>( nodeChildState[site_index] )->setState( 1 );
            static_cast<CharacterEventDiscrete*>( trunkParentState[site_index] )->setState( 1 );
            static_cast<CharacterEventDiscrete*>( budParentState[site_index] )->setState( 1 );

//            do_swap = false;
        }
    }
    
    
    if (perform_swap)
    {
        
        size_t first_site_index = swapCharacters[0];
        size_t second_site_index = swapCharacters[1];
        
        size_t first_state = static_cast<CharacterEventDiscrete*>( budParentState[first_site_index] )->getState();
        size_t second_state = static_cast<CharacterEventDiscrete*>( budParentState[second_site_index] )->getState();
        
        
        if (first_state == second_state) {
            
            ;
        }
        if (first_site_index == second_site_index) {
            
            ;
        }
        
        static_cast<CharacterEventDiscrete*>( budParentState[first_site_index] )->setState( second_state );
        static_cast<CharacterEventDiscrete*>( budParentState[second_site_index] )->setState( first_state );
        
//        static_cast<CharacterEventDiscrete*>( rightParentState[right_site_index] )->setState( left_state );
//        static_cast<CharacterEventDiscrete*>( rightParentState[left_site_index] )->setState( right_state );
    }
    
//    else {
//        size_t first_site_index = onCharacters[0];
//        size_t second_site_index = swapCharacters[1];
//        
//        size_t first_state = static_cast<CharacterEventDiscrete*>( budParentState[first_site_index] )->getState();
//        size_t second_state = static_cast<CharacterEventDiscrete*>( budParentState[second_site_index] )->getState();
//
//        static_cast<CharacterEventDiscrete*>( leftParentState[first_site_index] )->setState( second_state );
//        static_cast<CharacterEventDiscrete*>( rightParentState[second_site_index] )->setState( right_state );
//        static_cast<CharacterEventDiscrete*>( nodeChildState[ first_site_index ] )->setState( second_state );
//        
//        static_cast<CharacterEventDiscrete*>( leftParentState[second_site_index] )->setState( first_state );
//        
//        static_cast<CharacterEventDiscrete*>( rightParentState[first_site_index] )->setState( left_state );
//
//    }
    

//    std::cout << "AFTER  " << proposal_type << "\n";
//    printCladogeneticState();
    
    return lnP;
    
}

template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleJumpDispersalCharacters()
{
    double lnP = 0.0;
    
    return lnP;
}

template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleRootCharacters( void )
{
    // only apply move if this is an internal node
    if ( !node->isRoot() )
    {
        return 0.0;
    }
    
    // verify tree history object
    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    // states to update
    size_t node_index  = node->getIndex();
    std::vector<CharacterEvent*>& nodeParentState  = histories[ node_index ]->getParentCharacters();

    // sample states
//    std::set<size_t>::iterator it_s;
//    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
//    {
//        size_t site_index = *it_s;
//        size_t s = 0;
//        if (GLOBAL_RNG->uniform01() < 0.5) {
//            s = 1;
//        }
//        static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
//    }
    
    if ( proposal_type == "contract" ) {
        for (size_t i = 0; i < contractCharacters.size(); i++) {
            size_t s = 0;
            if (GLOBAL_RNG->uniform01() < 0.5) {
                s = 1;
            }

            static_cast<CharacterEventDiscrete*>(nodeParentState[i])->setState(s);
        }
    }
    else {
        for (size_t i = 0; i < expandCharacters.size(); i++) {
            size_t s = 0;
            if (GLOBAL_RNG->uniform01() < 0.5) {
                s = 1;
            }
            static_cast<CharacterEventDiscrete*>(nodeParentState[i])->setState(s);
        }
    }

    return 0.0;
}


template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleNodeCharacters( void )
{
    
    // only apply move if this is an internal node
    if ( node->isTip() )
    {
        return 0.0;
    }
    
    //    std::cout << "node_index = " << node->getIndex() << "\n";
    //    std::cout << "clado_state = " << storedCladogeneticEvent << "\n";
    //    printCladogeneticState();
    
    // compute reverse proposal probability
    double ln_bwd_prob = computeNodeLnProposal( storedCladogeneticEvent );
    
    // sample the cladogenetic event type
//    proposedCladogeneticEvent = sampleCladogenesisEventType();
    
    // do the proposal
    double ln_fwd_prob = 0.0;
    if (proposedCladogeneticEvent == "s") {
        ln_fwd_prob = sampleSympatryCharacters();
    } else if (proposedCladogeneticEvent == "a") {
        ln_fwd_prob = sampleAllopatryCharacters();
    } else if (proposedCladogeneticEvent == "j") {
        ln_fwd_prob = sampleJumpDispersalCharacters();
    } else {
        
        throw RbException("Unexpected proposedCladogeneticEvent type!");
    }
    
    
    //    std::cout << "node_index = " << node->getIndex() << "\n";
    //    std::cout << "clado_state = " << proposedCladogeneticEvent << "\n";
    //    printCladogeneticState();
    
    double lnP = ln_bwd_prob - ln_fwd_prob;
    
    return lnP;
    
}


/*
template<class charType>
double RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleExpandContractProposal(void)
{
    
    size_t n_on = onCharacters.size();
    size_t n_off = offCharacters.size();
    
    
    double p_fwd = 0.5;
    double p_bwd = 0.5;
    if (n_on == 1) {
        proposal_type = "expand";
        p_fwd = 1.0;
    }
    else if (n_off == 0) {
        proposal_type = "contract";
        p_fwd = 1.0;
    }
    else if (GLOBAL_RNG->uniform01() < 0.5) {
        proposal_type = "expand";
    }
    else {
        proposal_type = "contract";
    }
    
    // reverse move must be contract
    if (n_on == (numCharacters-1) && proposal_type == "expand") {
        p_bwd = 1.0;
    }
    
    // reverse move must be expand
    if (n_on == 2 && proposal_type == "contract") {
        p_bwd = 1.0;
    }
    
    // cladogenetic event type typically remains constant
    // unless edge case encountered
    proposedCladogeneticEvent = storedCladogeneticEvent;
    
    
    // cladogenetic transitions involving n_on==1
    
    bool expansion_from_1 = (n_on == 1 && proposal_type == "expand");
    bool contraction_from_2 = (n_on == 2 && proposal_type == "contract");
    
    if ( expansion_from_1 || contraction_from_2 ) {

        // get tree history
        GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
        if ( c == NULL )
        {
            throw RbException("Failed cast.");
        }
        
        // if range expands from size 1
        if ( expansion_from_1 ) {
            
            std::string clado_type = "";
            size_t n_clado = cladogeneticSamplingProbabilities.size();
            double p_clado = 1.0 / n_clado;
            double u_type = GLOBAL_RNG->uniform01() * n_clado;
            std::map<std::string, double>::iterator it;
            for (it = cladogeneticSamplingProbabilities.begin(); it != cladogeneticSamplingProbabilities.end(); it++)
            {
                u_type -= 1; //it->second;
                if (u_type <= 0.0) {
                    clado_type = it->first;
                    break;
                }
            }
            
            // assign new clado state
            c->setCladogeneticEvent( clado_type, node->getIndex() );
            proposedCladogeneticEvent = clado_type;

            // collect forward prob of choosing clado event
            p_fwd *= p_clado;
            
        }
        // if range contracts from size 2
        else if (contraction_from_2 ) {
            
            // get bwd expansion prob
            double p_clado = 1.0 / cladogeneticSamplingProbabilities.size(); //[ c->getCladogeneticEvent( node->getIndex() ) ];
            
            // only one possible cladogenetic event type
            proposedCladogeneticEvent = "s";
            
            // collect backward prob
            p_bwd *= p_clado;
        }
    }
    

    // account for probability of selecting random gain/lost area
    if (proposal_type == "expand") {
        p_fwd *= (1.0 / (n_off));
        p_bwd *= (1.0 / (n_on+1));
    } else if (proposal_type == "contract") {
        p_fwd *= (1.0 / (n_on));
        p_bwd *= (1.0 / (n_off+1));
    }
    
    
    if (n_on == 2)
    {
        
        ;
    }
    // return MH ratio for expand/contract
    return log( p_bwd / p_fwd );
    
}
*/

template<class charType>
std::string RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::sampleCladogenesisEventType(void)
{
    
    
    std::string clado_type = "";
    double u_type = GLOBAL_RNG->uniform01();
    std::map<std::string, double>::iterator it;
    for (it = cladogeneticSamplingProbabilities.begin(); it != cladogeneticSamplingProbabilities.end(); it++)
    {
        u_type -= it->second;
        if (u_type <= 0.0) {
            clado_type = it->first;
            break;
        }
    }
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    c->setCladogeneticEvent( clado_type, node->getIndex() );
    
    
    
    //    proposedCladogeneticEvent
    
    return clado_type;
}

template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::setProposalTuningParameter(double tp)
{
    lambda = tp;
}

template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::setSampledCharacters(const std::set<size_t>& s)
{
    sampledCharacters = s;
}


template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
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
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
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



template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::setSwapCharacters( void ) {
    
    swapCharacters.clear();
    leftCharacters.clear();
    rightCharacters.clear();
    
    n_left_fwd = 0;
    n_right_fwd = 0;
    n_left_bwd = 0;
    n_right_bwd = 0;
    
    if (node->isTip()) {
        throw RbException("Unexpected tip node encountered!");
    }
    
    // find swap characters for two easy special cases
    if (onCharacters.size() == 1) {
        swapCharacters.push_back( onCharacters[0] );
        return;
    }
    else if (onCharacters.size() == 2) {
        swapCharacters.push_back( onCharacters[0] );
        swapCharacters.push_back( onCharacters[1] );
        deprecated::random_shuffle( swapCharacters.begin(), swapCharacters.end());
        return;
    }
    
    // get history object
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // gather branch history variables
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    // collect node,left,right states
    size_t left_index = node->getChild(0).getIndex();
    size_t right_index = node->getChild(1).getIndex();
    const std::vector<CharacterEvent*>& leftParentState   = histories[ left_index  ]->getParentCharacters();
    const std::vector<CharacterEvent*>& rightParentState  = histories[ right_index  ]->getParentCharacters();
    
    // learn bud,trunk states
    for (size_t i = 0; i < leftParentState.size(); i++) {
        size_t s_left = static_cast< CharacterEventDiscrete* >( leftParentState[i] )->getState();
        size_t s_right = static_cast< CharacterEventDiscrete* >( rightParentState[i] )->getState();
        
        if (s_left == 1) {
            leftCharacters.push_back( i );
        }
        if (s_right == 1) {
            rightCharacters.push_back( i );
        }
    }
    
    
    // get different left/right states
    size_t s_left = 0;
    size_t s_right = 0;
    do {
        deprecated::random_shuffle( leftCharacters.begin(), leftCharacters.end());
        deprecated::random_shuffle( rightCharacters.begin(), rightCharacters.end());
        
        s_left = leftCharacters[0];
        s_right = rightCharacters[0];
        
    } while(s_left == s_right);
    
    swapCharacters.push_back( s_left );
    swapCharacters.push_back( s_right );
    deprecated::random_shuffle( swapCharacters.begin(), swapCharacters.end());
    
    return;
}


/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
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


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType>
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::tune( double rate )
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
void RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<charType>::undoProposal( void )
{
    
    GeneralTreeHistoryCtmc<charType>* c = dynamic_cast< GeneralTreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }
    size_t num_sites = c->getNumberOfSites();
    
    const std::vector<BranchHistory*>& histories = c->getHistories();
    
    // restore node state
    std::vector<CharacterEvent*>& nodeChildState   = histories[node->getIndex()]->getChildCharacters();
    std::vector<CharacterEvent*>& leftParentState  = histories[node->getChild(0).getIndex() ]->getParentCharacters();
    std::vector<CharacterEvent*>& rightParentState = histories[node->getChild(1).getIndex()]->getParentCharacters();
    
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
        std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()]->getParentCharacters();
        
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

#endif /* defined(__rb_mlandis__BiogeographicCladogeneticRejectionShiftProposal__) */
