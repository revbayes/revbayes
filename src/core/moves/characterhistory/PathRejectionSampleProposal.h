#ifndef PathRejectionSampleProposal_H
#define PathRejectionSampleProposal_H

#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "DeterministicNode.h"
#include "HomologousDiscreteCharacterData.h"
#include "DistributionBinomial.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateGeneratorSequence.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TreeChangeEventMessage.h"
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
    class PathRejectionSampleProposal : public Proposal {
    
    template<class ct>
    friend class NarrowExchangeCharacterHistoryProposal;

    template<class ct>
    friend class FixedNodeheightPruneAndRegraftCharacterHistoryProposal;

    template<class ct>
    friend class NodeTimeSlideUniformCharacterHistoryProposal;

    public:
        PathRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l=1.0, double r=0.234);   //!<  constructor

        // Basic utility functions
        void                                                        assignNode(TopologyNode* nd);
        virtual double                                              computeLnProposal(const TopologyNode& nd, const BranchHistory& bh);
        void                                                        cleanProposal(void);
        PathRejectionSampleProposal*                                clone(void) const;                                                              //!< Clone object
        double                                                      doProposal(void);                                                               //!< Perform proposal
        virtual const std::string&                                  getProposalName(void) const;                                                    //!< Get the name of the proposal for summary printing
        double                                                      getProposalTuningParameter(void) const;
        double                                                      getRootBranchLength(void);                                     //!< get the length of the root branch
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;                                   //!< Print the parameter summary
        void                                                        prepareProposal(void);                                                          //!< Prepare the proposal
        std::set<size_t>                                            sampleCharacters(double p);
        void                                                        setSampledCharacters(const std::set<size_t>& s);
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);                         //!< Set the rate generator.
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d);                 //!< Set the rate generator.
        void                                                        setProposalTuningParameter(double tp);
        void                                                        tune(double r);                                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                        undoProposal(void);                                                             //!< Reject the proposal

    protected:

        void                                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                 //!< Swap the DAG nodes on which the Proposal is working on
        void                                                        fillStateCounts(std::vector<CharacterEvent*> s, std::vector<size_t> &counts);
        double                                                      getBranchRate(size_t index) const;

        // parameters
        StochasticNode<AbstractHomologousDiscreteCharacterData>*    ctmc;
        const TypedDagNode<RateGenerator>*                          q_map_site;
        const TypedDagNode<RateGeneratorSequence>*                  q_map_sequence;

        std::multiset<CharacterEvent*,CharacterEventCompare>        stored_history;

        const TopologyNode*                                         node;

        double                                                      storedLnProb;
        double                                                      proposedLnProb;

        size_t                                                      numStates;
        size_t                                                      numCharacters;

        bool                                                        node_assigned;
        bool                                                        sampled_characters_assigned;
        double                                                      lambda;
        std::set<size_t>                                            sampledCharacters;
        std::set<size_t>                                            allCharacters;

    };
}


#include "DistributionExponential.h"
#include "TreeHistoryCtmc.h"


/**
 * Constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::PathRejectionSampleProposal<charType>::PathRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l, double r) : Proposal(r),
    ctmc(n),
    q_map_site( NULL ),
    q_map_sequence( NULL ),
    node(NULL),
    numCharacters(n->getValue().getNumberOfCharacters()),
    node_assigned(false),
    sampled_characters_assigned(false),
    lambda(l)
{

    addNode(ctmc);
    
    for (size_t i = 0; i < numCharacters; i++)
    {
        allCharacters.insert(i);
    }

}


template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::assignNode(TopologyNode* nd)
{
    node = nd;
    node_assigned = true;
}


template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::cleanProposal( void )
{

//    if ( node->isRoot() && !useTail )
//    {
//        return;
//    }
    
    // delete old events
    std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h;
    std::vector<CharacterEvent*> events;
    for (it_h = stored_history.rbegin(); it_h != stored_history.rend(); ++it_h)
    {
        if (lambda == 1.0)
        {
            events.push_back( *it_h );
        }
        else if (sampledCharacters.find( (*it_h)->getSiteIndex() ) != sampledCharacters.end())
        {
            events.push_back( *it_h );
        }
    }
    for ( size_t i=0; i<events.size(); ++i )
    {
        CharacterEvent* e = events[i];
        delete e;
    }
    stored_history.clear();
    sampledCharacters.clear();
    
    sampled_characters_assigned = false;
    node_assigned = false;
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::PathRejectionSampleProposal<charType>* RevBayesCore::PathRejectionSampleProposal<charType>::clone( void ) const
{
    return new PathRejectionSampleProposal( *this );
}


template<class charType>
double RevBayesCore::PathRejectionSampleProposal<charType>::computeLnProposal(const TopologyNode& nd, const BranchHistory& bh)
{
    
//    if ( nd.isRoot() && !useTail )
//    {
//        return 0.0;
//    }
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }

    double lnP = 0.0;
    
    std::vector<CharacterEvent*> currState = bh.getParentCharacters();
    const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh.getHistory();
    std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h;

    std::vector<size_t> counts(numStates,0);
    fillStateCounts(currState, counts);
    
    // get model parameters
    double currAge = 0.0;

    if ( nd.isRoot() )
    {
        currAge = nd.getAge() + p->getRootBranchLength();
    }
    else
    {
        currAge = nd.getParent().getAge();
    }

    // get sampling ratemap
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );

    // stepwise events
    double dt;
    double eventAge;
    double endAge = nd.getAge();
    double branchRate = getBranchRate(nd.getIndex());
    
    for (it_h = history.rbegin(); it_h != history.rend(); ++it_h)
    {
        // next event time
        double idx = (*it_h)->getSiteIndex();
        eventAge = (*it_h)->getAge();
        dt = currAge - eventAge;
        
        // get the new transition rate
        double tr = rm.getRate( static_cast<CharacterEventDiscrete*>(currState[ (*it_h)->getSiteIndex() ])->getState(), static_cast<CharacterEventDiscrete*>(*it_h)->getState(), currAge, branchRate);
        double sr = rm.getSumOfRates(currState, counts, currAge, branchRate);

        // lnP for stepwise events for p(x->y)
        lnP += log(tr) - (sr * dt);

        // update counts
        counts[ static_cast<CharacterEventDiscrete*>(currState[idx])->getState() ] -= 1;
        counts[ static_cast<CharacterEventDiscrete*>(*it_h)->getState() ] += 1;
        
        // update state
        currState[idx] = *it_h;
        currAge = eventAge;
        
    }
    
    // lnL for final non-event
    double sr = rm.getSumOfRates(currState, counts, currAge, branchRate);
    lnP -= sr * (currAge - endAge);
    
    return lnP;
}




template<class charType>
double RevBayesCore::PathRejectionSampleProposal<charType>::getProposalTuningParameter( void ) const
{
    return lambda;
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
double RevBayesCore::PathRejectionSampleProposal<charType>::doProposal( void )
{
        
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }

    // get model parameters
    double branch_length = node->getBranchLength();
    if ( node->isRoot() )
    {
        branch_length = p->getRootBranchLength();
    }
    
    
    // get sampling rates
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
    
    // rejection sample path history
    BranchHistory* bh = &p->getHistory(*node);
    std::vector<CharacterEvent*> parent_states = bh->getParentCharacters();
    std::vector<CharacterEvent*> child_states  = bh->getChildCharacters();
    std::multiset<CharacterEvent*,CharacterEventCompare> proposed_histories;
        
    // update histories for sites in sampledCharacters
    std::set<size_t>::iterator it_s;
    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
    {
        size_t site_index = *it_s;
        std::set<CharacterEvent*,CharacterEventCompare> tmpHistory;
        size_t currState = static_cast<CharacterEventDiscrete*>(parent_states[site_index])->getState();
        size_t endState  = static_cast<CharacterEventDiscrete*>(child_states[site_index])->getState();

        do
        {
            // delete previously rejected events
            tmpHistory.clear();
            
            // proceed with rejection sampling
            currState = static_cast<CharacterEventDiscrete*>(parent_states[site_index])->getState();
            double end_age = node->getAge();
            double age = end_age + branch_length;
            
            // repeated rejection sampling
            do
            {
                double r = 0.0;
                size_t nextState = 0;
                std::vector<double> rates(numStates,0.0);
                for (size_t i = 0; i < numStates; ++i)
                {
                    if (i == currState)
                    {
                        continue;
                    }
                    double v = rm.getRate(currState, i, age, getBranchRate(node->getIndex()));
                    rates[i] = v;
                    r += v;
                }
                double u = GLOBAL_RNG->uniform01() * r;
                for (size_t i = 0; i < numStates; ++i)
                {
                    u -= rates[i];
                    if (u <= 0.0)
                    {
                        nextState = i;
                        break;
                    }
                }
                // do not force valid time if event needed
                age -= RbStatistics::Exponential::rv(r, *GLOBAL_RNG);
                
                if (age > end_age)
                {
                    currState = nextState;
                    CharacterEvent* evt = new CharacterEventDiscrete(site_index, nextState, age);
                    tmpHistory.insert(evt);

                }
                else if (currState != endState)
                {
                    for (std::set<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h = tmpHistory.rbegin(); it_h != tmpHistory.rend(); it_h++)
                    {
                        delete *it_h;
                    }
                    
                }
            }
            while(age > end_age);
                
        }
        while (currState != endState);
        
        for (std::set<CharacterEvent*,CharacterEventCompare>::iterator it = tmpHistory.begin(); it != tmpHistory.end(); it++)
        {
            proposed_histories.insert(*it);
        }
    }
    
    // assign values back to model for likelihood
    if (lambda == 1.0) 
    {
        bh->setHistory(proposed_histories);
    }
    else
    {
        bh->updateHistory(proposed_histories, sampledCharacters);
    }
    
    // return hastings ratio
    proposedLnProb = computeLnProposal(*node, *bh);
    
    return storedLnProb - proposedLnProb;
}



template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::fillStateCounts(std::vector<CharacterEvent*> s, std::vector<size_t> &counts)
{
    for (size_t i = 0; i < s.size(); ++i)
    {
        counts[ static_cast<CharacterEventDiscrete*>(s[i])->getState() ] += 1;
    }
}


template<class charType>
double RevBayesCore::PathRejectionSampleProposal<charType>::getBranchRate(size_t index) const
{
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }

    double rate = p->getBranchRate(index);

    return rate;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::PathRejectionSampleProposal<charType>::getProposalName( void ) const
{
    static std::string name = "PathRejectionSampleProposal";

    return name;
}

template<class charType>
double RevBayesCore::PathRejectionSampleProposal<charType>::getRootBranchLength( void )
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
void RevBayesCore::PathRejectionSampleProposal<charType>::prepareProposal( void )
{
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // make sure the stored history is properly cleaned (no memory leaks)
    std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h;
    std::vector<CharacterEvent*> old_events;
    for (it_h = stored_history.rbegin(); it_h != stored_history.rend(); ++it_h)
    {
        if (lambda == 1.0)
        {
            old_events.push_back( *it_h );
        }
        else if (sampledCharacters.find( (*it_h)->getSiteIndex() ) != sampledCharacters.end())
        {
            old_events.push_back( *it_h );
        }
    }
    for ( size_t i=0; i<old_events.size(); ++i )
    {
        CharacterEvent* e = old_events[i];
        delete e;
    }
    stored_history.clear();
    
    storedLnProb = 0.0;
    proposedLnProb = 0.0;
    
    // only pick a random node if it wasn't assigned
    if ( node_assigned == false )
    {
        const Tree &tau = p->getTree();
        size_t num_nodes = tau.getNumberOfNodes();
        
        size_t node_index = 0;
        double root_branch_length = getRootBranchLength();
        if ( root_branch_length == 0 )
        {
            do {
            node_index = GLOBAL_RNG->uniform01() * num_nodes;
            node = &tau.getNode(node_index);
            } while ( node->isRoot() );
        }
        else
        {
            node_index = GLOBAL_RNG->uniform01() * num_nodes;
            node = &tau.getNode(node_index);
        }
    }
    
    
    BranchHistory* bh = &p->getHistory(*node);
    //    stored_history = history;
    const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh->getHistory();
    for (it_h = history.rbegin(); it_h != history.rend(); ++it_h)
    {
        if (lambda == 1.0)
        {
            stored_history.insert( (*it_h)->clone() );
        }
        else if (sampledCharacters.find( (*it_h)->getSiteIndex() ) != sampledCharacters.end())
        {
            stored_history.insert( (*it_h)->clone() );
        }
    }
    
    
    // determine sampled characters
    if (!sampled_characters_assigned)
    {
        sampledCharacters = sampleCharacters(lambda);
    }
    
    // flag node as dirty
    const_cast<TopologyNode*>(node)->fireTreeChangeEvent(RevBayesCore::TreeChangeEventMessage::CHARACTER_HISTORY);
    
    // compute backwards proposal probability
    double x = computeLnProposal(*node, *bh);
    storedLnProb = x;
    
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
void RevBayesCore::PathRejectionSampleProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
//    o << "lambda = " << lambda;
}

template<class charType>
std::set<size_t> RevBayesCore::PathRejectionSampleProposal<charType>::sampleCharacters(double p)
{

    if (p == 1.0)
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
void RevBayesCore::PathRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{

    q_map_site = d;
    numStates = q_map_site->getValue().getNumberOfStates();

}


template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{

    q_map_sequence = d;
    numStates = q_map_sequence->getValue().getNumberOfStates();

}


template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::setSampledCharacters(const std::set<size_t>& s)
{
    sampledCharacters = s;
    sampled_characters_assigned = true;
}


/**
 * Swap the current ctmc for a new one.
 *
 * \param[in]     oldN     The old ctmc that needs to be replaced.
 * \param[in]     newN     The new ctmc.
 */
template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if (oldN == ctmc)
    {
        ctmc = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN) ;
    }
    else if (oldN == q_map_site)
    {
        q_map_site = static_cast<TypedDagNode<RateGenerator>* >(newN);
    }
    else if (oldN == q_map_sequence)
    {
        q_map_sequence = static_cast<TypedDagNode<RateGeneratorSequence>* >(newN);
    }

}


template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType>
void RevBayesCore::PathRejectionSampleProposal<charType>::tune( double rate )
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
void RevBayesCore::PathRejectionSampleProposal<charType>::undoProposal( void )
{
    TreeHistoryCtmc<charType>* p = dynamic_cast< TreeHistoryCtmc<charType>* >( &ctmc->getDistribution() );
    if ( p == NULL )
    {
        throw RbException("Failed cast.");
    }
    
    // delete new events
    BranchHistory* bh = &p->getHistory(*node);
    
    std::multiset<CharacterEvent*,CharacterEventCompare> proposed_history = bh->getHistory();
    std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h;
    std::vector<CharacterEvent*> events;
    for (it_h = proposed_history.rbegin(); it_h != proposed_history.rend(); ++it_h)
    {
        if (lambda == 1.0)
        {
            events.push_back( *it_h );
        }
        else if (sampledCharacters.find( (*it_h)->getSiteIndex() ) != sampledCharacters.end())
        {
            events.push_back( *it_h );
        }
    }
    for ( size_t i=0; i<events.size(); ++i )
    {
        CharacterEvent* e = events[i];
        delete e;
    }
    
    // flag node as dirty
    const_cast<TopologyNode*>(node)->fireTreeChangeEvent(RevBayesCore::TreeChangeEventMessage::CHARACTER_HISTORY);
    
    // swap current value and stored value
    bh->setHistory(stored_history);

    // clear old histories
    proposed_history.clear();
    stored_history.clear();
    sampledCharacters.clear();
    
    sampled_characters_assigned = false;
    node_assigned = false;
        
}


#endif
