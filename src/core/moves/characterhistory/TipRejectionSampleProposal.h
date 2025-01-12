#ifndef TipRejectionSampleProposal_H
#define TipRejectionSampleProposal_H

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
#include "StochasticNode.h"
#include "TreeChangeEventMessage.h"
#include "TopologyNode.h"
#include "TypedDagNode.h"

#include <cmath>
#include <iostream>
#include <set>
#include <string>

namespace RevBayesCore {

    /**
     * Resample a tip-state and the branch history.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2025-01-10, version 1.2.6
     *
     */

    template<class charType>
    class TipRejectionSampleProposal : public Proposal {

    public:
        TipRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l=1.0, double r=0.234 );                                  //!<  constructor
        TipRejectionSampleProposal( const TipRejectionSampleProposal& p );                                                        //!<  constructor
        virtual                                                    ~TipRejectionSampleProposal(void);                              //!< Destructor
        
        TipRejectionSampleProposal&                                 operator=(const TipRejectionSampleProposal& p);

        // Basic utility functions
        void                                                        assignNode(TopologyNode* nd);
        void                                                        assignSiteIndexSet(const std::set<size_t>& s);
        TipRejectionSampleProposal*                                 clone(void) const;                                              //!< Clone object
        void                                                        cleanProposal(void);
//        virtual double                                              computeLnProposal();
        double                                                      doProposal(void);                                               //!< Perform proposal
        const std::string&                                          getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                                      getProposalTuningParameter(void) const;
        void                                                        printParameterSummary(std::ostream &o, bool name_only) const;                   //!< Print the parameter summary
        void                                                        prepareProposal(void);                                          //!< Prepare the proposal
        std::set<size_t>                                            chooseCharactersToSample(double p);
        void                                                        setSampledCharacters(const std::set<size_t>& s);
        void                                                        sampleTipCharacters(void);                                     //!< Sample the characters at the node
//        double                                                      sampleRootCharacters(void);                                     //!< Sample the characters at the root
        void                                                        setRateGenerator(const TypedDagNode<RateGenerator> *d);         //!< Set the rate generator.
        void                                                        setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d); //!< Set the rate generator.
        void                                                        setProposalTuningParameter(double tp);
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
        std::vector<size_t>                                         storedTipState;
//        std::vector<size_t>                                         storedSubrootState;

        TopologyNode*                                               node;
        double                                                      storedLnProb;
        double                                                      proposedLnProb;

        PathRejectionSampleProposal<charType>*                      tip_branch_proposal;

        TransitionProbabilityMatrix                                 tip_tp_matrix;

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
RevBayesCore::TipRejectionSampleProposal<charType>::TipRejectionSampleProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n, double l, double r ) : Proposal(r),
    ctmc(n),
    q_map_site( NULL ),
    q_map_sequence( NULL ),
    numCharacters(n->getValue().getNumberOfCharacters()),
    numStates(2),
    node( NULL ),
    tip_tp_matrix(2),
    lambda(l)
{

    addNode( ctmc );

    tip_branch_proposal  = new PathRejectionSampleProposal<charType>(n, l, r);
    
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
RevBayesCore::TipRejectionSampleProposal<charType>::TipRejectionSampleProposal( const TipRejectionSampleProposal& p ) : Proposal(p),
    ctmc( p.ctmc ),
    q_map_site( p.q_map_site ),
    q_map_sequence( p.q_map_sequence ),
    numCharacters( p.numCharacters ),
    numStates( p.numStates ),
    node( p.node ),
    tip_tp_matrix( p.tip_tp_matrix ),
    lambda( p.lambda )
{

    addNode( ctmc );
        
    tip_branch_proposal = p.tip_branch_proposal->clone();
    
    sampledCharacters   = p.sampledCharacters;
    allCharacters       = p.allCharacters;

    storedTipState     = p.storedTipState;
//    storedSubrootState  = p.storedSubrootState;

    storedLnProb        = p.storedLnProb;
    proposedLnProb      = p.proposedLnProb;

}


/**
  * Destructor
  *
  */
 template<class charType>
 RevBayesCore::TipRejectionSampleProposal<charType>::~TipRejectionSampleProposal( void )
 {
         
     delete tip_branch_proposal;

 }


/**
 * Copy constructor
 *
 * Here we simply allocate and initialize the Proposal object.
 */
template<class charType>
RevBayesCore::TipRejectionSampleProposal<charType>& RevBayesCore::TipRejectionSampleProposal<charType>::operator=( const TipRejectionSampleProposal& p )
{

    if (this != &p)
    {
        delete tip_branch_proposal;
        
        removeNode(ctmc);
        
        ctmc            = p.ctmc;
        q_map_site      = p.q_map_site;
        q_map_sequence  = p.q_map_sequence;
        numCharacters   = p.numCharacters;
        numStates       = p.numStates;
        node            = p.node;
        tip_tp_matrix   = p.tip_tp_matrix;
        lambda          = p.lambda;
    
        addNode( ctmc );
        
        tip_branch_proposal  = p.tip_branch_proposal->clone();
        
        
        sampledCharacters   = p.sampledCharacters;
        allCharacters       = p.allCharacters;

        storedTipState     = p.storedTipState;
//        storedSubrootState  = p.storedSubrootState;

        storedLnProb        = p.storedLnProb;
        proposedLnProb      = p.proposedLnProb;
    }

    
    return *this;
}


template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::assignNode(TopologyNode* nd)
{
    node = nd;
}


template<class charType>
std::set<size_t> RevBayesCore::TipRejectionSampleProposal<charType>::chooseCharactersToSample(double p)
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
void RevBayesCore::TipRejectionSampleProposal<charType>::cleanProposal( void )
{
    tip_branch_proposal->cleanProposal();

    storedTipState.clear();
//    storedSubrootState.clear();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class charType>
RevBayesCore::TipRejectionSampleProposal<charType>* RevBayesCore::TipRejectionSampleProposal<charType>::clone( void ) const
{
    return new TipRejectionSampleProposal( *this );
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
double RevBayesCore::TipRejectionSampleProposal<charType>::doProposal( void )
{
    proposedLnProb = 0.0;

    double proposedLnProbRatio = 0.0;

    // backward proposal probability
//    proposedLnProbRatio += computeLnProposal();

    // update node value
    sampleTipCharacters();

    // forward proposal probability
//    proposedLnProbRatio -= computeLnProposal();

    // update incident path
    proposedLnProbRatio += tip_branch_proposal->doProposal();

    
    return proposedLnProbRatio;
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class charType>
const std::string& RevBayesCore::TipRejectionSampleProposal<charType>::getProposalName( void ) const
{
    static std::string name = "TipRejectionSampleProposal";

    return name;
}


template<class charType>
double RevBayesCore::TipRejectionSampleProposal<charType>::getProposalTuningParameter( void ) const
{
    return lambda;
}


/**
 *
 */
template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::prepareProposal( void )
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
    } while ( node->isTip() == false );

    // prepare the path proposals
    tip_branch_proposal->assignNode(node);
    tip_branch_proposal->prepareProposal();


    // store node state values
    storedTipState.clear();

    size_t num_sites = p->getNumberOfSites();
    storedTipState.resize(num_sites,0);
    const std::vector<CharacterEvent*>& curr_tip_state = p->getHistory(*node).getChildCharacters();
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(curr_tip_state[site_index])->getState();
        storedTipState[site_index] = s;
    }

    
    // sample characters to be updated and pass to proposals
    sampledCharacters = chooseCharactersToSample(lambda);
    tip_branch_proposal->setSampledCharacters(sampledCharacters);

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
void RevBayesCore::TipRejectionSampleProposal<charType>::printParameterSummary(std::ostream &o, bool name_only) const
{
    o << "lambda = ";
    if (name_only == false)
    {
        o << lambda;
    }
}

template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::sampleTipCharacters( void )
{

    TreeHistoryCtmc<charType>* c = dynamic_cast< TreeHistoryCtmc<charType>* >(&ctmc->getDistribution());
    if ( c == NULL )
    {
        throw RbException("Failed cast.");
    }

    CharacterHistoryDiscrete& histories = c->getHistories();

    size_t node_index = node->getIndex();
    double node_age   = node->getAge();
    double node_rate  = c->getBranchRate( node_index );

    // states to update
    std::vector<CharacterEvent*>& nodeChildState   = histories[ node_index ].getChildCharacters();

    // we may also update this if it is the root state (not const)
    std::vector<CharacterEvent*>& nodeParentState = histories[node->getIndex()].getParentCharacters();
    
    // get parent age
    double parent_age = node->getParent().getAge();
    
    // get transition probs
    const RateGenerator& rm = ( q_map_sequence != NULL ? q_map_sequence->getValue() : q_map_site->getValue() );
    rm.calculateTransitionProbabilities(parent_age, node_age, node_rate, tip_tp_matrix);

    
    std::set<size_t>::iterator it_s;
    
    for (it_s = sampledCharacters.begin(); it_s != sampledCharacters.end(); it_s++)
    {
        size_t site_index = *it_s;
        size_t ancS  = static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->getState();

        double u = GLOBAL_RNG->uniform01();

        std::vector<double> state_probs(numStates,0.0);
        double sum = 0.0;
        size_t num_valid_states = 0;
        for ( size_t i=0; i<numStates; ++i )
        {
            if ( static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->isValidState(i) )
            {
                double p = tip_tp_matrix[ancS][i];
                sum += p;
                state_probs[i] = p;
                ++num_valid_states;
            }
        }
        if ( num_valid_states > 2 )
        {
            std::cerr << "Too many valid stated:\t" << num_valid_states << std::endl;
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
    }



}




template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGenerator> *d)
{

    q_map_site = d;
    numStates = q_map_site->getValue().getNumberOfStates();

    tip_branch_proposal->setRateGenerator( q_map_site );

    tip_tp_matrix = TransitionProbabilityMatrix(numStates);

}


template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::setRateGenerator(const TypedDagNode<RateGeneratorSequence> *d)
{

    q_map_sequence = d;
    numStates = q_map_sequence->getValue().getNumberOfStates();

    tip_branch_proposal->setRateGenerator( q_map_sequence );

    tip_tp_matrix = TransitionProbabilityMatrix(numStates);

}


template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::setSampledCharacters(const std::set<size_t>& s)
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
void RevBayesCore::TipRejectionSampleProposal<charType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
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

    tip_branch_proposal->swapNode(oldN, newN);
}


template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::setProposalTuningParameter(double tp)
{
    lambda = tp;
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 */
template<class charType>
void RevBayesCore::TipRejectionSampleProposal<charType>::tune( double rate )
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
void RevBayesCore::TipRejectionSampleProposal<charType>::undoProposal( void )
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
    
    for (size_t site_index = 0; site_index < num_sites; ++site_index)
    {
        size_t s = storedTipState[site_index];
        static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
    }
    
    
    // restore path state
    tip_branch_proposal->undoProposal();
    
    // clear sampled character set
    sampledCharacters.clear();
}

#endif /* defined(__rb_mlandis__TipRejectionSampleProposal__) */
