#ifndef GeneralTreeHistoryCtmcSiteIID_H
#define GeneralTreeHistoryCtmcSiteIID_H

#include "TreeHistoryCtmc.h"

#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "HomologousDiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "RbVector.h"
#include "Simplex.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

#include <cmath>
#include <cstring>

namespace RevBayesCore {

    
    template<class charType>
    class GeneralTreeHistoryCtmcSiteIID : public TreeHistoryCtmc<charType> {

    public:
        GeneralTreeHistoryCtmcSiteIID(const TypedDagNode< Tree > *t, size_t nChars, size_t nSites, bool useAmbigChar=false);
        virtual                                             ~GeneralTreeHistoryCtmcSiteIID(void);                                                //!< Virtual destructor

        // public member functions

        GeneralTreeHistoryCtmcSiteIID*                      clone(void) const;                                                           //!< Create an independent clone
        void                                                initializeTipValues(void);
        bool                                                drawInitValue(void);
        std::vector<double>                                 getRootFrequencies(void) const;
        void                                                redrawValue(void);
        virtual void                                        simulate(void);

        bool                                                samplePathStart(const TopologyNode& node);
        bool                                                samplePathEnd(const TopologyNode& node);
        bool                                                samplePathHistory(const TopologyNode& node);

        void                                                setRateGenerator(const TypedDagNode< RateGenerator > *rm);
        void                                                setRateGenerator(const TypedDagNode< RbVector< RateGenerator > > *rm);
        void                                                setPInv(const TypedDagNode< double > *);
        void                                                setRootFrequencies(const TypedDagNode< Simplex > *f);

    protected:

        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

        virtual double                                      computeRootLikelihood(const TopologyNode &n);
        virtual double                                      computeInternalNodeLikelihood(const TopologyNode &n);
        virtual double                                      computeTipLikelihood(const TopologyNode &node);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);

    private:

        // members
        const TypedDagNode< RateGenerator >*                homogeneousRateGenerator;
        const TypedDagNode< RbVector< RateGenerator > >*    heterogeneousRateGenerator;
        const TypedDagNode< Simplex >*                      rootFrequencies;
        const TypedDagNode< double >*                       pInv;
        

        // flags specifying which model variants we use
        bool                                                branchHeterogeneousSubstitutionMatrices;

        virtual void                                        simulate(const TopologyNode& node, BranchHistory* bh, std::vector< DiscreteTaxonData< charType > >& taxa);
        void                                                simulateHistory(const TopologyNode& node, BranchHistory* bh);
        std::vector<size_t>                                 computeCounts(const std::vector<CharacterEvent*>& s);

    };

}


#include "CharacterEventDiscrete.h"
#include "PathRejectionSampleProposal.h"
#include "RateMatrix_JC.h"

template<class charType>
RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::GeneralTreeHistoryCtmcSiteIID(const TypedDagNode<Tree> *tau, size_t nChars, size_t nSites, bool useAmbigChar) : TreeHistoryCtmc<charType>( tau, nChars, nSites, useAmbigChar )
{

    // initialize with default parameters
    homogeneousRateGenerator    = new ConstantNode<RateGenerator>("rateGenerator", new RateMatrix_JC( nChars ) );
    heterogeneousRateGenerator  = NULL;
    pInv                        = new ConstantNode<double>("pInv", new double(0.0) );
    rootFrequencies             = NULL;
    

    // flags specifying which model variants we use
    branchHeterogeneousSubstitutionMatrices     = false;

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( tau );
    this->addParameter( homogeneousRateGenerator );
    this->addParameter( heterogeneousRateGenerator );
    this->addParameter( rootFrequencies );
    this->addParameter( pInv );
    
}


template<class charType>
RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::~GeneralTreeHistoryCtmcSiteIID( void )
{
    ; // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


template<class charType>
RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>* RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::clone( void ) const
{

    return new GeneralTreeHistoryCtmcSiteIID<charType>( *this );
}

template<class charType>
std::vector<size_t> RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::computeCounts(const std::vector<CharacterEvent*>& s)
{
    std::vector<size_t> counts(this->num_states, 0);

    for (size_t i = 0; i < s.size(); i++)
    {
        counts[ static_cast<CharacterEventDiscrete*>(s[i])->getState() ] += 1;
    }

    return counts;
}

template<class charType>
double RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::computeRootLikelihood(const TopologyNode &n)
{
    double lnP = 0.0;

    const BranchHistory& bh = this->histories[n.getIndex()];
    const std::vector<CharacterEvent*>& rootState = bh.getChildCharacters();

    // if ( c->getRootBranchLength() - node->getAge() != 0 )
    // {
    //
    // }
    
    // get counts per state
    std::vector<int> counts(this->num_states, 0);
    for (size_t i = 0; i < rootState.size(); ++i)
    {
        ++counts[ static_cast<CharacterEventDiscrete*>(rootState[i])->getState() ];
    }

    // get log prob
    // PL comments: change below such that root frequencies can be obtained from stationary distribution when root branch is absent.
    std::vector<double> rf = getRootFrequencies();
    for (size_t i = 0; i < counts.size(); i++)
    {
        // if root branch is present
        lnP += counts[i] * log( rf[i] );
    }

    return lnP;
}

template<class charType>
double RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::computeInternalNodeLikelihood(const TopologyNode &node)
{

    if ( node.isRoot() == true )
    {
        return 0.0;
    }

    size_t node_index = node.getIndex();
    double branch_rate = this->getBranchRate(node_index);
    const RateGenerator& rm = homogeneousRateGenerator->getValue();

    // get the branch history
    BranchHistory& bh = this->histories[node_index];

    // check that node ages are consistent with character event ages
    if ( bh.areEventTimesValid(node) == false )
    {
        return RbConstants::Double::neginf;
    }

    // check parent and child states to make sure they match with the
    // ancestral and descendant branches; otherwise, return -Inf
    std::vector<CharacterEvent*> curr_state = bh.getParentCharacters();
    std::vector<CharacterEvent*> end_state  = bh.getChildCharacters();
    for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
    {
        const TopologyNode &child = node.getChild(i);
        size_t child_index = child.getIndex();
        BranchHistory& child_bh = this->histories[child_index];
        std::vector<CharacterEvent*> child_state = child_bh.getParentCharacters();
        for (size_t j = 0; j < this->num_sites; ++j)
        {
            if ( static_cast<CharacterEventDiscrete*>(end_state[j])->getState() != static_cast<CharacterEventDiscrete*>(child_state[j])->getState() )
            {
                return RbConstants::Double::neginf;
            }
        }
    }

    // we need the counts for faster computation
    std::vector<size_t> counts = computeCounts(curr_state);

    const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh.getHistory();
    std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h;

    // stepwise events
    double lnL = 0.0;
    double current_age = node.getParent().getAge();
    double end_age = node.getAge();
    double event_age;
    
    for (it_h = history.rbegin(); it_h != history.rend(); ++it_h)
    {
        CharacterEventDiscrete* char_event = static_cast<CharacterEventDiscrete*>(*it_h);

        // next event time
        double idx = char_event->getSiteIndex();
        event_age = char_event->getAge();
        size_t s = char_event->getState();

        // lnL for stepwise events for p(x->y)
        double tr = rm.getRate(static_cast<CharacterEventDiscrete*>(curr_state[idx])->getState(), char_event->getState(), current_age, branch_rate);
        double sr = rm.getSumOfRates(curr_state, counts) * branch_rate;
        lnL += log(tr) - sr * (current_age - event_age);
        
        // update counts
        counts[static_cast<CharacterEventDiscrete*>(curr_state[idx])->getState()] -= 1;
        counts[s] += 1;

        // update time and state
        curr_state[idx] = char_event;
        current_age = event_age;
    }

    // lnL that nothing else happens
    double sr = rm.getSumOfRates(curr_state) * branch_rate;
    lnL -= sr * (current_age - end_age);

    return lnL;
}



template<class charType>
double RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::computeTipLikelihood(const TopologyNode &node)
{
    double lnL = 0.0;
    return lnL;

}

template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::drawInitValue( void )
{

    if ( this->tipsInitialized == false )
    {
        initializeTipValues();
    }

    // sample node states
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        TopologyNode* nd = nodes[i];

        int samplePathEndCount = 0;
        do
        {
            samplePathEndCount++;
        } while (samplePathEnd(*nd) == false && samplePathEndCount < 100);

        int samplePathStartCount = 0;
        do
        {
            samplePathStartCount++;
        } while (samplePathStart(*nd) == false && samplePathStartCount < 100);

    }

    // sample paths
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        TopologyNode* nd = nodes[i];

        int samplePathHistoryCount = 0;
        do
        {
            ++samplePathHistoryCount;
        } while (samplePathHistory(*nd) == false && samplePathHistoryCount < 100);

        //        this->histories[i]->print();
    }

    double lnL = this->computeLnProbability();

    if (lnL == RbConstants::Double::neginf)
    {
        for (size_t i = 0; i < nodes.size(); i++)
        {
            this->fireTreeChangeEvent(*nodes[i]);
        }
        return false;
    }

    return true;
}

template<class charType>
std::vector<double> RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::getRootFrequencies( void ) const
{
    // PL comments: look at here for how to call the analytical function and root frequencies
    
    if ( branchHeterogeneousSubstitutionMatrices == true || rootFrequencies != NULL )
    {
        return rootFrequencies->getValue();
    }
    else
    {
        const RateMatrix *rm = dynamic_cast<const RateMatrix *>( &homogeneousRateGenerator->getValue() );
        if ( rm != NULL )
        {
            return rm->getStationaryFrequencies();
        }
        else
        {
            throw RbException("You either need to use a rate-matrix or specify root frequencies.");
        }
    }

}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::initializeTipValues( void )
{

    {
        std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            TopologyNode* node = nodes[i];
            if ( node->isTip() == true )
            {
                DiscreteTaxonData<charType>& d = static_cast< DiscreteTaxonData<charType>& >( this->value->getTaxonData( node->getName() ) );

                std::vector<CharacterEvent*> tipState;
                for (size_t j = 0; j < d.getNumberOfCharacters(); ++j)
                {
                    DiscreteCharacterState &state = d[j];
                    unsigned s = 0;
                    
                    if ( state.isGapState() == true || state.isMissingState() == true )
                    {
                        s = 0;
                        double u = GLOBAL_RNG->uniform01();
                        std::vector<double> state_probs(this->num_states, 1.0/this->num_states);
                        for ( size_t i=0; i<this->num_states; ++i )
                        {
                            u -= state_probs[i];
                            if ( u <= 0.0 )
                            {
                                break;
                            }
                            ++s;
                        }
                    }
                    else if ( state.isAmbiguous() == true )
                    {
                        RbBitSet obs_states = state.getState();
                        std::vector<double> state_probs(this->num_states, 0.0);
                        for ( size_t i=0; i<this->num_states; ++i )
                        {
                            if ( obs_states.test(i) )
                            {
                                state_probs[i] = 1.0/obs_states.count();
                            }
                        }
                        
                        double u = GLOBAL_RNG->uniform01();
                        s = 0;
                        for ( size_t i=0; i<this->num_states; ++i )
                        {
                            u -= state_probs[i];
                            if ( u <= 0.0 )
                            {
                                break;
                            }
                            ++s;
                        }
                    }
                    else
                    {
                        s = (unsigned) state.getStateIndex();
                    }
                    
                    CharacterEventDiscrete* evt = new CharacterEventDiscrete(j, s, 1.0);
                    
                    if ( state.isGapState() == true || state.isMissingState() == true )
                    {
                        evt->setMissingState( true );
                    }
                    else if ( state.isAmbiguous() == true )
                    {
                        RbBitSet obs_states = state.getState();
                        evt->setAmbiguousState( obs_states );
                    }
                    tipState.push_back( evt );
                }

                this->histories[node->getIndex()].setChildCharacters(tipState);
            }
        }

        this->tipsInitialized = true;
    }

}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::redrawValue( void )
{

    simulate();

}



template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::samplePathEnd(const TopologyNode& node)
{
    if ( node.isTip() == true )
    {
        return true; // do nothing
    }
    else
    {
        TransitionProbabilityMatrix leftTpMatrix(this->num_states);
        TransitionProbabilityMatrix rightTpMatrix(this->num_states);
        //        TransitionProbabilityGenerator ancTpGenerator(this->num_states);

        const RateGenerator& rm = homogeneousRateGenerator->getValue();

        double begin_age = node.getAge();

        size_t left_index  = node.getChild(0).getIndex();
        size_t right_index = node.getChild(1).getIndex();

        double left_branch_rate  = this->getBranchRate( left_index );
        double right_branch_rate = this->getBranchRate( right_index );

        rm.calculateTransitionProbabilities(begin_age, node.getChild(0).getAge(), left_branch_rate,  leftTpMatrix);
        rm.calculateTransitionProbabilities(begin_age, node.getChild(0).getAge(), right_branch_rate, rightTpMatrix);

        // for sampling probs
        const std::vector<CharacterEvent*>& leftChildState  = this->histories[left_index].getChildCharacters();
        const std::vector<CharacterEvent*>& rightChildState = this->histories[right_index].getChildCharacters();

        // to update
        std::vector<CharacterEvent*> nodeChildState = this->histories[node.getIndex()].getChildCharacters();
        for (size_t site_index=0; site_index<this->num_sites; ++site_index)
        {

            size_t desS1 = static_cast<CharacterEventDiscrete*>(leftChildState[site_index])->getState();
            size_t desS2 = static_cast<CharacterEventDiscrete*>(rightChildState[site_index])->getState();

            std::vector<double> state_probs(this->num_states, 0.0);
            double prob_sum = 0.0;
            for (size_t i = 0; i < this->num_states; ++i)
            {
                state_probs[i] = leftTpMatrix[i][desS1] * rightTpMatrix[i][desS2];
                prob_sum += state_probs[i];
            }

            double u = GLOBAL_RNG->uniform01() * prob_sum;
            unsigned int s = 0;
            for (unsigned int i = 0; i < this->num_states; ++i)
            {
                u -= state_probs[i];
                if (u <= 0.0)
                {
                    s = i;
                    break;
                }
            }

            static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->setState(s);
        }

    }

    return true;
}


template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::samplePathHistory(const TopologyNode& node)
{

    if ( node.isRoot() == true )
    {
        return true;
    }

    PathRejectionSampleProposal<charType> p( this->getStochasticNode() );
    p.setRateGenerator( homogeneousRateGenerator );

    p.assignNode(const_cast<TopologyNode*>(&node));

    p.prepareProposal();
    p.doProposal();
    p.cleanProposal();

    return true;
}


template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::samplePathStart(const TopologyNode& node)
{

    // ignore tips
    if ( node.isTip() == true )
    {
        return true;
    }

    // update daughter lineages' parent states with the current node's child state
    const std::vector<CharacterEvent*>& nodeChildState = this->histories[ node.getIndex() ].getChildCharacters();
    const std::vector<TopologyNode*>& children = node.getChildren();
    for (size_t i = 0; i < children.size(); ++i)
    {
        std::vector<CharacterEvent*> childParentState(nodeChildState.size());
        for (size_t j = 0; j < childParentState.size(); ++j)
        {
            childParentState[j] = new CharacterEventDiscrete(j, static_cast<CharacterEventDiscrete*>(nodeChildState[j])->getState(), 0.0);
        }
        this->histories[ children[i]->getIndex() ].setParentCharacters( childParentState );
    }

    return true;
}



template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::setPInv(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( pInv != NULL )
    {
        this->removeParameter( pInv );
        pInv = NULL;
    }

    // set the value
    pInv = r;

    // add the new parameter
    this->addParameter( pInv );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::setRateGenerator(const TypedDagNode< RateGenerator > *rm) {

    // remove the old parameter first
    if ( homogeneousRateGenerator != NULL )
    {
        this->removeParameter( homogeneousRateGenerator );
        homogeneousRateGenerator = NULL;
    }
    else // heterogeneousRateGenerator != NULL
    {
        this->removeParameter( heterogeneousRateGenerator );
        heterogeneousRateGenerator = NULL;
    }

    // set the value
    branchHeterogeneousSubstitutionMatrices = false;
    homogeneousRateGenerator = rm;

    // add the new parameter
    this->addParameter( homogeneousRateGenerator );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::setRateGenerator(const TypedDagNode< RbVector< RateGenerator > > *rm) {

    // remove the old parameter first
    if ( homogeneousRateGenerator != NULL )
    {
        this->removeParameter( homogeneousRateGenerator );
        homogeneousRateGenerator = NULL;
    }
    else // heterogeneousRateGenerator != NULL
    {
        this->removeParameter( heterogeneousRateGenerator );
        heterogeneousRateGenerator = NULL;
    }

    // set the value
    branchHeterogeneousSubstitutionMatrices = true;
    heterogeneousRateGenerator = rm;

    // add the new parameter
    this->addParameter( heterogeneousRateGenerator );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::setRootFrequencies(const TypedDagNode< Simplex > *f)
{

    // remove the old parameter first
    if ( rootFrequencies != NULL )
    {
        this->removeParameter( rootFrequencies );
        rootFrequencies = NULL;
    }

    if ( f != NULL )
    {
        // set the value
        rootFrequencies = f;
    }
    else
    {
        branchHeterogeneousSubstitutionMatrices = false;
    }

    // add the new parameter
    this->addParameter( rootFrequencies );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}



template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::simulate(void)
{
    TreeHistoryCtmc<charType>::simulate();
}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::simulateHistory(const TopologyNode& node, BranchHistory* bh)
{
    
    size_t branch_index = node.getIndex();
    double branch_rate = this->getBranchRate( branch_index );
    const RateGenerator& rm = homogeneousRateGenerator->getValue();

    // get the start and end times
    double branch_length = node.getBranchLength();
    double end_age = node.getAge();
    double start_age = end_age + branch_length;
    
    // get parent BranchHistory state
    std::vector<CharacterEvent*> currState = bh->getParentCharacters();
    std::vector<size_t> counts = computeCounts(currState);
    std::set<CharacterEvent*,CharacterEventCompare> history;

    // simulate path
    double t = start_age;
    double dt = 0.0;
    while (t - dt > end_age)
    {
        
        // sample next event time
        double sr = rm.getSumOfRates(currState, counts) * branch_rate;
        dt = RbStatistics::Exponential::rv(sr, *GLOBAL_RNG);
        if (t - dt > end_age)
        {
            // next event type
            CharacterEventDiscrete* evt = new CharacterEventDiscrete(0, 0, t - dt);
            double u = GLOBAL_RNG->uniform01() * sr;

            bool found = false;
            size_t i, s = 0;
            for (i = 0; !found && i < this->num_sites; ++i)
            {
                evt->setSiteIndex(i);
                for (s = 0; !found && s < this->num_states; ++s)
                {
                    // disregard virtual events (self-transitions)
                    if (s != static_cast<CharacterEventDiscrete*>(currState[i])->getState())
                    {
                        evt->setState(s);
                        //                        double r = rm.getRate(currState, evt, counts);
                        double r = rm.getRate(static_cast<CharacterEventDiscrete*>(currState[i])->getState(), evt->getState(), node.getAge(), branch_rate);

                        u -= r;
                        if (u <= 0.0)
                        {
                            found = true;
                            history.insert(evt);
                        }
                    }
                    if (found) break;
                }
                if (found) break;
            }

            // update counts
            counts[ static_cast<CharacterEventDiscrete*>(currState[i])->getState() ] -= 1;
            counts[s] += 1;

            // update history
            t -= dt;
            currState[i] = evt;
        }
    }

    bh->setHistory(history);

    for (size_t i = 0; i < this->num_sites; i++)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(currState[i])->getState();
        currState[i] = new CharacterEventDiscrete(i, s, 1.0);
    }

    bh->setChildCharacters(currState);
    
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::simulate(const TopologyNode& node, BranchHistory* bh, std::vector< DiscreteTaxonData< charType > >& taxa)
{


    //    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the sequence of this node
    size_t nodeIndex = node.getIndex();

    // get rate map for branch leading to node
    //    const RateGeneratorSequence& rm = homogeneousRateGeneratorSequence->getValue();

    // if root, set tail state
    if ( node.isRoot() == true )
    {
        std::vector<CharacterEvent*> childState;
        std::vector<CharacterEvent*> parentState;
        std::vector<double> rfs = getRootFrequencies();
        for (size_t i = 0; i < this->num_sites; i++)
        {
            unsigned s = 0;
            double u = GLOBAL_RNG->uniform01();
            for (; s < this->num_states; ++s)
            {
                u -= rfs[s];
                if (u <= 0.0)
                {
                    break;
                }

            }

            parentState.push_back(new CharacterEventDiscrete(i, s, 0.0));
            childState.push_back(new CharacterEventDiscrete(i, s, 1.0));
        }

        bh->setParentCharacters(parentState);
        bh->setChildCharacters(childState);

    }
    else
    {
        // simulate anagenic changes
        simulateHistory(node, bh);
    }

    const std::vector<CharacterEvent*>& childState = bh->getChildCharacters();
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        // create the character
        charType c = charType( this->template_state );
        c.setStateByIndex( static_cast<CharacterEventDiscrete*>(childState[i])->getState() );
        taxa[nodeIndex].addCharacter( c );
    }

    if ( node.isTip() == true )
    {
        taxa[nodeIndex].setTaxon( node.getTaxon() );
    }
    else
    {
        const std::vector<TopologyNode*>& children = node.getChildren();
        for (size_t i = 0; i < children.size(); ++i)
        {
            this->histories[ children[i]->getIndex() ] = BranchHistoryDiscrete(this->num_sites, this->num_states, children[i]->getIndex() );
            std::vector<CharacterEvent*> childParentCharacters = this->histories[children[i]->getIndex()].getParentCharacters();
            for (size_t j = 0; j < this->num_sites; ++j)
            {
                static_cast<CharacterEventDiscrete*>(childParentCharacters[j])->setState( static_cast<CharacterEventDiscrete*>(childState[j])->getState() );
            }
        }

        for (size_t i = 0; i < children.size(); ++i)
        {
            BranchHistory* bh_ch = &(this->histories[ children[i]->getIndex() ]);
            simulate( *children[i], bh_ch, taxa );
        }

    }

}


/** Swap a parameter of the distribution */
template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{

    if (oldP == homogeneousRateGenerator)
    {
        homogeneousRateGenerator = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    }
    else if (oldP == heterogeneousRateGenerator)
    {
        heterogeneousRateGenerator = static_cast<const TypedDagNode< RbVector< RateGenerator > >* >( newP );
    }
    else if (oldP == rootFrequencies)
    {
        rootFrequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else
    {
        TreeHistoryCtmc<charType>::swapParameterInternal(oldP,newP);
    }

}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmcSiteIID<charType>::touchSpecialization( const DagNode* affecter, bool touchAll )
{

    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == rootFrequencies )
    {
        const TopologyNode &root = this->tau->getValue().getRoot();
        this->flagNodeDirty( root );
    }
    else
    {
        TreeHistoryCtmc<charType>::touchSpecialization( affecter, touchAll );
    }

}

#endif
