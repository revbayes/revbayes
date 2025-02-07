#ifndef GeneralTreeHistoryCtmc_H
#define GeneralTreeHistoryCtmc_H

#include "TreeHistoryCtmc.h"
#include "AbstractCladogenicStateFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "ConstantNode.h"
#include "DeterministicNode.h"
#include "DiscreteCharacterState.h"
#include "HomologousDiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RateGeneratorSequence.h"
#include "RateGeneratorSequenceUsingMatrix.h"
#include "RbConstants.h"
#include "RbVector.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"

#include <cmath>
#include <cstring>
#include <set>
#include <vector>

namespace RevBayesCore {
    
    
    template<class charType>
    class GeneralTreeHistoryCtmc : public TreeHistoryCtmc<charType> {
        
    public:
        GeneralTreeHistoryCtmc(const TypedDagNode< Tree > *t, size_t nChars, size_t nSites, bool useAmbigChar=false, bool internal=false);
        GeneralTreeHistoryCtmc(const GeneralTreeHistoryCtmc &n);                                                                         //!< Copy constructor
        virtual                                            ~GeneralTreeHistoryCtmc(void);                                                //!< Virtual destructor
        
        // public member functions
        
        GeneralTreeHistoryCtmc*                             clone(void) const;                                                           //!< Create an independent clone
        virtual void                                        initializeTipValues(void);
        virtual bool                                        drawInitValue(void);
        virtual std::vector<double>                         getRootFrequencies(void) const;
        virtual std::vector<std::string>                    getCladogeneticEvents(void) const;
        virtual std::string                                 getCladogeneticEvent( size_t index ) const;
        const CladogeneticProbabilityMatrix&                getCladogeneticProbabilityMatrix( void ) const;
        virtual void                                        redrawValue(void);
        virtual void                                        reInitialized(void);
        virtual void                                        simulate(void);
        
        virtual bool                                        samplePathStart(const TopologyNode& node);
        virtual bool                                        samplePathEnd(const TopologyNode& node);
        virtual bool                                        samplePathHistory(const TopologyNode& node);
        
        void                                                setCladogeneticProbabilityMatrix(const DeterministicNode< CladogeneticProbabilityMatrix>* cp );
        void                                                setCladogeneticEvent(std::string s, size_t idx);
        
        void                                                setRateGenerator(const TypedDagNode< RateGeneratorSequence > *rm);
        void                                                setRateGenerator(const TypedDagNode< RbVector< RateGeneratorSequence > > *rm);
        void                                                setRootFrequencies(const TypedDagNode< Simplex > *f);
        void                                                setSiteRates(const TypedDagNode< RbVector< double > > *r);
        
    protected:
        
        // Parameter management functions
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter
        
        virtual double                                      computeRootLikelihood(const TopologyNode &n);
        virtual double                                      computeInternalNodeLikelihood(const TopologyNode &n);
        virtual double                                      computeTipLikelihood(const TopologyNode &node);
        virtual void                                        touchSpecialization(const DagNode *toucher, bool touchAll);
        
    private:
        
        virtual void                                        simulate(const TopologyNode& node, BranchHistory* bh, std::vector< DiscreteTaxonData< charType > >& taxa);
        void                                                simulateHistory(const TopologyNode& node, BranchHistory* bh);
        std::vector<size_t>                                 computeCounts(const std::vector<CharacterEvent*>& s);
        std::vector<std::set< size_t> >                     computeSitesWithStates(const std::vector<CharacterEvent*>& s);
        
        // members
        const TypedDagNode< RateGeneratorSequence >*                homogeneousRateGenerator;
        const TypedDagNode< RbVector< RateGeneratorSequence > >*    heterogeneousRateGenerator;
        const DeterministicNode< CladogeneticProbabilityMatrix >*   homogeneousCladogeneticProbabilityMatrix;
        const TypedDagNode< Simplex >*                              rootFrequencies;
        const TypedDagNode< RbVector< double > >*                   siteRates;
        
        // flags specifying which model variants we use
        bool                                                branchHeterogeneousSubstitutionMatrices;
        bool                                                rateVariationAcrossSites;
        bool                                                useCladogeneticEvents;
        
        // cladogenetic histories
        std::vector<std::string>                            cladogeneticEvents;
        
    };
    
}


#include "BranchHistoryDiscrete.h"
#include "CharacterEventDiscrete.h"
#include "PathRejectionSampleProposal.h"
#include "RateMatrix_JC.h"


template<class charType>
RevBayesCore::GeneralTreeHistoryCtmc<charType>::GeneralTreeHistoryCtmc(const TypedDagNode<Tree> *tau, size_t nChars, size_t nSites, bool useAmbigChar, bool internal) : TreeHistoryCtmc<charType>( tau, nChars, nSites, useAmbigChar, internal )
{
    
    // initialize with default parameters
    homogeneousRateGenerator    = new ConstantNode<RateGeneratorSequence>("rateGenerator", new RateGeneratorSequenceUsingMatrix( nChars, nSites ) );
    heterogeneousRateGenerator  = NULL;
    rootFrequencies             = new ConstantNode<Simplex>("rootFrequencies", new Simplex(nChars, 1.0/nChars));
    siteRates                   = NULL;
    homogeneousCladogeneticProbabilityMatrix = NULL;
    cladogeneticEvents          = std::vector<std::string>(tau->getValue().getNumberOfNodes(), "");
    
    // flags specifying which model variants we use
    branchHeterogeneousSubstitutionMatrices     = false;
    rateVariationAcrossSites                    = false;
    useCladogeneticEvents                       = false;
    
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do

    this->addParameter( tau );
    this->addParameter( homogeneousRateGenerator );
    this->addParameter( heterogeneousRateGenerator );
    this->addParameter( rootFrequencies );
    this->addParameter( siteRates );
    this->addParameter( homogeneousCladogeneticProbabilityMatrix );
    
    this->useDirtyNodes = true;
    
}


template<class charType>
RevBayesCore::GeneralTreeHistoryCtmc<charType>::GeneralTreeHistoryCtmc(const GeneralTreeHistoryCtmc &d) : TreeHistoryCtmc<charType>( d )
{
    // initialize with default parameters
    homogeneousRateGenerator    = d.homogeneousRateGenerator;
    heterogeneousRateGenerator  = d.heterogeneousRateGenerator;
    homogeneousCladogeneticProbabilityMatrix = d.homogeneousCladogeneticProbabilityMatrix;
    rootFrequencies             = d.rootFrequencies;
    siteRates                   = d.siteRates;
    cladogeneticEvents          = d.cladogeneticEvents;
    
    // flags specifying which model variants we use
    branchHeterogeneousSubstitutionMatrices     = d.branchHeterogeneousSubstitutionMatrices;
    rateVariationAcrossSites                    = d.rateVariationAcrossSites;
    useCladogeneticEvents                       = d.useCladogeneticEvents;
    
}


template<class charType>
RevBayesCore::GeneralTreeHistoryCtmc<charType>::~GeneralTreeHistoryCtmc( void )
{
    ; // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


template<class charType>
RevBayesCore::GeneralTreeHistoryCtmc<charType>* RevBayesCore::GeneralTreeHistoryCtmc<charType>::clone( void ) const
{
    
    return new GeneralTreeHistoryCtmc<charType>( *this );
}

template<class charType>
std::vector<size_t> RevBayesCore::GeneralTreeHistoryCtmc<charType>::computeCounts(const std::vector<CharacterEvent*>& s)
{
    std::vector<size_t> counts(this->num_states, 0);
    
    for (size_t i = 0; i < s.size(); i++)
    {
        counts[ static_cast<CharacterEventDiscrete*>(s[i])->getState() ] += 1;
    }
    
    return counts;
}

template<class charType>
std::vector<std::set<size_t> > RevBayesCore::GeneralTreeHistoryCtmc<charType>::computeSitesWithStates(const std::vector<CharacterEvent*>& s)
{
    std::vector<std::set<size_t> > sites_with_states(this->num_states, std::set<size_t>() );
    
    for (size_t i = 0; i < s.size(); i++)
    {
        sites_with_states[ static_cast<CharacterEventDiscrete*>(s[i])->getState() ].insert( static_cast<CharacterEventDiscrete*>(s[i])->getSiteIndex() );
        //counts[  ] += 1;
    }
    
    return sites_with_states;
}

template<class charType>
double RevBayesCore::GeneralTreeHistoryCtmc<charType>::computeRootLikelihood(const TopologyNode &n)
{
    
    double lnP = 0.0;
    
    BranchHistory* bh = this->histories[n.getIndex()];
    const std::vector<CharacterEvent*> rootState = bh->getParentCharacters();
    
    // get counts per state
    std::vector<int> counts(this->num_states, 0);
    for (size_t i = 0; i < rootState.size(); ++i)
    {
        ++counts[ static_cast<CharacterEventDiscrete*>(rootState[i])->getState() ];
    }
    
    // get log prob
    std::vector<double> rf = getRootFrequencies();
    for (size_t i = 0; i < counts.size(); i++)
    {
        // skip unused states
        if (counts[i] != 0) {
            lnP += counts[i] * log( rf[i] );
        }
    }

    return lnP;
}

template<class charType>
double RevBayesCore::GeneralTreeHistoryCtmc<charType>::computeInternalNodeLikelihood(const TopologyNode &node)
{
    

    size_t node_index = node.getIndex();
    double branch_rate = this->getBranchRate(node_index);
    const RateGeneratorSequence& rm = homogeneousRateGenerator->getValue();
    
    // get the branch history
    BranchHistory* bh = this->histories[node_index];
    std::vector<CharacterEvent*> curr_state = bh->getParentCharacters();
    std::vector<CharacterEvent*> end_state  = bh->getChildCharacters();

    // check that node ages are consistent with character event ages
    if ( bh->areEventTimesValid(this->tau->getValue(), node) == false )
    {
        return RbConstants::Double::neginf;
    }
    
    // check parent and child states to make sure they match with the
    // ancestral and descendant branches; otherwise, return -Inf
    if (!useCladogeneticEvents) {
        for (size_t i = 0; i < node.getNumberOfChildren(); ++i)
        {
            const TopologyNode &child = node.getChild(i);
            size_t child_index = child.getIndex();
            BranchHistory* child_bh = this->histories[child_index];
            std::vector<CharacterEvent*> child_state = child_bh->getParentCharacters();
            for (size_t j = 0; j < this->num_sites; ++j)
            {
                if (static_cast<CharacterEventDiscrete*>(end_state[j])->getState() != static_cast<CharacterEventDiscrete*>(child_state[j])->getState() )
                {
                    return RbConstants::Double::neginf;
                }
            }
        }
    }
    
    // we need the counts for faster computation
    std::vector<std::set<size_t> > sites_with_states = computeSitesWithStates(curr_state);
    
    // get branch history set and iterator
    const std::multiset<CharacterEvent*,CharacterEventCompare>& history = bh->getHistory();
    std::multiset<CharacterEvent*,CharacterEventCompare>::reverse_iterator it_h;
    
    // stepwise events
    double lnL = 0.0;
    double current_age = 0.0;
    if (node.isRoot()) {
        current_age = node.getAge() + this->getRootBranchLength();
    }
    else {
        current_age = node.getParent().getAge();
    }
    
    double end_age = node.getAge();
    double event_age;
    double sr = rm.getSumOfRates(curr_state, sites_with_states, current_age, branch_rate);

    for (it_h = history.rbegin(); it_h != history.rend(); ++it_h)
    {
        CharacterEventDiscrete* char_event = static_cast<CharacterEventDiscrete*>(*it_h);
        
        // next event time
        double idx = char_event->getSiteIndex();
        event_age = char_event->getAge();
        size_t s = char_event->getState();
        
        // lnL for stepwise events for p(x->y)
        double tr = rm.getRate(curr_state, char_event, sites_with_states, current_age, branch_rate);
        lnL += log(tr) - sr * (current_age - event_age);

        // path history is -Inf if any step is -Inf
        if (lnL == RbConstants::Double::neginf) return lnL;
        
        // update counts
        sites_with_states[ static_cast<CharacterEventDiscrete*>(curr_state[idx])->getState() ].erase( idx );
        sites_with_states[ s ].insert( idx );
        
        // update time and state
        curr_state[idx] = char_event;
        current_age = event_age;
        
        // update sum of rates
        sr = rm.getSumOfRates(curr_state, sites_with_states, current_age, branch_rate);

    }
    
    // lnL that nothing else happens
    lnL -= sr * (current_age - end_age);

    // account for the cladogenetic transition probabilites, if applicable
    if (useCladogeneticEvents && !node.isTip())
    {
        const AbstractCladogenicStateFunction* cf = dynamic_cast<const AbstractCladogenicStateFunction* >( &homogeneousCladogeneticProbabilityMatrix->getFunction() );
        size_t left_index = node.getChild(0).getIndex();
        size_t right_index = node.getChild(1).getIndex();
        double lnL_clado = cf->computeDataAugmentedCladogeneticLnProbability( this->histories, node_index, left_index, right_index );
        lnL += lnL_clado;
        
    }
    
    return lnL;
}



template<class charType>
double RevBayesCore::GeneralTreeHistoryCtmc<charType>::computeTipLikelihood(const TopologyNode &node)
{
    double lnL = 0.0;
    return lnL;
}


//template<class charType>
//void RevBayesCore::GeneralTreeHistoryCtmc<charType>::initializeTipValues( void )
//{
//    //    if (this->dag_node->isClamped())
//    {
//        std::vector<TopologyNode*> nodes = RevBayesCore::AbstractTreeHistoryCtmc<charType>::tau->getValue().getNodes();
//        for (size_t i = 0; i < nodes.size(); i++)
//        {
//            TopologyNode* node = nodes[i];
//            if (node->isTip())
//            {
//                DiscreteTaxonData<charType>& d = static_cast< DiscreteTaxonData<charType>& >( this->value->getTaxonData( node->getName() ) );
//                
//                std::vector<CharacterEvent*> tipState;
//                for (size_t j = 0; j < d.getNumberOfCharacters(); j++)
//                {
//                    unsigned s = 0;
////                    if (!this->using_ambiguous_characters)
//                    s = (unsigned)d[j].getStateIndex();
//                    
//                    CharacterEvent* evt = new CharacterEvent(j, s, 1.0);
//                    tipState.push_back( evt );
//                }
//                
//                this->histories[node->getIndex()]->setChildCharacters(tipState);
//                tipState.clear();
//            }
//        }
//        this->tipsInitialized = true;
//    }
//}

template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmc<charType>::drawInitValue( void )
{
    
    // convert the tip values of the data matrix into branch history objects
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
    
    // apply cladogenetic events, if applicable
    if (useCladogeneticEvents)
    {
        for (size_t i = 0; i < nodes.size(); i++)
        {
            TopologyNode* nd = nodes[i];
            if (nd->isTip() == false) {
                const AbstractCladogenicStateFunction* cf = dynamic_cast<const AbstractCladogenicStateFunction* >( &homogeneousCladogeneticProbabilityMatrix->getFunction() );
                size_t node_index = nd->getIndex();
                size_t left_index = nd->getChild(0).getIndex();
                size_t right_index = nd->getChild(1).getIndex();
                std::string clado_state = cf->simulateDataAugmentedCladogeneticState( this->histories, node_index, left_index, right_index );
                cladogeneticEvents[ node_index ] = clado_state;
            }
        }
    }
    
    // sample paths
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        TopologyNode* nd = nodes[i];
        
        int samplePathHistoryCount = 0;
        
        bool valid_sample = false;
        bool valid_likelihood = false;
        
        while ( !(valid_sample && valid_likelihood) && samplePathHistoryCount < 100)
        {
            samplePathHistoryCount++;
            valid_sample = samplePathHistory(*nd);
            double branch_lnL = computeInternalNodeLikelihood(*nd);
            
            valid_likelihood = branch_lnL != RbConstants::Double::neginf;
            
        }
    }
//    std::cout << "----\n";
//    std::cout << "Init\n";
//    std::cout << "----\n";
//    for (size_t i = 0; i < nodes.size(); i++)
//    {
//        TopologyNode* nd = nodes[i];
//        this->histories[ nd->getIndex() ]->print( nd );
//    }
//    
//    std::cout << "----\n";
    
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
std::vector<double> RevBayesCore::GeneralTreeHistoryCtmc<charType>::getRootFrequencies( void ) const
{
    
    if ( branchHeterogeneousSubstitutionMatrices == true || rootFrequencies != NULL )
    {
        return rootFrequencies->getValue();
    }
    else
    {
        const RateGeneratorSequence *rm = dynamic_cast<const RateGeneratorSequence *>( &homogeneousRateGenerator->getValue() );
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
std::vector<std::string> RevBayesCore::GeneralTreeHistoryCtmc<charType>::getCladogeneticEvents( void ) const
{
    return cladogeneticEvents;
}

template<class charType>
std::string RevBayesCore::GeneralTreeHistoryCtmc<charType>::getCladogeneticEvent( size_t index ) const
{
    return cladogeneticEvents[ index ];
}

template<class charType>
const RevBayesCore::CladogeneticProbabilityMatrix& RevBayesCore::GeneralTreeHistoryCtmc<charType>::getCladogeneticProbabilityMatrix( void ) const
{
    return homogeneousCladogeneticProbabilityMatrix->getValue();
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::initializeTipValues( void )
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
                    
                    s = (unsigned) state.getStateIndex();
                    
                    if ( state.isGapState() == true )
                    {
                        s = 0;
                    }
                    
                    CharacterEvent* evt = new CharacterEventDiscrete(j, s, 1.0);
                    tipState.push_back( evt );
                }
                
                this->histories[node->getIndex()]->setChildCharacters(tipState);
                //                tipState.clear();
            }
        }
        
        this->tipsInitialized = true;
    }
    
}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::redrawValue( void )
{
    
    simulate();
    
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::reInitialized( void ) {
    
    // redraw the internal history states
    drawInitValue();
    
}

template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmc<charType>::samplePathEnd(const TopologyNode& node)
{
    if ( node.isTip() == true )
    {
        return true; // do nothing
    }
    else
    {
        TransitionProbabilityMatrix leftTpMatrix(this->num_states);
        TransitionProbabilityMatrix rightTpMatrix(this->num_states);
        
        const RateGeneratorSequence& rm = homogeneousRateGenerator->getValue();
        
        double begin_age = node.getAge();
        
        size_t left_index  = node.getChild(0).getIndex();
        size_t right_index = node.getChild(1).getIndex();
        
        double left_branch_rate  = this->getBranchRate( left_index );
        double right_branch_rate = this->getBranchRate( right_index );
        
        rm.calculateTransitionProbabilities(begin_age, node.getChild(0).getAge(), left_branch_rate, leftTpMatrix);
        rm.calculateTransitionProbabilities(begin_age, node.getChild(1).getAge(), right_branch_rate, rightTpMatrix);
        
        // for sampling probs
        const std::vector<CharacterEvent*>& leftChildState  = this->histories[left_index]->getChildCharacters();
        const std::vector<CharacterEvent*>& rightChildState = this->histories[right_index]->getChildCharacters();
        
        // to update
        std::vector<CharacterEvent*> nodeChildState = this->histories[node.getIndex()]->getChildCharacters();
        std::vector<CharacterEvent*> nodeParentState = this->histories[node.getIndex()]->getParentCharacters();
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
            static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
        }
        
       
    }
    
    if ( node.isRoot() )
    {
        TransitionProbabilityMatrix nodeTpMatrix(this->num_states);
        
        const RateGeneratorSequence& rm = homogeneousRateGenerator->getValue();
        double begin_age = node.getAge() + this->getRootBranchLength();
        size_t node_index  = node.getIndex();
        double node_branch_rate = this->getBranchRate( node_index );
        rm.calculateTransitionProbabilities(begin_age, node.getAge(), node_branch_rate, nodeTpMatrix);
        
        // for sampling probs
        const std::vector<CharacterEvent*>& nodeChildState  = this->histories[node_index]->getChildCharacters();
        
        // to update
        std::vector<CharacterEvent*> nodeParentState = this->histories[node_index]->getParentCharacters();
        
        for (size_t site_index=0; site_index<this->num_sites; ++site_index)
        {
            
            size_t desS1 = static_cast<CharacterEventDiscrete*>(nodeChildState[site_index])->getState();
            
            std::vector<double> state_probs(this->num_states, 0.0);
            double prob_sum = 0.0;
            for (size_t i = 0; i < this->num_states; ++i)
            {
                state_probs[i] = nodeTpMatrix[i][desS1];
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
            
            static_cast<CharacterEventDiscrete*>(nodeParentState[site_index])->setState(s);
        }
    }
    
    return true;
}


template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmc<charType>::samplePathHistory(const TopologyNode& node)
{
    
//    if ( node.isRoot() == true )
//    {
//        return true;
//    }
//    
    PathRejectionSampleProposal<charType> p( this->getStochasticNode() );
    p.setRateGenerator( homogeneousRateGenerator );
    
    p.assignNode(const_cast<TopologyNode*>(&node));
    
    p.prepareProposal();
    p.doProposal();
    p.cleanProposal();
    
    return true;
}


template<class charType>
bool RevBayesCore::GeneralTreeHistoryCtmc<charType>::samplePathStart(const TopologyNode& node)
{
    
    // ignore tips
    if ( node.isTip() == true )
    {
        return true;
    }
    
    // update daughter lineages' parent states with the current node's child state
    const std::vector<CharacterEvent*>& nodeChildState = this->histories[ node.getIndex() ]->getChildCharacters();
    const std::vector<TopologyNode*>& children = node.getChildren();
    for (size_t i = 0; i < children.size(); ++i)
    {
        std::vector<CharacterEvent*> childParentState(nodeChildState.size());
        for (size_t j = 0; j < childParentState.size(); ++j)
        {
            childParentState[j] = new CharacterEventDiscrete(j, static_cast<CharacterEventDiscrete*>(nodeChildState[j])->getState(), 0.0);
        }
        this->histories[ children[i]->getIndex() ]->setParentCharacters( childParentState );
    }
    
    return true;
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::setCladogeneticProbabilityMatrix(const DeterministicNode< CladogeneticProbabilityMatrix > *cp) {
    
    // remove the old parameter first
    if ( homogeneousCladogeneticProbabilityMatrix != NULL )
    {
        this->removeParameter( homogeneousCladogeneticProbabilityMatrix );
        homogeneousCladogeneticProbabilityMatrix = NULL;
    }
    
    // set the value
    useCladogeneticEvents = true;
    homogeneousCladogeneticProbabilityMatrix = cp;
    
    // add the new parameter
    this->addParameter( homogeneousCladogeneticProbabilityMatrix );
    
    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
    
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::setCladogeneticEvent(std::string s, size_t idx)
{
    cladogeneticEvents[ idx ] = s;
}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::setRateGenerator(const TypedDagNode< RateGeneratorSequence > *rm) {
    
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
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::setRateGenerator(const TypedDagNode< RbVector< RateGeneratorSequence > > *rm) {
    
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
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::setRootFrequencies(const TypedDagNode< Simplex > *f)
{
    
    // remove the old parameter if f is a valid replacement
    if ( rootFrequencies != NULL && f != NULL )
    {
        this->removeParameter( rootFrequencies );
        rootFrequencies = NULL;
    }
    
    if ( f != NULL )
    {
        rootFrequencies = f;
        
        // add the parameter
        this->addParameter( rootFrequencies );
    }

    // redraw the current value
    if ( this->dag_node != NULL && !this->dag_node->isClamped() )
    {
        this->redrawValue();
    }
}


template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::setSiteRates(const TypedDagNode< RbVector< double > > *r) {
    
    // remove the old parameter first
    if ( siteRates != NULL )
    {
        this->removeParameter( siteRates );
        siteRates = NULL;
    }
    if ( r != NULL )
    {
        // set the value
        rateVariationAcrossSites = true;
        siteRates = r;
        this->num_site_rates = r->getValue().size();
        // MJL: not needed?
        // this->resizeLikelihoodVectors();
    }
    else
    {
        // set the value
        rateVariationAcrossSites = false;
        siteRates = NULL;
        this->num_site_rates = 1;
        // MJL: not needed?
        // this->resizeLikelihoodVectors();
        
    }
    
    // add the new parameter
    this->addParameter( siteRates );
    
    // redraw the current value
    if ( this->dag_node != NULL && !this->dag_node->isClamped() )
    {
        this->redrawValue();
    }
    
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::simulate(void)
{
    TreeHistoryCtmc<charType>::simulate();
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::simulate(const TopologyNode& node, BranchHistory* bh, std::vector< DiscreteTaxonData< charType > >& taxa)
{

    // get the sequence of this node
    size_t node_index = node.getIndex();
    double end_age = node.getAge();
    double start_age = 0.0;
    if (node.isRoot()) {
        start_age = node.getAge() + this->getRootBranchLength();
    }
    else {
        start_age = node.getParent().getAge();
    }
    

    // if root, set tail state
    if ( node.isRoot() == true )
    {
        std::vector<CharacterEvent*> childState;
        std::vector<CharacterEvent*> parentState;
        const std::vector<double>& rfs = rootFrequencies->getValue();
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
            
            parentState.push_back(new CharacterEventDiscrete(i, s, start_age ));
            childState.push_back(new CharacterEventDiscrete(i, s, end_age ));
        }
        
        // set parent states on the root's tail
        bh->setParentCharacters(parentState);
        
        simulateHistory(node, bh);
        
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
        taxa[node_index].addCharacter( c );
    }
    
    // bh->print( &node );
    
    if ( node.isTip() == true )
    {
        taxa[node_index].setTaxon( node.getTaxon() );
    }
    else
    {
        const std::vector<TopologyNode*>& children = node.getChildren();

        // simulate a cladogenetic state
        if (useCladogeneticEvents) {
            const AbstractCladogenicStateFunction* cf = dynamic_cast<const AbstractCladogenicStateFunction* >( &homogeneousCladogeneticProbabilityMatrix->getFunction() );
            size_t left_index = node.getChild(0).getIndex();
            size_t right_index = node.getChild(1).getIndex();
            cf->simulateDataAugmentedCladogeneticState( this->histories, node_index, left_index, right_index );
        }
        else {
            for (size_t i = 0; i < children.size(); ++i)
            {
                this->histories[ children[i]->getIndex() ] = new BranchHistoryDiscrete(this->num_sites, this->num_states, children[i]->getIndex() );
                std::vector<CharacterEvent*> childParentCharacters = this->histories[children[i]->getIndex()]->getParentCharacters();
                for (size_t j = 0; j < this->num_sites; ++j)
                {
                    static_cast<CharacterEventDiscrete*>(childParentCharacters[j])->setState( static_cast<CharacterEventDiscrete*>(childState[j])->getState() );
                }
            }
        }
        
        for (size_t i = 0; i < children.size(); ++i)
        {
            BranchHistory* bh_ch = this->histories[ children[i]->getIndex() ];
            simulate( *children[i], bh_ch, taxa );
        }
        
    }
    
    
}




template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::simulateHistory(const TopologyNode& node, BranchHistory* bh)
{
    
    size_t branch_index = node.getIndex();
    double branch_rate = this->getBranchRate( branch_index );
    const RateGeneratorSequence& rm = homogeneousRateGenerator->getValue();
    
    // get the start and end times
    
    double end_age = node.getAge();
    double start_age = 0.0;
    
    if (node.isRoot()) {
        start_age = node.getAge() + this->getRootBranchLength();
    }
    else {
        start_age = node.getParent().getAge();
    }
    double branch_length = start_age - end_age;
    
    // get parent BranchHistory state
    std::vector<CharacterEvent*> currState = bh->getParentCharacters();
    std::vector<size_t> counts = computeCounts(currState);
    std::set<CharacterEvent*,CharacterEventCompare> history;
    
    // get start sum of rates
    // NOTE: cannot call sum of rates using counts because it may be non-iid evolution
    double sr = rm.getSumOfRates(currState, start_age, branch_rate);
    
    // simulate path
    double curr_age = start_age;
    double dt = 0.0;
    while (curr_age - dt > end_age)
    {
        
        // sample next event time
        dt = RbStatistics::Exponential::rv(sr, *GLOBAL_RNG);
        if (curr_age - dt > end_age)
        {
    
            // next event type
            CharacterEventDiscrete* evt = new CharacterEventDiscrete(0, 0, curr_age - dt);
            double u = GLOBAL_RNG->uniform01() * sr;
            
            bool found = false;
            size_t i;
            size_t s = 0;
            for (i = 0; !found && i < this->num_sites; i++)
            {
                
                evt->setSiteIndex(i);
                size_t curr_state = static_cast<CharacterEventDiscrete*>(currState[i])->getState();
                for (s = 0; !found && s < this->num_states; ++s)
                {
                    // disregard virtual events (self-transitions)
                    if (s != curr_state)
                    {
                        evt->setState(s);
                        double r = rm.getRate(currState, evt, curr_age, branch_rate);
                        
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
            curr_age -= dt;
            currState[i] = evt;
            
            // update sum of rates
            sr = rm.getSumOfRates(currState, curr_age, branch_rate);
        }
    }
    
    bh->setHistory(history);
    
    
    for (size_t i = 0; i < this->num_sites; i++)
    {
        size_t s = static_cast<CharacterEventDiscrete*>(currState[i])->getState();
        currState[i] = new CharacterEventDiscrete(i, s, end_age);
    }
    
    bh->setChildCharacters(currState);
    
}


/** Swap a parameter of the distribution */
template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    
    if (oldP == homogeneousRateGenerator)
    {
        homogeneousRateGenerator = static_cast<const TypedDagNode< RateGeneratorSequence >* >( newP );
    }
    else if (oldP == heterogeneousRateGenerator)
    {
        heterogeneousRateGenerator = static_cast<const TypedDagNode< RbVector< RateGeneratorSequence > >* >( newP );
    }
    else if (oldP == homogeneousCladogeneticProbabilityMatrix)
    {
        homogeneousCladogeneticProbabilityMatrix = static_cast<const DeterministicNode< CladogeneticProbabilityMatrix >* >( newP );
    }
    else if (oldP == rootFrequencies)
    {
        rootFrequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == siteRates)
    {
        siteRates = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else
    {
        TreeHistoryCtmc<charType>::swapParameterInternal(oldP,newP);
    }
    
}

template<class charType>
void RevBayesCore::GeneralTreeHistoryCtmc<charType>::touchSpecialization( const DagNode* affecter, bool touchAll )
{
    
    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == rootFrequencies || affecter == rootFrequencies )
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
