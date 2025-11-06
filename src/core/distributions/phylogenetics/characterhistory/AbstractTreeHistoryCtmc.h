#ifndef AbstractTreeHistoryCtmc_H
#define AbstractTreeHistoryCtmc_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "BranchHistory.h"
#include "ContinuousCharacterData.h"
#include "DiscreteTaxonData.h"
#include "HomologousDiscreteCharacterData.h"
#include "DiscreteCharacterState.h"
#include "DnaState.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMatrix.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

#include <cmath>

namespace RevBayesCore {

    template<class charType>
    class AbstractTreeHistoryCtmc : public TypedDistribution< AbstractHomologousDiscreteCharacterData >, public TreeChangeEventListener {

    public:
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractTreeHistoryCtmc(const TypedDagNode<Tree> *t, size_t nChars, size_t nSites, bool useAmbigChar=false);
        AbstractTreeHistoryCtmc(const AbstractTreeHistoryCtmc &n);                                                                              //!< Copy constructor
        virtual                                                            ~AbstractTreeHistoryCtmc(void);                                      //!< Virtual destructor

        // public member functions
        // pure virtual
        virtual AbstractTreeHistoryCtmc*                                    clone(void) const = 0;                                              //!< Create an independent clone
        virtual void                                                        redrawValue(void) = 0;
        virtual bool                                                        drawInitValue(void) = 0;
        virtual void                                                        initializeTipValues(void) = 0;
        virtual bool                                                        samplePathStart(const TopologyNode& node, const std::set<size_t>& indexSet) = 0;
        virtual bool                                                        samplePathEnd(const TopologyNode& node, const std::set<size_t>& indexSet) = 0;
        virtual bool                                                        samplePathHistory(const TopologyNode& node, const std::set<size_t>& indexSet) = 0;

        // non-virtual
        double                                                              computeLnProbability(void);
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                         //!< The tree has changed and we want to know which part.
        BranchHistory&                                                      getHistory(const TopologyNode& nd);
        const BranchHistory&                                                getHistory(const TopologyNode& nd) const;
        std::vector<BranchHistory*>                                         getHistories(void);
        const std::vector<BranchHistory*>&                                  getHistories(void) const;
        void                                                                setHistory(const BranchHistory& bh, const TopologyNode& nd);
        void                                                                setHistories(const std::vector<BranchHistory*>& bh);
        void                                                                setValue(AbstractHomologousDiscreteCharacterData *v, bool f=false);           //!< Set the current value, e.g. attach an observation (clamp)
        void                                                                setTipProbs(const HomologousCharacterData* tp);

        virtual const std::vector<double>&                                  getTipProbs(const TopologyNode& nd);
        virtual const std::vector<std::vector<double> >&                    getTipProbs(void);

        virtual void                                                        simulate(void);


    protected:

        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

        // helper method for this and derived classes
        void                                                                flagNodeDirty(const TopologyNode& n);
        double                                                              computeBranchRate(size_t index);

        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);

        // pure virtual methods
        virtual double                                                      computeRootLikelihood(const TopologyNode &nd) = 0;
        virtual double                                                      computeInternalNodeLikelihood(const TopologyNode &nd) = 0;
        virtual double                                                      computeTipLikelihood(const TopologyNode &nd) = 0;
        virtual const std::vector<double>&                                  getRootFrequencies(void) = 0;

        // members
        double                                                              lnProb;
        const size_t                                                        num_chars;
        size_t                                                              num_sites;
        size_t                                                              num_site_rates;
        const TypedDagNode<Tree>*                                           tau;
        
        // the likelihoods
        std::vector<size_t>                                                 activeLikelihood;
        std::vector<std::vector<double> >                                   historyLikelihoods;

        // the data
        std::vector<std::vector<std::uint64_t> >                            char_matrix;
        std::vector<std::vector<bool> >                                     gap_matrix;
        std::vector<BranchHistory*>                                         histories;
        std::vector<std::vector<double> >                                   tipProbs;

        // convenience variables available for derived classes too
        std::vector<bool>                                                   changed_nodes;
        std::vector<bool>                                                   dirty_nodes;
        
        // flags
        bool                                                                using_ambiguous_characters;
        bool                                                                treatUnknownAsGap;
        bool                                                                treatAmbiguousAsGaps;
        bool                                                                tipsInitialized;

    private:
        // private methods
        void                                                                fillLikelihoodVector(const TopologyNode &n);
        void                                                                initializeHistoriesVector(void);
        virtual void                                                        simulate(const TopologyNode& node, BranchHistory* bh, std::vector<DiscreteTaxonData< charType > >& taxa) = 0;
    };
}

template<class charType>
RevBayesCore::AbstractTreeHistoryCtmc<charType>::AbstractTreeHistoryCtmc(const TypedDagNode<Tree> *t, size_t nChars, size_t nSites, bool useAmbigChar) : TypedDistribution< AbstractHomologousDiscreteCharacterData >(  new HomologousDiscreteCharacterData<charType>() ),
num_chars( nChars ),
num_sites( nSites ),
num_site_rates( 1 ),
tau( t ),
activeLikelihood( std::vector<size_t>(tau->getValue().getNumberOfNodes(), 0) ),
historyLikelihoods(),
char_matrix(),
gap_matrix(),
histories(),
changed_nodes( std::vector<bool>(tau->getValue().getNumberOfNodes(),false) ),
dirty_nodes( std::vector<bool>(tau->getValue().getNumberOfNodes(), true) ),
using_ambiguous_characters( useAmbigChar ),
treatUnknownAsGap( true ),
treatAmbiguousAsGaps( true ),
tipsInitialized( false )
{
    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( tau );

    // We don't want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );

    // initialize histories
    initializeHistoriesVector();

}


template<class charType>
RevBayesCore::AbstractTreeHistoryCtmc<charType>::AbstractTreeHistoryCtmc(const AbstractTreeHistoryCtmc &n) : TypedDistribution< AbstractHomologousDiscreteCharacterData >( n ),
num_chars( n.num_chars ),
num_sites( n.num_sites ),
num_site_rates( n.num_site_rates ),
tau( n.tau ),
activeLikelihood( n.activeLikelihood ),
historyLikelihoods( n.historyLikelihoods ),
char_matrix( n.char_matrix ),
gap_matrix( n.gap_matrix ),
histories( n.histories ),
tipProbs( n.tipProbs ),
changed_nodes( n.changed_nodes ),
dirty_nodes( n.dirty_nodes ),
using_ambiguous_characters( n.using_ambiguous_characters ),
treatUnknownAsGap( n.treatUnknownAsGap ),
treatAmbiguousAsGaps( n.treatAmbiguousAsGaps ),
tipsInitialized( n.tipsInitialized )
{
    // We don'e want tau to die before we die, or it can't remove us as listener
    tau->getValue().getTreeChangeEventHandler().addListener( this );

}


/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
template<class charType>
RevBayesCore::AbstractTreeHistoryCtmc<charType>::~AbstractTreeHistoryCtmc( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

    // remove myself from the tree listeners
    if ( tau != NULL )
    {

        tau->getValue().getTreeChangeEventHandler().removeListener( this );

    }

}


//template<class charType>
//RevBayesCore::AbstractTreeHistoryCtmc<charType>* RevBayesCore::AbstractTreeHistoryCtmc<charType>::clone( void ) const
//{
//    
//    return new AbstractTreeHistoryCtmc<charType>( *this );
//}

template<class charType>
double RevBayesCore::AbstractTreeHistoryCtmc<charType>::computeLnProbability( void )
{

    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        dirty_nodes = std::vector<bool>(tau->getValue().getNumberOfNodes(), true);
    }

    this->lnProb = 0.0;

    const std::vector<TopologyNode*>& nodes = tau->getValue().getNodes();

    for (size_t i = 0; i < nodes.size(); i++)
    {
        const TopologyNode& nd = *nodes[i];
        size_t node_index = nd.getIndex();
        fillLikelihoodVector(nd);
        double nodeLnProb = historyLikelihoods[ activeLikelihood[node_index] ][node_index];
//        std::cout << nodeLnProb << "\n";
        this->lnProb += nodeLnProb;
    }

    return this->lnProb;
}


template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::fillLikelihoodVector(const TopologyNode &node)
{
    size_t node_index = node.getIndex();
    if (!dirty_nodes[node_index])
        return;

    // compute
    double lnL = computeInternalNodeLikelihood(node);

    if (node.isTip())
        lnL += computeTipLikelihood(node);

    if (node.isRoot())
        lnL += computeRootLikelihood(node);
    
    historyLikelihoods[ activeLikelihood[node_index] ][node_index] = lnL;
    
    // mark as computed
    dirty_nodes[node_index] = false;

}


template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::fireTreeChangeEvent( const RevBayesCore::TopologyNode &n, const unsigned& m ) {
    
    // call a recursive flagging of all node above (closer to the root) and including this node
    flagNodeDirty(n);

//    size_t idx = n.getIndex();
//    std::cout << "fireTreeChangeEvent() " << idx << "  " << (changed_nodes[idx] ? "1" : "0") << (dirty_nodes[idx] ? "1" : "0") << "\n";
}

template<class charType>
const RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType>::getHistory(const TopologyNode& nd) const
{
    return histories[nd.getIndex()];
}

template<class charType>
RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType>::getHistory(const TopologyNode& nd)
{
    return *histories[nd.getIndex()];
}


//template<class charType>
//const RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType>::getHistory(size_t idx) const
//{
//
//    return histories[idx];
//}
//
//template<class charType>
//RevBayesCore::BranchHistory&  RevBayesCore::AbstractTreeHistoryCtmc<charType>::getHistory(size_t idx)
//{
//    return *histories[idx];
//}

template<class charType>
const std::vector<RevBayesCore::BranchHistory*>& RevBayesCore::AbstractTreeHistoryCtmc<charType>::getHistories(void) const
{
    return histories;
}

template<class charType>
std::vector<RevBayesCore::BranchHistory*> RevBayesCore::AbstractTreeHistoryCtmc<charType>::getHistories(void)
{
    return histories;
}


template<class charType>
const std::vector<std::vector<double> >& RevBayesCore::AbstractTreeHistoryCtmc<charType>::getTipProbs(void)
{
    return tipProbs;
}

template<class charType>
const std::vector<double>& RevBayesCore::AbstractTreeHistoryCtmc<charType>::getTipProbs(const TopologyNode& nd)
{
    return tipProbs[nd.getIndex()];
}

template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::setTipProbs(const HomologousCharacterData* tp)
{

    tipProbs.clear();

    size_t num_taxa = tp->getNumberOfTaxa();
    size_t numCharacters = tp->getNumberOfCharacters();

    const std::vector<TopologyNode*>& nodes = this->tau->getValue().getNodes();

    tipProbs.resize(num_taxa);
    const ContinuousCharacterData* ccdp = static_cast<const ContinuousCharacterData*>(tp);
    for (size_t i = 0; i < nodes.size(); i++)
    {
        TopologyNode* nd = nodes[i];
        if (!nd->isTip())
            continue;

        const ContinuousTaxonData* cd = &ccdp->getTaxonData(nd->getName());
        for (size_t j = 0; j < numCharacters; j++)
        {
            double v = cd->getCharacter(j);
            //tipProbs[nd->getIndex()].push_back(1-v);
            tipProbs[nd->getIndex()].push_back(v);
        }
    }
}



template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::initializeHistoriesVector( void ) {

    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
    histories.resize(nodes.size());
    for (size_t i = 0; i < nodes.size(); i++)
    {
        TopologyNode* nd = nodes[i];
        histories[nd->getIndex()] = new BranchHistory(num_sites,num_chars,nd->getIndex());
    }

    historyLikelihoods.resize(2);
    for (size_t i = 0; i < 2; i++)
        historyLikelihoods[i].resize(nodes.size(), 0.0);

}


template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::keepSpecialization( const DagNode* affecter ) {

    // reset all flags
    for (std::vector<bool>::iterator it = this->dirty_nodes.begin(); it != this->dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
    for (std::vector<bool>::iterator it = this->changed_nodes.begin(); it != this->changed_nodes.end(); ++it)
    {
        (*it) = false;
    }

}



template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::flagNodeDirty( const RevBayesCore::TopologyNode &n ) {

    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();

    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( !dirty_nodes[index] )
    {
        // set the flag
        dirty_nodes[index] = true;
        
        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( !changed_nodes[index] )
        {
            activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
            //activeLikelihood[index] = 0;
            changed_nodes[index] = true;
        }

    }

}

template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::restoreSpecialization( const DagNode* affecter ) {

    // reset the flags
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = false;
    }

    //std::cout << "affecter " << affecter->getName() << "\n";

    // restore the active likelihoods vector
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changed_nodes[index] )
        {

            activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
            //activeLikelihood[index] = 0;
//            if (affecter->getName() == "ctmc") std::cout << index << " " << activeLikelihood[index] << "\n";
        }

        // set all flags to false
        changed_nodes[index] = false;
    }


    return;
}

template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::setHistory(const BranchHistory& bh, const TopologyNode& nd)
{
    histories[ nd.getIndex() ] = new BranchHistory(bh);
}


template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::setHistories(const std::vector<BranchHistory*>& bh)
{
    for (size_t i = 0; i < bh.size(); i++)
        histories[i] = bh[i];

}

template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::setValue(AbstractHomologousDiscreteCharacterData *v, bool force)
{

    // delegate to the parent class
    TypedDistribution< AbstractHomologousDiscreteCharacterData >::setValue(v, force);


    drawInitValue();

    if ( this->dag_node != NULL )
    {
        this->dag_node->getLnProbability();
    }

}


template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::simulate(void)
{

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new HomologousDiscreteCharacterData<charType>();

    // create a vector of taxon data
    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData< charType > >( tau->getValue().getNumberOfNodes(), DiscreteTaxonData<charType>( Taxon("") ) );

    // recursively simulate, starting with the root heading tipwards
    const TopologyNode& nd = tau->getValue().getRoot();
    histories[ nd.getIndex() ] = new BranchHistory(num_sites, num_chars, nd.getIndex());
    BranchHistory* bh = histories[ nd.getIndex() ];

    simulate(nd, bh, taxa);

    // add the taxon data to the character data
    for (size_t i = 0; i < tau->getValue().getNumberOfTips(); ++i)
    {
        this->value->addTaxonData( taxa[i] );
//        this->value->getTaxonData( tau->getValue().getNodes()[i]->getName() );
    }

    TypedDistribution< AbstractHomologousDiscreteCharacterData >::setValue(this->value);
}


/**
 * Swap a parameter of the distribution. We receive this call just before being replaced by a variable,
 * in which case the variable deletes the old parameter. We also receive this call during the cloning of
 * a DAG. Also in that case it is safe to leave the memory management of the tau parameter to others,
 * namely to the destructor of the original distribution owning the old parameter.
 */
template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP) {

    // we only have the topology here as the parameter
    if (oldP == tau)
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
        tau = static_cast<const TypedDagNode<Tree>* >( newP );
        tau->getValue().getTreeChangeEventHandler().addListener( this );
    }

}


template<class charType>
void RevBayesCore::AbstractTreeHistoryCtmc<charType>::touchSpecialization( const DagNode* affecter, bool touchAll )
{

    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == this->dag_node )
    {
        // do nothing, assume tree events have been fired
        ;
    }
    else if ( affecter != tau )
    {
        
        for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
        {
            (*it) = true;
        }

        // flip the active likelihood pointers
        for (size_t index = 0; index < changed_nodes.size(); ++index)
        {
            if ( !changed_nodes[index] )
            {
                activeLikelihood[index] = (activeLikelihood[index] == 0 ? 1 : 0);
                //activeLikelihood[index] = 0;
                changed_nodes[index] = true;
            }
        }
    }
}

#endif /* defined(__rb_mlandis__TreeHistory__) */
