#ifndef PhyloCTMCClado_H
#define PhyloCTMCClado_H

#include "AbstractCladogenicStateFunction.h"
#include "CharacterHistory.h"
#include "ChromosomesCladogenicStateFunction.h"
#include "CladogeneticProbabilityMatrix.h"
#include "AbstractPhyloCTMCSiteHomogeneous.h"
#include "BiogeographicCladoEvent.h"
#include "DistributionExponential.h"
#include "RateMatrix.h"
#include "RbBitSet.h"
#include "RbException.h"
#include "RbVector.h"
#include "Simplex.h"
#include "MatrixReal.h"
#include "Taxon.h"
#include "Tree.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "TypedDistribution.h"
#include "RandomNumberGenerator.h"

namespace RevBayesCore {

    template<class charType>
    class PhyloCTMCClado : public AbstractPhyloCTMCSiteHomogeneous<charType> {

    public:
//        AbstractPhyloCTMCSiteHomogeneous(const TypedDagNode<Tree> *t, size_t nChars, size_t nMix, bool c, size_t nSites, bool amb, bool wd = false, bool internal = false, bool gapmatch = true );
        PhyloCTMCClado(const TypedDagNode< Tree > *t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch);
        PhyloCTMCClado(const PhyloCTMCClado &n);
        virtual                                            ~PhyloCTMCClado(void);                                                                   //!< Virtual destructor

        // public member functions
        PhyloCTMCClado*                                     clone(void) const;                                                                          //!< Create an independent clone
        virtual double                                      computeLnProbability(void);
        virtual std::vector<charType>						drawAncestralStatesForNode(const TopologyNode &n);
        virtual void                                        drawJointConditionalAncestralStates(std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates);
        virtual void                                        recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates, const std::vector<size_t>& sampledSiteRates);

        virtual void                                        redrawValue(void);
        void                                                setCladogenesisMatrix(const TypedDagNode< CladogeneticProbabilityMatrix > *r);
        void                                                setCladogenesisMatrix(const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* r);
        void                                                setCladogenesisTimes(const TypedDagNode< RbVector< RbVector< double > > >* rm);

    protected:

        virtual void                                        resizeLikelihoodVectors(void);

        void                                                computeRootLikelihood(size_t root, size_t l, size_t r);
        void                                                computeRootLikelihood(size_t root, size_t l, size_t r, size_t m);
        void                                                computeInternalNodeLikelihood(const TopologyNode &n, size_t nIdx, size_t l, size_t r);
        void                                                computeTipLikelihood(const TopologyNode &node, size_t nIdx);
        void                                                updateTransitionProbabilities(size_t node_idx);

        virtual void                                        computeMarginalNodeLikelihood(size_t node_idx, size_t parentIdx);
        virtual void                                        computeMarginalRootLikelihood();
        virtual std::vector< std::vector< double > >        sumMarginalLikelihoods(size_t node_index);
        
        virtual void                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);            //!< Swap a parameter

        // the likelihoods
        std::vector<double>                                 cladoPartialLikelihoods;
        std::vector<double>                                 cladoMarginalLikelihoods;

        // offsets for nodes
        size_t                                              cladoActiveLikelihoodOffset;
        size_t                                              cladoNodeOffset;
        size_t                                              cladoMixtureOffset;
        size_t                                              cladoSiteOffset;


    private:
        virtual void                                            simulateClado(const TopologyNode& node, std::vector< DiscreteTaxonData< charType > > &t, const std::vector<size_t> &perSiteRates);
        virtual double                                          sumRootLikelihood( void );
        void                                                    updateTransitionProbabilityMatrices(void);
        
        const TypedDagNode< CladogeneticProbabilityMatrix >*                       homogeneousCladogenesisMatrix;
        const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >*           heterogeneousCladogenesisMatrices;
        const TypedDagNode< RbVector< RbVector< double > > >*   cladogenesisTimes;

        bool useObservedCladogenesis = false;
        bool useSampledCladogenesis = false;
        bool branchHeterogeneousCladogenesis = false;
        bool store_internal_nodes;
        bool gap_match_clamped;
    };

}

#include "AbstractCharacterHistoryBirthDeathProcess.h"
#include "BranchHistory.h"
#include "ConstantNode.h"
#include "StochasticNode.h"
#include "DECCladogeneticStateFunction.h"
#include "CladogeneticProbabilityMatrixFunction.h"
#include "DiscreteCharacterState.h"
#include "RateMatrix_JC.h"
#include "RandomNumberFactory.h"

#include <cmath>
#include <cstring>
#include <map>
#include <vector>

//        AbstractPhyloCTMCSiteHomogeneous(const TypedDagNode<Tree> *t, size_t nChars, size_t nMix, bool c, size_t nSites, bool amb, bool wd = false, bool internal = false, bool gapmatch = true );

template<class charType>
RevBayesCore::PhyloCTMCClado<charType>::PhyloCTMCClado(const TypedDagNode<Tree> *t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch)
    : AbstractPhyloCTMCSiteHomogeneous<charType>(  t, nChars, 1, c, nSites, amb, false, false, true ),
      store_internal_nodes(internal),
      gap_match_clamped(gapmatch)
{
//    unsigned numReducedChar = (unsigned)( log( nChars ) / log( 2 ) );
//    std::vector<std::string> et;
//    et.push_back("s");
//    et.push_back("a");
    
//    const TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >* events, const TypedDagNode<RevBayesCore::RbVector<double> >* probs, int n_states 
    
    // create a dummy matrix of identical cladogenetic inheritance triplets
    RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> >* clado_events_mtx_tmp = new RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> >();
    for (size_t i = 0; i < nChars; i++) {
        clado_events_mtx_tmp->push_back( RevBayesCore::RbVector<std::int64_t>(3, i) );
    }
    
    // populate a dummy node with those events
    TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<std::int64_t> > >* clado_events_tmp;
    clado_events_tmp = new RevBayesCore::ConstantNode< RevBayesCore::RbVector< RevBayesCore::RbVector<std::int64_t> > >(".cladogenetic_events", clado_events_mtx_tmp);
    
    // populate a dummy event probs vector where each event has prob = 1
    TypedDagNode< RevBayesCore::RbVector<double> >* clado_probs_tmp;
    clado_probs_tmp = new RevBayesCore::ConstantNode<RevBayesCore::RbVector<double> >(".probabilities", new RevBayesCore::RbVector<double>(clado_events_mtx_tmp->size(), 1.0));
    
    // build a dummy cladogenetic probability matrix function
    CladogeneticProbabilityMatrixFunction* clado_func_tmp = new CladogeneticProbabilityMatrixFunction( clado_events_tmp, clado_probs_tmp, (int)nChars );
    
    // place that function into a dummy deterministic node
    DeterministicNode<CladogeneticProbabilityMatrix>* clado_node_tmp = new DeterministicNode<CladogeneticProbabilityMatrix>( "cladogenesisMatrix", clado_func_tmp );
    
    // finally, assign our dummy node to the model's cladogenetic probability matrix variable
    homogeneousCladogenesisMatrix = clado_node_tmp;
   
    /*
    homogeneousCladogenesisMatrix            = new DeterministicNode<CladogeneticProbabilityMatrix>( "cladogenesisMatrix",
                                               new DECCladogeneticStateFunction(
                                                    new ConstantNode<Simplex>( "", new Simplex(2, 0.5)),
                                                    new ConstantNode<RbVector<RbVector<double> > >("", new RbVector<RbVector<double> >(nChars, RbVector<double>(nChars, 1))),
                                                    new ConstantNode<RbVector<RbVector<double> > >("", new RbVector<RbVector<double> >(nChars, RbVector<double>(nChars, 1))),
                                                    numReducedChar,
                                                    2,
                                                    et)
                                                );
     */
    heterogeneousCladogenesisMatrices        = NULL;
    cladogenesisTimes                        = NULL;
 
    
    cladoActiveLikelihoodOffset      =  this->num_nodes*this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    cladoNodeOffset                  =                  this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    cladoMixtureOffset               =                                       this->num_patterns*this->num_chars*this->num_chars;
    cladoSiteOffset                  =                                                          this->num_chars*this->num_chars;
   
    // check if the tree is a stochastic node before getting its distribution
    if ( this->tau->isStochastic() )
    {
        if ( dynamic_cast<const AbstractCharacterHistoryBirthDeathProcess* >( &this->tau->getDistribution() ) != NULL )
            useSampledCladogenesis = true;
    }

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( homogeneousCladogenesisMatrix );
    this->addParameter( heterogeneousCladogenesisMatrices );
    this->addParameter( cladogenesisTimes );
    
}

template<class charType>
RevBayesCore::PhyloCTMCClado<charType>::PhyloCTMCClado(const PhyloCTMCClado &n) :
    AbstractPhyloCTMCSiteHomogeneous<charType>( n ),

    cladoPartialLikelihoods( n.cladoPartialLikelihoods ),
    cladoMarginalLikelihoods( n.cladoMarginalLikelihoods ),

    useObservedCladogenesis(n.useObservedCladogenesis),
    useSampledCladogenesis(n.useSampledCladogenesis),
    store_internal_nodes(n.store_internal_nodes),
    gap_match_clamped(n.gap_match_clamped),
    branchHeterogeneousCladogenesis(n.branchHeterogeneousCladogenesis)
{
    // initialize with default parameters
    homogeneousCladogenesisMatrix       = n.homogeneousCladogenesisMatrix;
    heterogeneousCladogenesisMatrices   = n.heterogeneousCladogenesisMatrices;
    cladogenesisTimes                   = n.cladogenesisTimes;
    
    cladoActiveLikelihoodOffset      =  this->num_nodes*this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    cladoNodeOffset                  =                  this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    cladoMixtureOffset               =                                       this->num_patterns*this->num_chars*this->num_chars;
    cladoSiteOffset                  =                                                          this->num_chars*this->num_chars;
}



template<class charType>
RevBayesCore::PhyloCTMCClado<charType>::~PhyloCTMCClado( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


template<class charType>
RevBayesCore::PhyloCTMCClado<charType>* RevBayesCore::PhyloCTMCClado<charType>::clone( void ) const {

    return new PhyloCTMCClado<charType>( *this );
}

template<class charType>
double RevBayesCore::PhyloCTMCClado<charType>::computeLnProbability( void )
{

    
    // if we are not in MCMC mode, then we need to (temporarily) allocate memory
    if ( this->in_mcmc_mode == false )
    {
        cladoPartialLikelihoods.resize(2*this->num_nodes*this->num_site_rates*this->num_sites*this->num_chars*this->num_chars);
    }
    
    double lnL = RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeLnProbability();
    
    // if we are not in MCMC mode, then we need to (temporarily) free memory
    if ( this->in_mcmc_mode == false )
    {
        // free the partial likelihoods
        cladoPartialLikelihoods.clear();
    }
    
    return lnL;
}



template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::computeRootLikelihood( size_t root, size_t left, size_t right)
{
    // get the root frequencies
    const std::vector<double>& f = this->getRootFrequencies();
    const TopologyNode& node = this->tau->getValue().getRoot();
    const std::map<std::vector<unsigned>, double> eventMapProbs =
        branchHeterogeneousCladogenesis
        ? heterogeneousCladogenesisMatrices->getValue()[root].getEventMap(node.getAge())
        : homogeneousCladogenesisMatrix->getValue().getEventMap(node.getAge());

    // bypass cladogenetic probs if it's a sampled ancestor
    const bool has_sampled_ancestor_child =
        node.getChild(0).isSampledAncestorTip() || node.getChild(1).isSampledAncestorTip();

    // get the pointers to the partial likelihoods of the left and right subtree
    auto& pl_left  = this->getPartialLikelihoodsForNode(left);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    auto& pl_root  = this->createEmptyPartialLikelihoodsForNode(root, pl_left.dims());
    assert(pl_root.dims() == pl_left.dims());
    assert(pl_root.dims() == pl_right.dims());

    // cache loop-invariant members as locals
    const size_t num_chars       = this->num_chars;
    const size_t num_site_rates  = this->num_site_rates;
    const size_t num_patterns    = this->num_patterns;
    const size_t mixture_offset  = this->mixtureOffset;
    const size_t site_offset     = this->siteOffset;

    double*       __restrict__ const p_root  = pl_root.likelihoods.data();
    const double* __restrict__ const p_left  = pl_left.likelihoods.data();
    const double* __restrict__ const p_right = pl_right.likelihoods.data();

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
    {
        // compute the per site probabilities
        for (size_t site = 0; site < num_patterns; ++site)
        {
            const size_t offset = mixture * mixture_offset + site * site_offset;
            double*       __restrict__ const p_site       = p_root  + offset;
            const double* __restrict__ const p_site_left  = p_left  + offset;
            const double* __restrict__ const p_site_right = p_right + offset;

            // zero the accumulator
            for (size_t i = 0; i < num_chars; ++i)
            {
                p_site[i] = 0.0;
            }

            if (!has_sampled_ancestor_child)
            {
                // cladogenetic probs for bifurcations
                for (const auto& entry : eventMapProbs)
                {
                    const std::vector<unsigned>& idx = entry.first;
                    const size_t c1 = idx[0];
                    const size_t c2 = idx[1];
                    const size_t c3 = idx[2];
                    const double pcl = entry.second;

                    p_site[c1] += p_site_left[c2] * p_site_right[c3] * pcl;
                }
            }
            else
            {
                // no cladogenetic probs for sampled ancestors
                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    p_site[c1] += p_site_left[c1] * p_site_right[c1];
                }
            }

            // multiply by root frequencies
            for (size_t i = 0; i < num_chars; ++i)
            {
                p_site[i] *= f[i];
            }
        }
    }

    this->scale(root, left, right);
}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle)
{
    throw RbException()<<"PhyloCTMCClado::computeRootLikelihood(root, left, right, middle): the root should never have three children for PhyloCTMCClado!";
}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::computeInternalNodeLikelihood(const TopologyNode &node, size_t node_index, size_t left, size_t right)
{
    const auto& eventMapProbs = branchHeterogeneousCladogenesis
        ? heterogeneousCladogenesisMatrices->getValue()[node_index].getEventMap(node.getAge())
        : homogeneousCladogenesisMatrix->getValue().getEventMap(node.getAge());

    // bypass cladogenetic probs if it's a sampled ancestor
    const bool has_sampled_ancestor_child =
        node.getChild(0).isSampledAncestorTip() || node.getChild(1).isSampledAncestorTip();

    // compute the transition probability matrix
    this->updateTransitionProbabilityMatrix(node_index);

    // get the pointers to the partial likelihoods for this node and the two descendant subtrees
    auto& pl_left  = this->getPartialLikelihoodsForNode(left);
    auto& pl_right = this->getPartialLikelihoodsForNode(right);
    auto& pl_node  = this->createEmptyPartialLikelihoodsForNode(node_index, pl_left.dims());
    assert(pl_node.dims() == pl_left.dims());
    assert(pl_node.dims() == pl_right.dims());

    // cache loop-invariant members as locals
    const size_t num_chars            = this->num_chars;
    const size_t num_site_rates       = this->num_site_rates;
    const size_t num_patterns         = this->num_patterns;
    const size_t mixture_offset       = this->mixtureOffset;
    const size_t site_offset          = this->siteOffset;
    const size_t clado_mixture_offset = this->cladoMixtureOffset;
    const size_t clado_site_offset    = this->cladoSiteOffset;

    double*       __restrict__ const p_node_base  = pl_node.likelihoods.data();
    const double* __restrict__ const p_left_base  = pl_left.likelihoods.data();
    const double* __restrict__ const p_right_base = pl_right.likelihoods.data();
    double*       __restrict__ const p_clado_base =
        this->cladoPartialLikelihoods.data()
        + this->activeLikelihood[node_index] * this->cladoActiveLikelihoodOffset
        + node_index * this->cladoNodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* __restrict__ const tp = this->pmatrices[node_index][mixture].theMatrix;

        // compute the per site probabilities
        for (size_t site = 0; site < num_patterns; ++site)
        {
            const size_t offset       = mixture * mixture_offset       + site * site_offset;
            const size_t clado_offset = mixture * clado_mixture_offset + site * clado_site_offset;

            double*       __restrict__ const p_site       = p_node_base  + offset;
            const double* __restrict__ const p_site_left  = p_left_base  + offset;
            const double* __restrict__ const p_site_right = p_right_base + offset;
            double*       __restrict__ const p_clado_site = p_clado_base + clado_offset;

            // zero the cladogenetic accumulator
            for (size_t i = 0; i < num_chars; ++i)
            {
                p_clado_site[i] = 0.0;
            }

            if (!has_sampled_ancestor_child)
            {
                // cladogenetic probs for bifurcations
                for (const auto& entry : eventMapProbs)
                {
                    const std::vector<unsigned>& idx = entry.first;
                    const size_t c1 = idx[0];
                    const size_t c2 = idx[1];
                    const size_t c3 = idx[2];
                    const double pcl = entry.second;

                    p_clado_site[c1] += p_site_left[c2] * p_site_right[c3] * pcl;
                }
            }
            else
            {
                // no cladogenetic probs for sampled ancestors
                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    p_clado_site[c1] += p_site_left[c1] * p_site_right[c1];
                }
            }

            // anagenetic step: transition from older-end state c0 to younger-end state c1
            for (size_t c0 = 0; c0 < num_chars; ++c0)
            {
                double sum_ana = 0.0;
                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    sum_ana += tp[c0 * num_chars + c1] * p_clado_site[c1];
                }
                p_site[c0] = sum_ana;
            }
        }
    }

    this->scale(node_index, left, right);

}


template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::computeMarginalNodeLikelihood( size_t node_index, size_t parentnode_index )
{
    // get cladogenic transition probs
    const TopologyNode& node = this->tau->getValue().getNode(node_index);
    const auto& eventMapProbs = branchHeterogeneousCladogenesis
        ? heterogeneousCladogenesisMatrices->getValue()[node_index].getEventMap(node.getAge())
        : homogeneousCladogenesisMatrix->getValue().getEventMap(node.getAge());

    // compute the transition probability matrix
    this->updateTransitionProbabilityMatrix(node_index);

    // cache loop-invariant members as locals
    const size_t num_chars            = this->num_chars;
    const size_t num_site_rates       = this->num_site_rates;
    const size_t num_patterns         = this->num_patterns;
    const size_t mixture_offset       = this->mixtureOffset;
    const size_t site_offset          = this->siteOffset;
    const size_t clado_mixture_offset = this->cladoMixtureOffset;
    const size_t clado_site_offset    = this->cladoSiteOffset;

    // get the pointers to the partial likelihoods and the marginal likelihoods
    const double* __restrict__ const p_node                       = this->getPartialLikelihoodsForNode(node_index).likelihoods.data();
    const double* __restrict__ const p_parent_node_marginal       = this->getMarginalLikelihoodsForNode(parentnode_index);
    double*       __restrict__ const p_node_marginal              = this->getMutableMarginalLikelihoodsForNode(node_index);

    const double* __restrict__ const p_clado_node =
        this->cladoPartialLikelihoods.data()
        + this->activeLikelihood[node_index] * this->cladoActiveLikelihoodOffset
        + node_index * this->cladoNodeOffset;
    const double* __restrict__ const p_clado_parent_node_marginal =
        this->cladoMarginalLikelihoods.data() + parentnode_index * this->cladoNodeOffset;
    double*       __restrict__ const p_clado_node_marginal =
        this->cladoMarginalLikelihoods.data() + node_index      * this->cladoNodeOffset;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* __restrict__ const tp = this->pmatrices[node_index][mixture].theMatrix;

        // iterate over all sites
        for (size_t site = 0; site < num_patterns; ++site)
        {
            const size_t offset       = mixture * mixture_offset       + site * site_offset;
            const size_t clado_offset = mixture * clado_mixture_offset + site * clado_site_offset;

            const double* __restrict__ const p_site                = p_node                       + offset;
            const double* __restrict__ const p_parent_site_marginal = p_parent_node_marginal      + offset;
            double*       __restrict__ const p_site_marginal       = p_node_marginal              + offset;
            // clado pointers are set up but currently unused; see note below
            // const double* __restrict__ const p_clado_site                = p_clado_node                    + clado_offset;
            // const double* __restrict__ const p_clado_parent_site_marginal = p_clado_parent_node_marginal   + clado_offset;
            // double*       __restrict__ const p_clado_site_marginal       = p_clado_node_marginal           + clado_offset;

            // anagenetic step: marginalize over start state k for each end state j
            for (size_t j = 0; j < num_chars; ++j)
            {
                double sum = 0.0;
                for (size_t k = 0; k < num_chars; ++k)
                {
                    sum += p_parent_site_marginal[k] * tp[k * num_chars + j];
                }
                p_site_marginal[j] = p_site[j] * sum;
            }

            // TODO: cladogenetic step
            // iterate over all (X_L, X_R) states, after cladogenesis
            for (const auto& entry : eventMapProbs)
            {
                (void)entry;
                // unimplemented
            }
        }
    }

    throw RbException()<<"computeMarginalNodeLikelihood: unimplemented!";
}



template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::computeMarginalRootLikelihood( void )
{
    // get the root node
    const TopologyNode& root = this->tau->getValue().getRoot();
    const size_t node_index = root.getIndex();

    // get the root frequencies
    const std::vector<double>& f = this->getRootFrequencies();

    // cache loop-invariant members as locals
    const size_t num_chars       = this->num_chars;
    const size_t num_site_rates  = this->num_site_rates;
    const size_t num_patterns    = this->num_patterns;
    const size_t mixture_offset  = this->mixtureOffset;
    const size_t site_offset     = this->siteOffset;

    assert(f.size() == num_chars);

    // get the pointers to the partial likelihoods and the marginal likelihoods
    const double* __restrict__ const p_node          = this->getPartialLikelihoodsForNode(node_index).likelihoods.data();
    double*       __restrict__ const p_node_marginal = this->getMutableMarginalLikelihoodsForNode(node_index);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
    {
        // iterate over all sites
        for (size_t site = 0; site < num_patterns; ++site)
        {
            const size_t offset = mixture * mixture_offset + site * site_offset;

            for (size_t j = 0; j < num_chars; ++j)
            {
                p_node_marginal[offset + j] = p_node[offset + j] * f[j];
            }
        }
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::computeTipLikelihood(const TopologyNode &node, size_t node_index)
{
    double* p_node = this->createEmptyPartialLikelihoodsForNode(
        node_index, {this->num_site_mixtures, this->pattern_block_size, this->num_chars}
        ).likelihoods.data();

    // get the current correct tip index in case the whole tree changed (after performing an empiricalTree Proposal)
    const size_t data_tip_index = this->taxon_name_2_tip_index_map[node.getName()];
    const std::vector<bool>&              gap_node      = this->gap_matrix[data_tip_index];
    const std::vector<std::uint64_t>&     char_node     = this->char_matrix[data_tip_index];
    const std::vector<RbBitSet>&          amb_char_node = this->ambiguous_char_matrix[data_tip_index];

    // compute the transition probabilities
    this->updateTransitionProbabilityMatrix(node_index);

    // cache loop-invariant members as locals
    const size_t num_chars          = this->num_chars;
    const size_t num_site_mixtures  = this->num_site_mixtures;
    const size_t block_size         = this->pattern_block_size;
    const size_t mixture_offset     = this->mixtureOffset;
    const size_t site_offset        = this->siteOffset;
    const bool   using_ambiguous    = this->using_ambiguous_characters;
    const bool   using_weighted     = this->using_weighted_characters;

    double* __restrict__ const p_node_r = p_node;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double* __restrict__ const tp = this->pmatrices[node_index][mixture].theMatrix;

        // iterate over all sites
        for (size_t site = 0; site < block_size; ++site)
        {
            const size_t offset = mixture * mixture_offset + site * site_offset;
            double* __restrict__ const p_site = p_node_r + offset;

            // is this site a gap?
            if (gap_node[site])
            {
                // since this is a gap we assume the actual state could have been any state
                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    p_site[c1] = 1.0;
                }
            }
            else if (using_ambiguous && !using_weighted)
            {
                const RbBitSet& val = amb_char_node[site];
                assert(val.size() == num_chars);

                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    const double* __restrict__ const d = tp + c1 * num_chars;
                    double tmp = 0.0;
                    for (size_t i = 0; i < num_chars; ++i)
                    {
                        if (val.test(i))
                        {
                            tmp += d[i];
                        }
                    }
                    p_site[c1] = tmp;
                }
            }
            else if (using_weighted)
            {
                const RbBitSet& val = amb_char_node[site];
                assert(val.size() == num_chars);
                const auto& weights = this->value->getCharacter(node_index, site).getWeights();

                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    const double* __restrict__ const d = tp + c1 * num_chars;
                    double tmp = 0.0;
                    for (size_t i = 0; i < num_chars; ++i)
                    {
                        if (val.test(i))
                        {
                            tmp += d[i] * weights[i];
                        }
                    }
                    p_site[c1] = tmp;
                }
            }
            else // no ambiguous characters in use
            {
                const std::uint64_t org_val = char_node[site];
                for (size_t c1 = 0; c1 < num_chars; ++c1)
                {
                    p_site[c1] = tp[c1 * num_chars + org_val];
                }
            }
        }
    }

    this->scale(node_index);
}


/**
 * Draw a vector of ancestral states from the marginal distribution (non-conditional of the other ancestral states).
 * Here we assume that the marginal likelihoods have been updated.
 */
template<class charType>
std::vector<charType> RevBayesCore::PhyloCTMCClado<charType>::drawAncestralStatesForNode(const TopologyNode &node)
{
	
    size_t node_index = node.getIndex();
	
    // get the marginal likelihoods
    const auto marginals = sumMarginalLikelihoods(node_index);
    
    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector< charType > ancestralSeq = std::vector<charType>();
	
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        size_t pattern = i;
        // if the matrix is compressed use the pattern for this site
        if (this->compressed)
        {
            pattern = this->site_pattern[i];
        }

        // create the character
        charType c = charType( this->template_state );
        
        // sum the likelihoods for each character state
        auto& siteMarginals = marginals[pattern];
        double sumMarginals = 0.0;
        for (int j = 0; j < siteMarginals.size(); j++)
        {
            sumMarginals += siteMarginals[j];
        }

        double u = rng->uniform01();
        if (sumMarginals == 0.0)
        {

            // randomly draw state if all states have 0 probability
            c.setStateByIndex((size_t)(u*c.getNumberOfStates()));

        }
        else
        {

            // the marginals don't add up to 1, so rescale u
            u *= sumMarginals;

            // draw the character state
            size_t stateIndex = 0;
            while ( true )
            {

                u -= siteMarginals[stateIndex];

                if ( u > 0.0 )
                {
                    stateIndex++;

                    if ( stateIndex == c.getNumberOfStates() )
                    {
                        stateIndex = 0;
                        c.setToFirstState();
                    }

                    else
                    {
                        c++;
                    }
                }
                else
                {
                    break;
                }
            }
        }

        // add the character to the sequence
        ancestralSeq.push_back( c );
    }

    return ancestralSeq;
}


/**
 * Draw a vector of ancestral states from the joint-conditional distribution of states.
 */
template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::drawJointConditionalAncestralStates(std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    this->sampled_site_mixtures.resize(this->num_sites);

    const TopologyNode& root = this->tau->getValue().getRoot();
    const size_t node_index = root.getIndex();
    const size_t right      = root.getChild(0).getIndex();
    const size_t left       = root.getChild(1).getIndex();

    // get working variables
    const std::vector<double>& f = this->getRootFrequencies();
    std::vector<double> siteProbVector(1, 1.0);
    if (this->site_rates_probs != NULL)
    {
        siteProbVector = this->site_rates_probs->getValue();
    }

    // get cladogenesis values
    const auto& eventMapProbs = branchHeterogeneousCladogenesis
        ? heterogeneousCladogenesisMatrices->getValue()[node_index].getEventMap(root.getAge())
        : homogeneousCladogenesisMatrix->getValue().getEventMap(root.getAge());

    // cache loop-invariant members as locals
    const size_t num_sites      = this->num_sites;
    const size_t num_site_rates = this->num_site_rates;
    const size_t site_offset    = this->siteOffset;
    const size_t mixture_offset = this->mixtureOffset;

    const double* __restrict__ const p_node  = this->getPartialLikelihoodsForNode(node_index).likelihoods.data();
    const double* __restrict__ const p_left  = this->getPartialLikelihoodsForNode(left).likelihoods.data();
    const double* __restrict__ const p_right = this->getPartialLikelihoodsForNode(right).likelihoods.data();

    std::vector<size_t> sampledSiteRates(num_sites, 0);

    // scratch buffer reused across sites to avoid reallocating the map each iteration
    std::map<std::vector<unsigned>, double> sampleProbs;

    for (size_t i = 0; i < num_sites; ++i)
    {
        sampleProbs.clear();
        double sum = 0.0;

        // if the matrix is compressed use the pattern for this site
        const size_t pattern = this->compressed ? this->site_pattern[i] : i;
        const size_t site_base = pattern * site_offset;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
        {
            const size_t offset = site_base + mixture * mixture_offset;
            const double mixture_prob = siteProbVector[mixture];

            // iterate over all cladogenetic (A, L, R) triplets
            for (const auto& [v, pcl] : eventMapProbs)
            {
                const double prob = p_node [offset + v[0]]
                    * p_left [offset + v[1]]
                    * p_right[offset + v[2]]
                    * f[v[0]]
                    * mixture_prob
                    * pcl;

                std::vector<unsigned> key = v;
                key.push_back(static_cast<unsigned>(mixture));
                sampleProbs[std::move(key)] = prob;
                sum += prob;
            }
        }

        // sample from the cumulative distribution
        double u = rng->uniform01() * sum;
        for (const auto& entry : sampleProbs)
        {
            u -= entry.second;
            if (u < 0.0)
            {
                const std::vector<unsigned>& v = entry.first;

                charType ca(this->template_state);
                charType cl(this->template_state);
                charType cr(this->template_state);
                ca += v[0];
                cl += v[1];
                cr += v[2];

                endStates  [node_index][i] = ca;
                startStates[node_index][i] = ca;
                startStates[left       ][i] = cl;
                startStates[right      ][i] = cr;
                sampledSiteRates[i] = v[3];
                break;
            }
        }
    }

    // recurse
    std::vector<TopologyNode*> children = root.getChildren();
    for (size_t i = 0; i < children.size(); ++i)
    {
        if (!children[i]->isTip())
        {
            recursivelyDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampledSiteRates);
        }
        else
        {
            AbstractPhyloCTMCSiteHomogeneous<charType>::tipDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampledSiteRates);
        }
    }
    
}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates, const std::vector<size_t>& sampledSiteRates)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get working variables
    size_t node_index = node.getIndex();
    size_t left = node.getChild(0).getIndex();
    size_t right = node.getChild(1).getIndex();
    
    std::map<std::vector<unsigned>, double> eventMapProbs = ( branchHeterogeneousCladogenesis ? heterogeneousCladogenesisMatrices->getValue()[node_index].getEventMap(node.getAge()) : homogeneousCladogenesisMatrix->getValue().getEventMap(node.getAge()) );
    
    std::map<std::vector<unsigned>, double> sampleProbs;
    std::map<std::vector<unsigned>, double>::iterator it_s;
    std::map<std::vector<unsigned>, double>::iterator it_p;

    // get transition probabilities
    this->updateTransitionProbabilityMatrix( node_index );
    
    // get the pointers to the partial likelihoods and the marginal likelihoods
    const double*   p_left  = this->getPartialLikelihoodsForNode(left).likelihoods.data();
    const double*   p_right = this->getPartialLikelihoodsForNode(right).likelihoods.data();

    // sample characters conditioned on start states, going to end states
    std::vector<double> p(this->num_chars, 0.0);
    for (size_t i = 0; i < this->num_sites; i++)
    {
        size_t cat = sampledSiteRates[i];
        size_t k = startStates[node_index][i].getStateIndex();
        
        // sum to sample
        double sum = 0.0;

		// if the matrix is compressed use the pattern for this site
        size_t pattern = i;
		if (this->compressed)
        {
			pattern = this->site_pattern[i];
		}

        const double* p_left_site_mixture  = p_left  + cat * this->mixtureOffset + pattern * this->siteOffset;
        const double* p_right_site_mixture = p_right + cat * this->mixtureOffset + pattern * this->siteOffset;

        // iterate over possible end-anagenesis states for each site given start-anagenesis state
        for (it_p = eventMapProbs.begin(); it_p != eventMapProbs.end(); it_p++)
        {
            // triplet of (A,L,R) states
            const std::vector<unsigned>& v = it_p->first;

            const double* p_left_site_mixture_j  = p_left_site_mixture  + v[1];
            const double* p_right_site_mixture_j = p_right_site_mixture + v[2];

            // anagenesis prob
            size_t j = v[0];
            double tp_kj = this->pmatrices[node_index][cat][k][j];
            
            // anagenesis + cladogenesis prob
            sampleProbs[ it_p->first ] = it_p->second * tp_kj * *p_left_site_mixture_j * *p_right_site_mixture_j;
            sum += sampleProbs[ it_p->first ];

        }

        // sample char from p
        charType ca = charType( this->template_state );
        charType cl = charType( this->template_state );
        charType cr = charType( this->template_state );
        double u = rng->uniform01() * sum;
        for (it_s = sampleProbs.begin(); it_s != sampleProbs.end(); it_s++)
        {
            u -= it_s->second;
            if (u < 0.0)
            {
                const std::vector<unsigned>& v = it_s->first;
                ca += v[0];
                cl += v[1];
                cr += v[2];
                endStates[node_index][i] = ca;
                startStates[left][i] = cl;
                startStates[right][i] = cr;
                break;
            }
        }
    }

    // recurse
    std::vector<TopologyNode*> children = node.getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        // recurse towards tips
        if (!children[i]->isTip() == true )
        {
            recursivelyDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampledSiteRates);
        }
        else
        {
            AbstractPhyloCTMCSiteHomogeneous<charType>::tipDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampledSiteRates);
        }
        
    }
    
}


template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::resizeLikelihoodVectors( void )
{
    // call base resize
    RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::resizeLikelihoodVectors();

    size_t n = this->num_nodes*this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    
    // only do this if we are in MCMC mode. This will safe memory
    if ( this->in_mcmc_mode == true )
    {
        
        // we resize the partial likelihood vectors to the new dimensions
        cladoPartialLikelihoods.clear();
        
        cladoPartialLikelihoods.resize(2*n);
        
        // reinitialize likelihood vectors
        for (size_t i = 0; i < 2*n; i++)
        {
            cladoPartialLikelihoods[i] = 0.0;
        }
        
    }
    
    if ( this->useMarginalLikelihoods == true )
    {
        // we resize the partial likelihood vectors to the new dimensions
        cladoMarginalLikelihoods.clear();
        
        cladoMarginalLikelihoods.resize(n);
        
        // reinitialize likelihood vectors
        for (size_t i = 0; i < n; i++)
        {
            cladoMarginalLikelihoods[i] = 0.0;
        }
        
    }
	
    // set the offsets for easier iteration through the likelihood vector
    cladoActiveLikelihoodOffset      =  this->num_nodes*this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    cladoNodeOffset                  =                  this->num_site_rates*this->num_patterns*this->num_chars*this->num_chars;
    cladoMixtureOffset               =                                       this->num_patterns*this->num_chars*this->num_chars;
    cladoSiteOffset                  =                                                          this->num_chars*this->num_chars;

}


template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::setCladogenesisMatrix(const TypedDagNode< CladogeneticProbabilityMatrix > *cm) {

    // remove the old parameter first
    if ( homogeneousCladogenesisMatrix != NULL )
    {
        this->removeParameter( homogeneousCladogenesisMatrix );
        homogeneousCladogenesisMatrix = NULL;
    }
    else
    {
        this->removeParameter( heterogeneousCladogenesisMatrices );
        heterogeneousCladogenesisMatrices = NULL;
    }

    // set the value
    branchHeterogeneousCladogenesis = false;
    useObservedCladogenesis = true;
    homogeneousCladogenesisMatrix = cm;

    // add the new parameter
    this->addParameter( homogeneousCladogenesisMatrix );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::setCladogenesisMatrix(const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > > *cm) {
    
    // remove the old parameter first
    if ( homogeneousCladogenesisMatrix != NULL )
    {
        this->removeParameter( homogeneousCladogenesisMatrix );
        homogeneousCladogenesisMatrix = NULL;
    }
    else
    {
        this->removeParameter( heterogeneousCladogenesisMatrices );
        heterogeneousCladogenesisMatrices = NULL;
    }

    // set the value
    branchHeterogeneousCladogenesis = true;
    useObservedCladogenesis = true;
    heterogeneousCladogenesisMatrices = cm;

    // add the new parameter
    this->addParameter( heterogeneousCladogenesisMatrices );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::setCladogenesisTimes(const TypedDagNode< RbVector< RbVector< double > > > *ct)
{

    if (cladogenesisTimes != NULL)
    {
        this->removeParameter( cladogenesisTimes );
        cladogenesisTimes = NULL;
    }

    // set the value
    useSampledCladogenesis = true;
    cladogenesisTimes = ct;

    // add the new parameter
    this->addParameter( cladogenesisTimes );

    // redraw the current value
    if ( this->dagNode == NULL || !this->dag_node->isClamped() )
    {
        this->redrawValue();
    }

}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::simulateClado( const TopologyNode &node, std::vector< DiscreteTaxonData< charType > > &taxa, const std::vector<size_t> &perSiteRates)
{
//    std::cout << "SIMULATE\n";
    // first simulate cladogenic changes
    if (node.getNumberOfChildren() > 2)
    {
        throw RbException( "The tree is not bifurcating. Cannot simulate cladogenic evolution." );
    }
    
    // get cladogenesis event map (sparse transition probability matrix)
    std::map<std::vector<unsigned>, double> eventMapProbs = homogeneousCladogenesisMatrix->getValue().getEventMap(node.getAge());
    
    // get the character state of this node before cladogenic change
    size_t node_index = node.getIndex();
    
    const DiscreteTaxonData< charType > &parent = taxa[ node_index ];
    DiscreteTaxonData< charType > *left = new DiscreteTaxonData<charType>( Taxon("") );
    DiscreteTaxonData< charType > *right = new DiscreteTaxonData<charType>( Taxon("") );
    
    // simulate the left and right character states after cladogenic change
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        // this is the parent's state before clado change
        std::uint64_t parentState = parent.getCharacter( i ).getStateIndex();
        
        // simulate left and right states after clado changes
        charType cl = charType( this->num_chars );
        charType cr = charType( this->num_chars );
        cl.setToFirstState();
        cr.setToFirstState();
        double u = rng->uniform01();
        std::map<std::vector<unsigned>, double>::iterator it;
        for (it = eventMapProbs.begin(); it != eventMapProbs.end(); it++)
        {
            const std::vector<unsigned>& states = it->first;
            
            if ( parentState == states[0] )
            {
                u -= it->second;
                if (u < 0.0)
                {
//                    std::cout << states[0] << " -> " << states[1] << " | " << states[2] << "\n";
                    cl += states[1];
                    cr += states[2];
                    left->addCharacter( cl );
                    right->addCharacter( cr );
                    break;
                }
            }
        }
        if (left->getNumberOfCharacters() == 0 && right->getNumberOfCharacters() == 0)
        {
            cl += (int) parentState;
            cr += (int) parentState;
            left->addCharacter( cl );
            right->addCharacter( cr );
        }
    }
    
    // now simulate anagenic changes
    
    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();
    bool first_child = true;
    
    // simulate the sequence for each child
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        updateTransitionProbabilities( child.getIndex() );

        DiscreteTaxonData< charType > &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < this->num_sites; ++i )
        {
            // get the parent's state after clado change
            std::uint64_t parentState;
            if (first_child)
            {
                parentState = left->getCharacter( i ).getStateIndex();
            }
            else
            {
                parentState = right->getCharacter( i ).getStateIndex();
            }

            // use the parent's end state to calculate anagenetic changes
            const double *freqs = this->pmatrices[ child.getIndex() ][ perSiteRates[i] ][ parentState ];
            
            // create the character
            charType c = charType( this->num_chars );
            c.setToFirstState();
            
            // draw the state
            double u = rng->uniform01();
            size_t stateIndex = 0;
            while ( true )
            {
                u -= *freqs;
                ++stateIndex;
                
                if ( u > 0.0 && stateIndex < this->num_chars)
                {
                    ++c;
                    ++freqs;
                }
                else
                {
                    break;
                }
                
            }
            
            // add the character to the sequence
            taxon.addCharacter( c );
        }

        if ( child.isTip() )
        {
            taxon.setTaxon( child.getTaxon() );
        }
        else
        {
            // recursively simulate the sequences
            simulateClado( child, taxa, perSiteRates );
        }
        first_child = false;
    }
    
    // clean up
    delete left;
    delete right;
}

template<class charType>
std::vector< std::vector<double> > RevBayesCore::PhyloCTMCClado<charType>::sumMarginalLikelihoods( size_t node_index )
{
    // cache loop-invariant members as locals
    const size_t num_chars       = this->num_chars;
    const size_t num_site_rates  = this->num_site_rates;
    const size_t num_patterns    = this->num_patterns;
    const size_t mixture_offset  = this->mixtureOffset;
    const size_t site_offset     = this->siteOffset;

    std::vector<std::vector<double>> per_mixture_Likelihoods(num_patterns, std::vector<double>(num_chars, 0.0));

    // get the pointer to the marginal likelihoods
    const double* __restrict__ const p_node_marginal = this->getMarginalLikelihoodsForNode(node_index);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
    {
        // iterate over all sites
        for (size_t site = 0; site < num_patterns; ++site)
        {
            const size_t offset = mixture * mixture_offset + site * site_offset;
            std::vector<double>& site_likelihoods = per_mixture_Likelihoods[site];

            for (size_t j = 0; j < num_chars; ++j)
            {
                site_likelihoods[j] += p_node_marginal[offset + j];
            }
        }
    }

    return per_mixture_Likelihoods;
}


template<class charType>
double RevBayesCore::PhyloCTMCClado<charType>::sumRootLikelihood( void )
{
    // get the root node
    const TopologyNode& root = this->tau->getValue().getRoot();
    const size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods
    const double* __restrict__ const p_node = this->getPartialLikelihoodsForNode(node_index).likelihoods.data();
    auto& scale_node = this->getPartialLikelihoodsForNode(node_index).scale;

    // cache loop-invariant members as locals
    const size_t num_chars       = this->num_chars;
    const size_t num_site_rates  = this->num_site_rates;
    const size_t num_patterns    = this->num_patterns;
    const size_t mixture_offset  = this->mixtureOffset;
    const size_t site_offset     = this->siteOffset;

    // create a vector for the per mixture likelihoods
    std::vector<double> per_mixture_Likelihoods(num_patterns, 0.0);

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < num_site_rates; ++mixture)
    {
        // iterate over all sites
        for (size_t site = 0; site < num_patterns; ++site)
        {
            const size_t offset = mixture * mixture_offset + site * site_offset;

            // sum the likelihoods across starting states
            double tmp = 0.0;
            for (size_t i = 0; i < num_chars; ++i)
            {
                tmp += p_node[offset + i];
            }

            per_mixture_Likelihoods[site] += tmp;
        }
    }    

    // sum the log-likelihoods for all sites together
    double sumPartialProbs = 0.0;
    // get the root frequencies
    const std::vector<double> &f = this->getRootFrequencies();
    
    double p_inv = this->p_inv->getValue();
    double oneMinusPInv = 1.0 - p_inv;
    std::vector< size_t >::const_iterator patterns = this->pattern_counts.begin();
    if ( p_inv > 0.0 )
    {
        for (size_t site = 0; site < this->num_patterns; ++site, ++patterns)
        {
            
            if ( RbSettings::userSettings().getUseScaling() == true )
            {
                
                if ( this->site_invariant[site] )
                {
                    double ftotal = 0.0;
                    for ( size_t c = 0; c < this->invariant_site_index[site].size(); c++ )
                    {
                        ftotal += f[this->invariant_site_index[site][c]];
                    }

                    sumPartialProbs += log( p_inv * ftotal * exp(scale_node[site] * log_scale_factor) + oneMinusPInv * per_mixture_Likelihoods[site] / this->num_site_rates ) * *patterns;
                }
                else
                {
                    sumPartialProbs += log( oneMinusPInv * per_mixture_Likelihoods[site] / this->num_site_rates ) * *patterns;
                }
                sumPartialProbs -= scale_node[site] * log_scale_factor * *patterns;
                
            }
            else // no scaling
            {
                
                if ( this->site_invariant[site] )
                {
                    double ftotal = 0.0;
                    for ( size_t c = 0; c < this->invariant_site_index[site].size(); c++ )
                    {
                        ftotal += f[this->invariant_site_index[site][c]];
                    }

                    sumPartialProbs += log( p_inv * ftotal + oneMinusPInv * per_mixture_Likelihoods[site] / this->num_site_rates ) * *patterns;
                }
                else
                {
                    sumPartialProbs += log( oneMinusPInv * per_mixture_Likelihoods[site] / this->num_site_rates ) * *patterns;
                }

            }
        }
    }
    else
    {
        
        for (size_t site = 0; site < this->num_patterns; ++site, ++patterns)
        {
            
            sumPartialProbs += log( per_mixture_Likelihoods[site] / this->num_site_rates ) * *patterns;
            
            if ( RbSettings::userSettings().getUseScaling() == true )
            {
                
                sumPartialProbs -= scale_node[site] * log_scale_factor * *patterns;
            }

        }


    }

    return sumPartialProbs;
}



/** Swap a parameter of the distribution */
template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == homogeneousCladogenesisMatrix)
    {
        homogeneousCladogenesisMatrix = static_cast<const TypedDagNode< CladogeneticProbabilityMatrix >* >( newP );
    }
    else if (oldP == heterogeneousCladogenesisMatrices)
    {
        heterogeneousCladogenesisMatrices = static_cast<const TypedDagNode< RbVector< CladogeneticProbabilityMatrix > >* >( newP );
    }
    else if (oldP == cladogenesisTimes)
    {
        cladogenesisTimes = static_cast<const TypedDagNode< RbVector< RbVector< double > > >* >( newP );
    }
    else
    {
        RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::swapParameterInternal(oldP, newP);
    }
}

template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::updateTransitionProbabilities(size_t node_idx)
{
    if (not this->pmatrices.is_dirty(node_idx)) return;

    // get cladogenesis event map (sparse transition probability matrix)
    const TopologyNode* node = this->tau->getValue().getNodes()[node_idx];
    if (node->isRoot()) throw RbException("ERROR: dnPhyloCTMC called updateTransitionProbabilities for the root node\n");

    // FIXME: Don't copy this!
    std::map<std::vector<unsigned>, double> eventMapProbs = homogeneousCladogenesisMatrix->getValue().getEventMap(node->getAge());
 
    // first, get the rate matrix for this branch
    RateMatrix_JC jc(this->num_chars);
    const RateGenerator *rm = &jc;

    if ( this->branch_heterogeneous_substitution_matrices == true )
    {
        if (this->heterogeneous_rate_matrices != NULL) {
            rm = &this->heterogeneous_rate_matrices->getValue()[node_idx];
        }
    }
    else
    {
        if (this->homogeneous_rate_matrix != NULL) {
            rm = &this->homogeneous_rate_matrix->getValue();
        }
    }
    
    // second, get the clock rate for the branch
    double rate = 1.0;
    if ( this->branch_heterogeneous_clock_rates == true )
    {
        if (this->heterogeneous_clock_rates != NULL) {
            rate = this->heterogeneous_clock_rates->getValue()[node_idx];
        }
    }
    else
    {
        if (this->homogeneous_clock_rate != NULL) {
            rate = this->homogeneous_clock_rate->getValue();
        }
    }

    // and finally compute the per site rate transition probability matrix
    double startAge = node->getParent().getAge();
    double endAge = node->getAge();
    
    // if the tree is not a time tree, then the age will be not a number
    if ( RbMath::isFinite(endAge) == false )
    {
        // we assume by default that the end is at time 0
        endAge = 0.0;
    }
    
    TransitionProbabilityMatrix tmpMatrix(this->num_chars);

    // get sampled cladogenic events for this branch
    if (useSampledCladogenesis)
    {
        
        // convert underlying tree type
        const AbstractCharacterHistoryBirthDeathProcess* dist = dynamic_cast<const AbstractCharacterHistoryBirthDeathProcess* >( &this->tau->getDistribution() );
        
        // get history information
        const CharacterHistory &tree_history = dist->getCharacterHistory();
        const BranchHistory &branch_history = tree_history[node_idx];
        const std::multiset<CharacterEvent*,CharacterEventCompare>& events = branch_history.getHistory();
        
        if (events.size() == 0)
        {
            RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateTransitionProbabilityMatrix(node_idx);
        }
        else
        {
            // get cladogenetic transition probs (assumes they are constant w/r/t age)
            TransitionProbabilityMatrix cp(this->num_chars);
            for (size_t i = 0; i < this->num_chars; i++)
                cp[i][i] = 0.0;
            cp[0][0] = 1.0;
            
            // first compute clado probs at younger end of branch
            for (std::map<std::vector<unsigned>, double>::iterator it = eventMapProbs.begin(); it != eventMapProbs.end(); ++it)
            {
                // sparse elements from map
                const std::vector<unsigned>& idx = it->first;
                const size_t i = idx[0];
                const size_t j = idx[1];
                double p_clado = it->second;
                cp[i][j] += p_clado;
            }
            
            
            TransitionProbabilityMatrix tp(this->num_chars);
            for (size_t i = 0; i < this->num_chars; i++)
                tp[i][i] = 1.0;
            
            // for each interval between events, go from present to past
            double t = 0.0;
            double dt = 0.0;
            double event_age = startAge;
            bool first_event = true;
            for (auto it = events.rbegin(); it != events.rend(); it++)
            {
                t += dt;
                dt = (*it)->getAge() - t;
                event_age = event_age - dt;
                
                // anagenetic changes occurring between (event_age, event_age-dt)
                rm->calculateTransitionProbabilities(event_age+dt, event_age, rate, tmpMatrix );
                
                if (first_event)
                {
                    tp = tmpMatrix;
                    first_event = false;
                }
                else
                {
                    tp *= tmpMatrix;
                }

                // cladogenetic component
                tp *= cp;
            }

            // last interval
            rm->calculateTransitionProbabilities( event_age, endAge,  rate, tmpMatrix );
            tp *= tmpMatrix;
            
            auto& pmat_mixture = this->pmatrices.init_for_writing(node_idx);
            pmat_mixture.resize(1, tp);
            pmat_mixture[0] = tp;
        }
    }
    else
    {
        RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateTransitionProbabilityMatrix(node_idx);
    }
}


template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::updateTransitionProbabilityMatrices( void )
{
    // doing nothing here as this function has only been implemented for the classes that
    // do not have its own updateTransitionProbabilities function.
    // making this an empty function for this class for now so no redundant computation would be incurred when
    // updateTransitionProbabilityMatrices gets called in computeLnProbability in AbstractPhyloCTMCSiteHomogeneous
}


template<class charType>
void RevBayesCore::PhyloCTMCClado<charType>::redrawValue( void )
{
    
    bool do_mask = this->dag_node != NULL && this->dag_node->isClamped() && gap_match_clamped;

    std::vector<std::vector<bool> > mask = std::vector<std::vector<bool> >(this->tau->getValue().getNumberOfTips(), std::vector<bool>());
    // we cannot use the stored gap matrix because it uses the pattern compression
    // therefore we create our own mask
    if ( do_mask == true )
    {
        // set the gap states as in the clamped data
        for (size_t i = 0; i < this->tau->getValue().getNumberOfTips(); ++i)
        {
            // create a temporary variable for the taxon
            std::vector<bool> taxon_mask = std::vector<bool>(this->num_sites,false);
            
            const std::string &taxon_name = this->tau->getValue().getNode( i ).getName();
            AbstractDiscreteTaxonData& taxon = this->value->getTaxonData( taxon_name );
            
            for ( size_t site=0; site<this->num_sites; ++site)
            {
                taxon_mask[site] = taxon.getCharacter( site ).isGapState();
            }
            
            mask[i] = taxon_mask;
        }
    }
    
    // delete the old value first
    delete this->value;
    
    // create a new character data object
    this->value = new HomologousDiscreteCharacterData<charType>();
    
    // create a vector of taxon data
    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData< charType > >( this->num_nodes, DiscreteTaxonData<charType>( Taxon("") ) );
    
    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector<size_t> perSiteRates = std::vector<size_t>(this->num_sites,0);
    std::vector<bool> inv = std::vector<bool>(this->num_sites,false);
    double prob_invariant = this->getPInv();
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        // draw if this site is invariant
        double u = rng->uniform01();
        if ( u < prob_invariant )
        {
            // this site is invariant
            inv[i] = true;
            
        }
        else if ( this->num_site_rates  > 1 )
        {
            // draw the rate for this site
            u = rng->uniform01();
            size_t rateIndex = size_t(u*this->num_site_rates);
            perSiteRates[i] = rateIndex;
            
        }
        else
        {
            // there is only a single site rate so this is 1.0
            perSiteRates[i] = 0;
            
        }
        
    }
    
    // simulate the root sequence
    std::vector<std::vector<double> > freqs;
    this->getRootFrequencies(freqs);
    // const std::vector< double > &stationary_freqs = this->getRootFrequencies();
    this->getRootFrequencies();
    DiscreteTaxonData< charType > &root = taxa[ this->tau->getValue().getRoot().getIndex() ];
    for ( size_t i = 0; i < this->num_sites; ++i )
    {
        const std::vector< double > &stationary_freqs = freqs[perSiteRates[i] % freqs.size()];
        
        // create the character
        charType c = charType( this->template_state );
        c.setToFirstState();

        // draw the state
        double u = rng->uniform01();
        std::vector< double >::const_iterator freq = stationary_freqs.begin();
        while ( true )
        {
            u -= *freq;
            
            if ( u > 0.0 )
            {
                 if (c.getStateIndex() + 1 >= c.getNumberOfStates())
                 {
                     c.setToFirstState();
                 }
                 else
                 {
                    ++c;
                 }
                ++freq;
            }
            else
            {
                break;
            }
            
        }
        
        // add the character to the sequence
        root.addCharacter( c );
    }

    // recursively simulate the sequences
    root.setTaxon( Taxon("Root") );
    
    // recursively simulate the sequences
    simulateClado( this->tau->getValue().getRoot(), taxa, perSiteRates );
    
    // add the taxon data to the character data
    for (size_t i = 0; i < this->tau->getValue().getNumberOfNodes(); ++i)
    {
        
        const TopologyNode& node = this->tau->getValue().getNode(i);
        size_t node_index = node.getIndex();
        
        if (node.isTip() == true)
        {
            this->value->addTaxonData( taxa[node_index] );
        }
        else if (store_internal_nodes == true)
        {
            std::stringstream ss;
            ss << "Index_" << node_index + 1;
            taxa[node_index].setTaxon( Taxon(ss.str()) );
            this->value->addTaxonData( taxa[node_index] );
        }
        
    }

    if ( do_mask == true )
    {
        // set the gap states as in the clamped data
        for (size_t i = 0; i < this->tau->getValue().getNumberOfTips(); ++i)
        {
            const std::string &taxon_name = this->tau->getValue().getNode( i ).getName();
            AbstractDiscreteTaxonData& taxon = this->value->getTaxonData( taxon_name );
            
            for ( size_t site=0; site<this->num_sites; ++site)
            {
                DiscreteCharacterState &c = taxon.getCharacter(site);
                if ( mask[i][site] == true )
                {
                    c.setGapState( true );
                }
            }
            
        }
        
    }
    
    // compress the data and initialize internal variables
    this->compress();
    
    this->markAllPartialLikelihoodsDirty();
    
    // flip the active likelihood pointers
    for (size_t index = 0; index < this->changed_nodes.size(); ++index)
    {
        if ( this->changed_nodes[index] == false )
        {
            this->activeLikelihood[index] = (this->activeLikelihood[index] == 0 ? 1 : 0);
            this->changed_nodes[index] = true;
        }
    }
    
}

#endif /* defined(__revbayes_proj__PhyloCTMCClado__) */
