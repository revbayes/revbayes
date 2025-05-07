#ifndef AbstractPhyloCTMCSiteHomogeneous_H
#define AbstractPhyloCTMCSiteHomogeneous_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "ConstantNode.h"
#include "DiscreteTaxonData.h"
#include "DnaState.h"
#include "MatrixReal.h"
#include "EigenSystem.h"
#include "MemberObject.h"
#include "RbConstants.h"
#include "RbMathLogic.h"
#include "RbSettings.h"
#include "RbVector.h"
#include "RateGenerator.h"
#include "Simplex.h"
#include "TopologyNode.h"
#include "TransitionProbabilityMatrix.h"
#include "Tree.h"
#include "TreeChangeEventListener.h"
#include "TypedDistribution.h"
#include "SiteMixtureModel.h"

#include <memory.h>

#if defined( RB_BEAGLE )
#include "RlUserInterface.h"
#include "BeagleInstance.h"
#include "BeagleUtilities.h"
#include "libhmsbeagle/beagle.h"

#if defined ( RB_USE_EIGEN3 )
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif /* RB_USE_EIGEN3 */

#endif /* RB_BEAGLE */


namespace RevBayesCore {

    /**
     * @brief Homogeneous distribution of character state evolution along a tree class (PhyloCTMC).
     *
     * This file contains the distribution class for a character state evolving along a tree.
     * This abstract base class can be derived for any character evolution model with homogeneous mixture sites. A
     * homogeneous mixture model over sites is a model where all sites are drawn from the same distribution and the
     * specific instance of the per site parameter is integrated over.
     * The per site parameter could be a rate scaler (e.g. the + gamma models) or different rate matrices
     * or anything else.
     *
     * The pruning algorithm is implemented in this base class and calls some few pure virtual methods.
     * The important functions you have to override are:
     * - computeRootLikelihood(size_t root, size_t l, size_t r, size_t m)
     * - computeInternalNodeLikelihoodBranchWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r)
     * - computeTipLikelihood(const TopologyNode &node, size_t nIdx)
     * - getRootFrequencies()
     * - updateTransitionProbabilities()
     *
     *
     * The data are stored for convenience in this class in a matrix (std::vector<std::vector< unsigned > >) 
     * and can be compressed.
     *
     * The partial likelihoods are stored in a c-style array called partial_branch_likelihoods. The dimension are
     * partial_branch_likelihoods[active][node_index][siteRateIndex][siteIndex][charIndex], however, since this is a one-dimensional c-style array,
     * you have to access the partial_branch_likelihoods via
     * partial_branch_likelihoods[active*num_nodes*num_site_mixtures*pattern_block_size*num_states +
     *                    node_index*num_site_mixtures*pattern_block_size*num_states +
     *                    siteRateIndex*pattern_block_size*num_states +
     *                    siteIndex*num_states +
     *                    charIndex]
     * Since this is a bit complex, we have some offset variables for convenience:
     * active_branch_likelihood_offset  =  num_nodes*num_site_mixtures*pattern_block_size*num_states;
     * node_offset                      =  num_site_mixtures*pattern_block_size*num_states;
     * mixture_offset                   =  pattern_block_size*num_states;
     * site_offset                      =  num_states;
     * This gives the more convenient access via
     * partial_branch_likelihoods[active*active_branch_likelihood_offset + node_index*node_offset + siteRateIndex*mixture_offset + siteIndex*site_offset + charIndex]
     *
     * Our implementation of the partial likelihoods means that we can store the partial likelihood of a node,
     * but not for site rates.
     * We also use twice as much memory because we store the partial likelihood along each branch and not only 
     * for each internal node.
     * This gives us a speed improvement during MCMC proposal in the order of a factor 2.
     *
     * The transition probability matrices are stored in a c-style array called partial_branch_likelihoods. The dimension are
     * pmatrices[active][node_index][siteMixtureIndex], however, since this is a one-dimensional c-style array,
     * you have to access the partial_branch_likelihoods via
     * pmatrices[active * num_nodes * num_site_mixtures +
     *                    node_index * num_site_mixtures +
     *                    site_mixture_index]
     * Since this is a bit complex, we have some offset variables for convenience:
     * active_P_matrix_offset         =  num_nodes * num_site_mixtures;
     * node_offset                  =  num_site_mixtures;
     * This gives the more convenient access via
     * pmatrices[active * active_P_matrix_offset + node_index * node_offset + site_mixture_index]
     *
     */
    template<class charType>
    class AbstractPhyloCTMCSiteHomogeneous : public TypedDistribution< AbstractHomologousDiscreteCharacterData >, public MemberObject< RbVector<double> >, public MemberObject < MatrixReal >, public TreeChangeEventListener {

	virtual void                                                       checkInvariants() const;
    public:
        
        enum PARTIAL_LIKELIHOOD_STORING { BRANCH, NODE, BOTH };
        
        // Note, we need the size of the alignment in the constructor to correctly simulate an initial state
        AbstractPhyloCTMCSiteHomogeneous(const TypedDagNode<Tree> *t, size_t nChars, size_t nMix, bool c, size_t nSites, bool amb, bool wd = false, bool internal = false, bool gapmatch = true );
        AbstractPhyloCTMCSiteHomogeneous(const AbstractPhyloCTMCSiteHomogeneous &n);                                                                                    //!< Copy constructor
        virtual                                                            ~AbstractPhyloCTMCSiteHomogeneous(void);                                                     //!< Virtual destructor

        // public member functions
        // pure virtual
        virtual AbstractPhyloCTMCSiteHomogeneous*                           clone(void) const = 0;                                                                      //!< Create an independent clone

        // non-virtual
        void                                                                bootstrap(void);
        virtual double                                                      computeLnProbability(void);
        virtual std::vector<charType>                                       drawAncestralStatesForNode(const TopologyNode &n);
        virtual void                                                        drawJointConditionalAncestralStates(std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates); //!< Simulate ancestral states for each node and each site
        virtual void                                                        drawSiteMixtureAllocations(); //!< For site mixture models (rates and/or matrices), sample the allocation of each site among the mixture categories
        virtual void                                                        drawStochasticCharacterMap(std::vector<std::string>& character_histories, size_t site, bool use_simmap_default=true); //!< Simulate the history of evolution along each branch for each site
        void                                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;     //!< Map the member methods to internal function calls
        void                                                                executeMethod(const std::string &n, const std::vector<const DagNode*> &args, MatrixReal &rv) const;     //!< Map the member methods to internal function calls
        void                                                                fireTreeChangeEvent(const TopologyNode &n, const unsigned& m=0);                                                 //!< The tree has changed and we want to know which part.
        virtual void                                                        recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates, const std::vector<size_t>& sampledSiteRates); //!< Simulate the ancestral states for a given node, conditional on its ancestor's state and the tip data
        virtual bool                                                        recursivelyDrawStochasticCharacterMap(const TopologyNode &node, std::vector<std::string>& character_histories, std::vector<std::vector<charType> >& start_states, std::vector<std::vector<charType> >& end_states, size_t site, bool use_simmap_default); //!< Simulate the history of evolution for a given site on a given branch, conditional on start and end states
        virtual void                                                        redrawValue(void);
        void                                                                reInitialized(void);
        void                                                                setMcmcMode(bool tf);                                                                       //!< Change the likelihood computation to or from MCMC mode.
        void                                                                setValue(AbstractHomologousDiscreteCharacterData *v, bool f=false);                         //!< Set the current value, e.g. attach an observation (clamp)
        virtual void                                                        tipDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates, const std::vector<size_t>& sampledSiteRates);
        void	                                                            updateMarginalNodeLikelihoods(void);
        const TypedDagNode<Tree>*                                           getTree(void);

        void                                                                setClockRate(const TypedDagNode< double > *r);
        void                                                                setClockRate(const TypedDagNode< RbVector< double > > *r);
        void                                                                setPInv(const TypedDagNode< double > *);
        void                                                                setMixtureModel(const TypedDagNode< SiteMixtureModel > *mm);
        void                                                                setRateMatrix(const TypedDagNode< RateGenerator > *rm);
        void                                                                setRateMatrix(const TypedDagNode< RbVector< RateGenerator > > *rm);
        void                                                                setRootFrequencies(const TypedDagNode< Simplex > *f);
        void                                                                setSiteRates(const TypedDagNode< RbVector< double > > *r);
        void                                                                setSiteRatesProbs(const TypedDagNode< Simplex > *rp);
        void                                                                setUseMarginalLikelihoods(bool tf);
        void                                                                setUseSiteMatrices(bool sm, const TypedDagNode< Simplex > *s = NULL);
        void                                                                swap_taxon_name_2_tip_index(std::string tip1, std::string tip2);

        bool                                                                hasSiteRateMixture();
        bool                                                                hasSiteMatrixMixture();
        void                                                                getSampledMixtureComponents(size_t &site_index, size_t &rate_component, size_t &matrix_component );

    protected:

        // helper method for this and derived classes
        void                                                                recursivelyFlagBranchDirty(const TopologyNode& n);
        void                                                                recursivelyFlagNodeDirty(const TopologyNode& n);
        void                                                                flagNodeDirtyPmatrix(size_t node_idx);
        //-- reinitialize beagle in this function
        virtual void                                                        resizeLikelihoodVectors(void);
        virtual void                                                        setActivePIDSpecialized(size_t i, size_t n);                                                          //!< Set the number of processes for this distribution.
        virtual void                                                        updateTransitionProbabilities(size_t node_idx);
        virtual void                                                        updateTransitionProbabilityMatrices(void);
        virtual std::vector<double>                                         getRootFrequencies( size_t mixture = 0 ) const;
        virtual void                                                        getRootFrequencies( std::vector<std::vector<double> >& ) const;
        virtual std::vector<double>                                         getMixtureProbs( void ) const;
        virtual double                                                      getPInv(void) const;


        // Parameter management functions.
        virtual void                                                        swapParameterInternal(const DagNode *oldP, const DagNode *newP);                             //!< Swap a parameter


        // virtual methods that may be overwritten, but then the derived class should call this methods
        virtual void                                                        keepSpecialization(const DagNode* affecter);
        virtual void                                                        restoreSpecialization(const DagNode *restorer);
        virtual void                                                        touchSpecialization(const DagNode *toucher, bool touchAll);

        // pure virtual methods
        virtual void                                                        computeInternalNodeLikelihoodBranchWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r) = 0;
        virtual void                                                        computeTipLikelihood(const TopologyNode &node, size_t nIdx) = 0;
        virtual void                                                        computeRootLikelihood( size_t root, size_t left, size_t right) = 0;
        virtual void                                                        computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle) = 0;

        virtual void                                                        computeInternalNodeLikelihoodNodeWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r) { throw RbException("Missing implementation of computeInternalNodeLikelihoodNodeWise."); }
        virtual void                                                        computeRootLikelihoodNode( size_t root, size_t left, size_t right) { throw RbException("Missing implementation of computeRootLikelihoodNode."); }
        virtual void                                                        computeRootLikelihoodNode( size_t root, size_t left, size_t right, size_t middle) { throw RbException("Missing implementation of computeRootLikelihoodNode."); }

        virtual void                                                        computeInternalNodeLikelihoodBranchNodeWise(const TopologyNode &n, size_t nIdx, size_t l, size_t r) { throw RbException("Missing implementation of computeInternalNodeLikelihoodBranchNodeWise."); }
        virtual void                                                        computeRootLikelihoodBranchNode( size_t root, size_t left, size_t right) { throw RbException("Missing implementation of computeRootLikelihoodBranchNode."); }
        virtual void                                                        computeRootLikelihoodBranchNode( size_t root, size_t left, size_t right, size_t middle) { throw RbException("Missing implementation of computeRootLikelihoodBranchNode."); }

        // virtual methods that you may want to overwrite
        virtual void                                                        compress(void);
        virtual void                                                        computeMarginalNodeLikelihood(size_t node_idx, size_t parentIdx);
        virtual void                                                        computeMarginalRootLikelihood();
        virtual void                                                        computeRootLikelihoods( std::vector< double > &rv ) const;
        virtual void                                                        computeRootLikelihoodsPerSiteMixture( MatrixReal &rv ) const;
        virtual void                                                        computeRootLikelihoodsPerSiteRate( MatrixReal &rv ) const;
        virtual void                                                        fillTipLikelihoods(void);
        virtual double                                                      sumRootLikelihood( void );
        virtual std::vector< std::vector< double > >*                       sumMarginalLikelihoods(size_t node_index);
        virtual std::vector<size_t>                                         getIncludedSiteIndices();

        virtual void                                                        scale(size_t i);
        virtual void                                                        scale(size_t i, size_t l, size_t r);
        virtual void                                                        scale(size_t i, size_t l, size_t r, size_t m);
        
        // members
        double                                                              lnProb = 0;
        std::optional<double>                                               storedLnProb;
        PARTIAL_LIKELIHOOD_STORING                                          partial_likelihood_storing_approach;
        size_t                                                              num_matrices = 1;
        size_t                                                              num_nodes;
        size_t                                                              num_sites;
        const size_t                                                        num_states;
        size_t                                                              num_site_rates;
        size_t                                                              num_site_mixtures;
        size_t                                                              num_tips;
        const TypedDagNode<Tree>*                                           tau;
        
        std::vector<TransitionProbabilityMatrix>                            transition_prob_matrices;
        std::vector<TransitionProbabilityMatrix>                            pmatrices;
        std::vector<size_t>                                                 active_pmatrices;
        std::vector<bool>                                                   pmat_changed_nodes;
        mutable std::vector<bool>                                           pmat_dirty_branches;
        std::optional<std::vector<bool>>                                    prev_pmat_dirty_branches;
        
        // offsets for nodes
        size_t                                                              active_P_matrix_offset;
        size_t                                                              pmat_node_offset;

        // the likelihoods
        mutable double*                                                     partial_branch_likelihoods = nullptr;
        mutable double*                                                     partial_node_likelihoods = nullptr;
        mutable double*                                                     tip_likelihoods = nullptr;
        std::vector<size_t>                                                 active_branch_likelihood;
        std::vector<size_t>                                                 active_node_likelihood;
        double*                                                             marginalLikelihoods = nullptr;
        std::optional<std::vector<bool>>                                    prev_dirty_nodes;
        std::optional<std::vector<bool>>                                    prev_dirty_branches;

        std::vector< std::vector< std::vector<double> > >                   per_node_site_log_scaling_factors;
        std::vector< std::vector< std::vector<std::vector<double> > > >     per_node_site_mixture_log_scaling_factors;

        // the data
        std::vector<std::vector<RbBitSet> >                                 ambiguous_char_matrix;
        std::vector<std::vector<std::uint64_t> >                            char_matrix;
        std::vector<std::vector<bool> >                                     gap_matrix;
        std::vector<size_t>                                                 pattern_counts;
        std::vector<bool>                                                   site_invariant;
        std::vector<std::vector<size_t> >                                   invariant_site_index;
        size_t                                                              num_patterns;
        bool                                                                compressed;
        std::vector<size_t>                                                 site_pattern;    // an array that keeps track of which pattern is used for each site
        std::map<std::string,size_t>                                        taxon_name_2_tip_index_map;

        // flags for likelihood recomputation
        bool                                                                touched = false;
        std::vector<bool>                                                   changed_branches;
        std::vector<bool>                                                   changed_nodes;
        mutable std::vector<bool>                                           dirty_branches;
        mutable std::vector<bool>                                           dirty_nodes;

#if defined ( RB_BEAGLE )
        mutable std::vector<int>                                            active_eigen_system;
        mutable std::vector<bool>                                           touched_eigen_system;
#endif /* RB_BEAGLE */

        // offsets for nodes
        size_t                                                              active_branch_likelihood_offset;
        size_t                                                              active_node_likelihood_offset;
        size_t                                                              node_offset;
        size_t                                                              tip_offset;
        size_t                                                              mixture_offset;
        size_t                                                              site_offset;

        // flags
        bool                                                                using_ambiguous_characters;
        bool                                                                treatUnknownAsGap = true;
        bool                                                                treatAmbiguousAsGaps = false;

        bool                                                                using_weighted_characters;

        bool                                                                useMarginalLikelihoods = false;
        mutable bool                                                        in_mcmc_mode = false;

        // members
        const TypedDagNode< double >*                                       homogeneous_clock_rate = nullptr;
        const TypedDagNode< RbVector< double > >*                           heterogeneous_clock_rates = nullptr;
        const TypedDagNode< SiteMixtureModel >*                     mixture_model = nullptr;
        const TypedDagNode< RateGenerator >*                                homogeneous_rate_matrix = nullptr;
        const TypedDagNode< RbVector< RateGenerator > >*                    heterogeneous_rate_matrices = nullptr;
        const TypedDagNode< Simplex >*                                      root_frequencies = nullptr;
        const TypedDagNode< RbVector< double > >*                           site_rates = nullptr;
        const TypedDagNode< Simplex >*                                      site_matrix_probs = nullptr;
        const TypedDagNode< Simplex >*                                      site_rates_probs = nullptr;
        const TypedDagNode< double >*                                       p_inv = nullptr;


        // flags specifying which model variants we use
        bool                                                                branch_heterogeneous_clock_rates = false;
        bool                                                                branch_heterogeneous_substitution_matrices = true;
        bool                                                                rate_variation_across_sites = false;

        // MPI variables
        size_t                                                              pattern_block_start;
        size_t                                                              pattern_block_end;
        size_t                                                              pattern_block_size;

        bool                                                                store_internal_nodes;
        bool                                                                gap_match_clamped;

        charType                                                            template_state;                                 //!< Template state used for ancestral state estimation. This makes sure that the state labels are preserved.

        // containers for ancestral state/stochastic mapping functions
        bool                                                                has_ancestral_states = false;
        std::vector<size_t>                                                 sampled_site_mixtures;
        size_t                                                              sampled_site_rate_component = 0;
        size_t                                                              sampled_site_matrix_component = 0;

#if defined( RB_BEAGLE )
        bool                                                                b_initialized = false;
        BeagleInstance*                                                     beagle_instance;

        std::vector<double>                                                 b_inPatternWeights;
        std::vector<double>                                                 b_inCategoryRates;
        std::vector<double>                                                 b_inCategoryWeights;
        size_t                                                              num_mixtures;
                
        //-- Initialize a BEAGLE instance per model
        void                                                                initializeBeagleInstances      ( void );
        //-- Free all BEAGLE instances
        void                                                                freeBeagleInstances            ( void );

        //-- Set the tree tip states and base frequencies for the BEAGLE instance.
        void                                                                initializeBeagleTips           ( void );

        double                                                              calculateBranchLength          ( const TopologyNode &node
                                                                                                           , size_t node_index
                                                                                                           );
        
        //-- BEAGLE updaters
        void                                                                updateBeagleEigensystems       ( void );
        void                                                                updateBeagleSiteRates          ( void );
#endif /* RB_BEAGLE */
        
    private:

        // private methods
        void                                                                fillLikelihoodVector(const TopologyNode &n, size_t nIdx);
        void                                                                recursiveMarginalLikelihoodComputation(size_t nIdx);
        virtual void                                                        simulate(const TopologyNode& node, std::vector< DiscreteTaxonData< charType > > &t, const std::vector<bool> &inv, const std::vector<size_t> &perSiteRates);
        virtual void                                                        updateTransitionProbabilityMatrix(size_t node_idx);
        
    };

}


#include "DiscreteCharacterState.h"
#include "DistributionExponential.h"
#include "HomologousDiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RateMatrix_JC.h"
#include "StochasticNode.h"

#include <cmath>

#ifdef RB_MPI
#include <mpi.h>
#endif

#if defined( RB_BEAGLE )
#include "BeagleInstance.h"
#endif /* RB_BEAGLE */


template<class charType>
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::AbstractPhyloCTMCSiteHomogeneous(const TypedDagNode<Tree> *t, size_t nChars, size_t nMix, bool c, size_t nSites, bool amb, bool internal, bool gapmatch, bool wd) :
TypedDistribution< AbstractHomologousDiscreteCharacterData >(  NULL ),
num_nodes( t->getValue().getNumberOfNodes() ),
num_sites( nSites ),
num_states( nChars ),
num_site_rates( nMix ),
num_site_mixtures( nMix ),
num_tips( t->getValue().getNumberOfTips() ),
tau( t ),
transition_prob_matrices( std::vector<TransitionProbabilityMatrix>(num_site_mixtures, TransitionProbabilityMatrix(num_states) ) ),
partial_branch_likelihoods( NULL ),
partial_node_likelihoods( NULL ),
tip_likelihoods( NULL ),
active_branch_likelihood( std::vector<size_t>(num_nodes, 0) ),
active_node_likelihood( std::vector<size_t>(num_nodes, 0) ),
per_node_site_mixture_log_scaling_factors( std::vector<std::vector< std::vector< std::vector<double> > > >(2, std::vector<std::vector< std::vector<double> > >(num_nodes, std::vector< std::vector<double> >( num_site_mixtures, std::vector<double>(num_sites, 0.0) ) ) ) ),
per_node_site_log_scaling_factors( std::vector<std::vector< std::vector< double > > >(2, std::vector<std::vector< double > >(num_nodes, std::vector< double >( num_sites, 0.0 ) ) ) ),
site_invariant( num_sites, false ),
invariant_site_index( num_sites ),
num_patterns( num_sites ),
compressed( c ),
site_pattern( std::vector<size_t>(num_sites, 0) ),
taxon_name_2_tip_index_map(),
changed_branches( std::vector<bool>(num_nodes, false) ),
changed_nodes( std::vector<bool>(num_nodes, false) ),
dirty_branches( std::vector<bool>(num_nodes, true) ),
dirty_nodes( std::vector<bool>(num_nodes, true) ),
#if defined( RB_BEAGLE )
active_eigen_system( std::vector<int>(1, 0) ),
touched_eigen_system( std::vector<bool>(1, false) ),
beagle_instance( NULL ),
#endif /* RB_BEAGLE */
using_ambiguous_characters( amb ),
using_weighted_characters( wd ),
pattern_block_start( 0 ),
pattern_block_end( num_patterns ),
pattern_block_size( num_patterns ),
store_internal_nodes( internal ),
gap_match_clamped( gapmatch )
{

    tau->getValue().getTreeChangeEventHandler().addListener( this );

    active_branch_likelihood_offset =  num_nodes*num_site_mixtures*pattern_block_size*num_states;
    active_node_likelihood_offset   =  num_tips*num_site_mixtures*pattern_block_size*num_states;
    node_offset                     =  num_site_mixtures*pattern_block_size*num_states;
    tip_offset                      =  pattern_block_size*num_states;
    mixture_offset                  =  pattern_block_size*num_states;
    site_offset                     =  num_states;
    
    active_P_matrix_offset          =  num_nodes * num_site_mixtures;
    pmat_node_offset                =  num_site_mixtures;
    active_pmatrices                =  std::vector<size_t>(num_nodes, 0);
    pmat_changed_nodes              =  std::vector<bool>(num_nodes, false);
    pmat_dirty_branches             =  std::vector<bool>(num_nodes, true);
    pmatrices                       =  std::vector<TransitionProbabilityMatrix>(active_P_matrix_offset * 2, TransitionProbabilityMatrix(num_states));

    // add the parameters to our set (in the base class)
    // in that way other class can easily access the set of our parameters
    // this will also ensure that the parameters are not getting deleted before we do
    this->addParameter( tau );
    this->addParameter( homogeneous_clock_rate );
    this->addParameter( heterogeneous_clock_rates );
    this->addParameter( mixture_model );
    this->addParameter( homogeneous_rate_matrix );
    this->addParameter( heterogeneous_rate_matrices );
    this->addParameter( root_frequencies );
    this->addParameter( site_rates );
    this->addParameter( site_matrix_probs );
    this->addParameter( site_rates_probs );
    this->addParameter( p_inv );

    // initially we use only a single processor until someone else tells us otherwise
    this->setActivePID( this->pid, 1 );
}


template<class charType>
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::AbstractPhyloCTMCSiteHomogeneous(const AbstractPhyloCTMCSiteHomogeneous &n) :
TypedDistribution< AbstractHomologousDiscreteCharacterData >( n ),
lnProb( n.lnProb ),
storedLnProb( n.storedLnProb ),
num_nodes( n.num_nodes ),
num_sites( n.num_sites ),
num_states( n.num_states ),
num_site_rates( n.num_site_rates ),
num_site_mixtures( n.num_site_mixtures ),
num_matrices( n.num_matrices ),
num_tips( n.num_tips ),
tau( n.tau ),
transition_prob_matrices( n.transition_prob_matrices ),
partial_branch_likelihoods( NULL ),
partial_node_likelihoods( NULL ),
tip_likelihoods( NULL ),
active_branch_likelihood( n.active_branch_likelihood ),
active_node_likelihood( n.active_node_likelihood ),
marginalLikelihoods( NULL ),
per_node_site_mixture_log_scaling_factors( n.per_node_site_mixture_log_scaling_factors ),
per_node_site_log_scaling_factors( n.per_node_site_log_scaling_factors ),
ambiguous_char_matrix( n.ambiguous_char_matrix ),
char_matrix( n.char_matrix ),
gap_matrix( n.gap_matrix ),
pattern_counts( n.pattern_counts ),
site_invariant( n.site_invariant ),
invariant_site_index( n.invariant_site_index ),
num_patterns( n.num_patterns ),
compressed( n.compressed ),
site_pattern( n.site_pattern ),
taxon_name_2_tip_index_map( n.taxon_name_2_tip_index_map ),
touched( false ),
changed_branches( n.changed_branches ),
changed_nodes( n.changed_nodes ),
dirty_branches( n.dirty_branches ),
dirty_nodes( n.dirty_nodes ),
#if defined( RB_BEAGLE )
active_eigen_system( n.active_eigen_system ),
touched_eigen_system( n.touched_eigen_system ),
beagle_instance( NULL ),
#endif /* RB_BEAGLE */
prev_dirty_nodes( n.prev_dirty_nodes ),
prev_dirty_branches( n.prev_dirty_branches ),
using_ambiguous_characters( n.using_ambiguous_characters ),
treatUnknownAsGap( n.treatUnknownAsGap ),
treatAmbiguousAsGaps( n.treatAmbiguousAsGaps ),
using_weighted_characters( n.using_weighted_characters ),
useMarginalLikelihoods( n.useMarginalLikelihoods ),
in_mcmc_mode( n.in_mcmc_mode ),
pattern_block_start( n.pattern_block_start ),
pattern_block_end( n.pattern_block_end ),
pattern_block_size( n.pattern_block_size ),
store_internal_nodes( n.store_internal_nodes ),
gap_match_clamped( n.gap_match_clamped ),
template_state( n.template_state ),
has_ancestral_states( n.has_ancestral_states ),
sampled_site_rate_component( n.sampled_site_rate_component ),
sampled_site_matrix_component( n.sampled_site_matrix_component )
{

    // initialize with default parameters
    homogeneous_clock_rate          = n.homogeneous_clock_rate;
    heterogeneous_clock_rates       = n.heterogeneous_clock_rates;
    mixture_model                = n.mixture_model;
    homogeneous_rate_matrix         = n.homogeneous_rate_matrix;
    heterogeneous_rate_matrices     = n.heterogeneous_rate_matrices;
    root_frequencies                = n.root_frequencies;
    site_rates                      = n.site_rates;
    site_matrix_probs               = n.site_matrix_probs;
    site_rates_probs                = n.site_rates_probs;
    p_inv                           = n.p_inv;

    active_branch_likelihood_offset =  n.active_branch_likelihood_offset;
    active_node_likelihood_offset   =  n.active_node_likelihood_offset;
    node_offset                     =  n.node_offset;
    tip_offset                      =  n.tip_offset;
    mixture_offset                  =  n.mixture_offset;
    site_offset                     =  n.site_offset;
    
    active_P_matrix_offset          =  n.active_P_matrix_offset;
    pmat_node_offset                =  n.pmat_node_offset;
    active_pmatrices                =  n.active_pmatrices;
    pmat_changed_nodes              =  n.pmat_changed_nodes;
    pmat_dirty_branches             =  n.pmat_dirty_branches;
    prev_pmat_dirty_branches           =  n.prev_pmat_dirty_branches;
    pmatrices                       =  n.pmatrices;

    // flags specifying which model variants we use
    branch_heterogeneous_clock_rates               = n.branch_heterogeneous_clock_rates;
    branch_heterogeneous_substitution_matrices     = n.branch_heterogeneous_substitution_matrices;
    rate_variation_across_sites                    = n.rate_variation_across_sites;

    tau->getValue().getTreeChangeEventHandler().addListener( this );

    // copy the partial likelihoods if necessary
    if ( in_mcmc_mode == true )
    {
#if defined( RB_BEAGLE )
        //-- If we are in MCMC mode, we need to make a new BEAGLE instance for our clone.
        charType tmp_char;
        if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
        {
            this->initializeBeagleInstances();
        }
        else
        {
            if ( n.partial_branch_likelihoods != NULL )
            {
                partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
                memcpy(partial_branch_likelihoods, n.partial_branch_likelihoods, 2*active_branch_likelihood_offset*sizeof(double));
            }
            else
            {
                partial_branch_likelihoods = NULL;
            }
            
            if ( n.partial_node_likelihoods != NULL )
            {
                partial_node_likelihoods = new double[2*active_node_likelihood_offset];
                memcpy(partial_node_likelihoods, n.partial_node_likelihoods, 2*active_node_likelihood_offset*sizeof(double));
            }
            else
            {
                partial_node_likelihoods = NULL;
            }
            
            if ( n.tip_likelihoods != NULL )
            {
                tip_likelihoods = new double[num_tips*pattern_block_size*num_states];
                memcpy(tip_likelihoods, n.tip_likelihoods, num_tips*pattern_block_size*num_states*sizeof(double));
            }
            else
            {
                tip_likelihoods = NULL;
            }
        }
#else
        if ( n.partial_branch_likelihoods != NULL )
        {
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
            memcpy(partial_branch_likelihoods, n.partial_branch_likelihoods, 2*active_branch_likelihood_offset*sizeof(double));
        }
        else
        {
            partial_branch_likelihoods = NULL;
        }
        
        if ( n.partial_node_likelihoods != NULL )
        {
            partial_node_likelihoods = new double[2*active_node_likelihood_offset];
            memcpy(partial_node_likelihoods, n.partial_node_likelihoods, 2*active_node_likelihood_offset*sizeof(double));
        }
        else
        {
            partial_node_likelihoods = NULL;
        }
        
        if ( n.tip_likelihoods != NULL )
        {
            tip_likelihoods = new double[num_tips*pattern_block_size*num_states];
            memcpy(tip_likelihoods, n.tip_likelihoods, num_tips*pattern_block_size*num_states*sizeof(double));
        }
        else
        {
            tip_likelihoods = NULL;
        }
#endif /* RB_BEAGLE */
    }

    // copy the marginal likelihoods if necessary
    if ( useMarginalLikelihoods == true && RbSettings::userSettings().getUseBeagle() != true )
    {
        marginalLikelihoods = new double[active_branch_likelihood_offset];
        memcpy(marginalLikelihoods, n.marginalLikelihoods, active_branch_likelihood_offset*sizeof(double));
    }
}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::checkInvariants( void ) const
{
    // 1. Invariant: if a node is dirty then its parent is also dirty
    auto tree_nodes = tau->getValue().getNodes();
    for(auto node: tree_nodes)
    {
	if (dirty_nodes[node->getIndex()] and not node->isRoot())
	    assert(dirty_nodes[node->getParent().getIndex()]);
    }

    // 2. Invariant: if the P-matrix for a node is dirty, then the conditional likelihoods should be dirty
    for(auto node: tree_nodes)
    {
	int index = node->getIndex();
	if (pmat_dirty_branches[index])
	    assert(dirty_branches[index]);
    }
}

/**
 * Destructor. Because we added ourselves as a reference to tau when we added a listener to its
 * TreeChangeEventHandler, we need to remove ourselves as a reference and possibly delete tau
 * when we die. All other parameters are handled by others.
 */
template<class charType>
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::~AbstractPhyloCTMCSiteHomogeneous( void )
{
    // We don't delete the params, because they might be used somewhere else too. The model needs to do that!

    // remove myself from the tree listeners
    if ( tau != NULL )
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );
    }
    
#ifdef RB_BEAGLE
    charType tmp_char;
    if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
    {
        this->freeBeagleInstances();
    }
    else
    {
        // If BEAGLE is not used, we still need to free the partial likelihoods
        delete[] partial_branch_likelihoods;
        delete[] partial_node_likelihoods;
        delete[] tip_likelihoods;
        delete[] marginalLikelihoods;
    }
#else
    // free the partial likelihoods
    delete[] partial_branch_likelihoods;
    delete[] partial_node_likelihoods;
    delete[] tip_likelihoods;
    delete[] marginalLikelihoods;
#endif
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::bootstrap( void )
{

    // first we re-compress the data
    compress();

    RandomNumberGenerator *rng = GLOBAL_RNG;

    std::vector<size_t> bootstrapped_pattern_counts = std::vector<size_t>(num_patterns,0);

    for (size_t i = 0; i<num_sites; ++i)
    {
        double u = rng->uniform01() * num_sites;
        size_t pattern_index = 0;
        while ( u > double(pattern_counts[pattern_index]) )
        {
            u -= double(pattern_counts[pattern_index]);
            ++pattern_index;
        }

        ++bootstrapped_pattern_counts[pattern_index];

    }

    pattern_counts = bootstrapped_pattern_counts;

}

namespace RevBayesCore
{
inline void mark_ambiguous_and_missing_as_gap(AbstractHomologousDiscreteCharacterData& data, const vector<size_t>& site_indices, std::vector<TopologyNode*> nodes)
{
    for (auto& node: nodes)
    {
        if ( node->isTip() )
        {
            AbstractDiscreteTaxonData& taxon_data = data.getTaxonData( node->getName() );
            for (auto site_index: site_indices)
            {
                DiscreteCharacterState &c = taxon_data.getCharacter(site_index);

                if ( c.isAmbiguous() or  c.isMissingState() )
                {
                    c.setGapState( true );
                }
            }
        }
    }
}

inline void mark_unknown_as_gap(AbstractHomologousDiscreteCharacterData& data, const vector<size_t>& site_indices, std::vector<TopologyNode*> nodes)
{
    for (auto& node: nodes)
    {
        if ( node->isTip() )
        {
            AbstractDiscreteTaxonData& taxon_data = data.getTaxonData( node->getName() );
            for (auto site_index: site_indices)
            {
                DiscreteCharacterState &c = taxon_data.getCharacter(site_index);

                if ( c.getNumberOfStates() == c.getNumberObservedStates() or c.isMissingState())
                {
                    c.setGapState( true );
                }
            }
        }
    }
}

inline bool has_ambiguous_nongap_characters(AbstractHomologousDiscreteCharacterData& data, const vector<size_t>& site_indices, std::vector<TopologyNode*> nodes)
{
    for (auto& node: nodes)
    {
        if ( node->isTip() )
        {
            AbstractDiscreteTaxonData& taxon_data = data.getTaxonData( node->getName() );
            for (auto site_index: site_indices)
            {
                DiscreteCharacterState &c = taxon_data.getCharacter(site_index);

                if ( not c.isGapState() and (c.isAmbiguous() or c.isMissingState()) )
                    return true;
            }
        }
    }

    return false;
}

inline bool has_weighted_characters(AbstractHomologousDiscreteCharacterData& data, const vector<size_t>& site_indices, std::vector<TopologyNode*> nodes)
{
    for (auto& node: nodes)
    {
        if ( node->isTip() )
        {
            AbstractDiscreteTaxonData& taxon_data = data.getTaxonData( node->getName() );
            for (auto site_index: site_indices)
            {
                DiscreteCharacterState &c = taxon_data.getCharacter(site_index);

                if ( c.isWeighted() ) return true;
            }
        }
    }

    return false;
}

}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::compress( void )
{
    // only if the value has been set
    if ( this->value == NULL )
    {
        return;
    }

    ambiguous_char_matrix.clear();
    char_matrix.clear();
    gap_matrix.clear();
    pattern_counts.clear();
    num_patterns = 0;

    // resize the matrices
    size_t tips = tau->getValue().getNumberOfTips();
    ambiguous_char_matrix.resize(tips);
    char_matrix.resize(tips);
    gap_matrix.resize(tips);

    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> site_indices = getIncludedSiteIndices();

    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();

    if (treatAmbiguousAsGaps)
    {
        mark_ambiguous_and_missing_as_gap(*value, site_indices, nodes);
    }
    
    if (treatUnknownAsGap)
    {
        mark_unknown_as_gap(*value, site_indices, nodes);
    }
    
    // set the global variable if we use ambiguous characters (besides gaps)
    using_ambiguous_characters = has_ambiguous_nongap_characters(*value, site_indices, nodes);

    // set the global variable if we use weighted characters
    using_weighted_characters = has_weighted_characters(*value, site_indices, nodes);

    std::vector<bool> unique(num_sites, true);
    std::vector<size_t> indexOfSitePattern;

    // compress the character matrix if we're asked to
    if ( compressed == true )
    {
        // find the unique site patterns and compute their respective frequencies
        std::map<std::string,size_t> patterns;
        for (size_t site = 0; site < num_sites; ++site)
        {
            bool all_missing = true;
            // create the site pattern
            std::string pattern = "";
            for (auto& node: nodes)
            {
                if ( node->isTip() )
                {
                    AbstractDiscreteTaxonData& taxon = value->getTaxonData( node->getName() );
                    CharacterState &c = taxon.getCharacter(site_indices[site]);
                    all_missing &= (c.isMissingState() || c.isGapState());
                    pattern += c.getStringValue();
                }
            }
            
            // only add this pattern if not all are missing
            if ( all_missing == false )
            {
                // check if we have already seen this site pattern
                std::map<std::string, size_t>::const_iterator index = patterns.find( pattern );
                if ( index != patterns.end() )
                {
                    // we have already seen this pattern
                    // increase the frequency counter
                    pattern_counts[ index->second ]++;

                    // obviously this site isn't unique nor the first encounter
                    unique[site] = false;

                    // remember which pattern this site uses
                    site_pattern[site] = index->second;
                }
                else
                {
                    // create a new pattern frequency counter for this pattern
                    pattern_counts.push_back(1);
                    
                    // insert this pattern with the corresponding index in the map
                    patterns.insert( std::pair<std::string,size_t>(pattern,num_patterns) );
                    
                    // remember which pattern this site uses
                    site_pattern[site] = num_patterns;
                    
                    // increase the pattern counter
                    num_patterns++;
                    
                    // add the index of the site to our pattern-index vector
                    indexOfSitePattern.push_back( site );
                    
                    // flag that this site is unique (or the first occurence of this pattern)
                    unique[site] = true;
                }
            }
        }
    }
    else
    {
        // we do not compress
        num_patterns = num_sites;
        pattern_counts     = std::vector<size_t>(num_sites,1);
        indexOfSitePattern = std::vector<size_t>(num_sites,1);
        for (size_t i = 0; i < this->num_sites; i++)
        {
            indexOfSitePattern[i] = i;
        }
    }

    // compute which block of the data this process needs to compute
    pattern_block_start = size_t(floor( (double(pid-active_PID)   / num_processes ) * num_patterns) );
    pattern_block_end   = size_t(floor( (double(pid+1-active_PID) / num_processes ) * num_patterns) );
    pattern_block_size  = pattern_block_end - pattern_block_start;


    std::vector<size_t> process_pattern_counts = std::vector<size_t>(pattern_block_size,0);
    taxon_name_2_tip_index_map.clear();
    // allocate and fill the cells of the matrices
    
    for (auto& the_node: nodes)
    {
        if ( the_node->isTip() )
        {
            size_t node_index = the_node->getIndex();
            taxon_name_2_tip_index_map.insert( std::pair<std::string,size_t>(the_node->getName(), node_index) );
            AbstractDiscreteTaxonData& taxon = value->getTaxonData( the_node->getName() );

            // resize the column
            if ( using_ambiguous_characters == true )
            {
                ambiguous_char_matrix[node_index].resize(pattern_block_size);
            }
            else
            {
                char_matrix[node_index].resize(pattern_block_size);
            }
            gap_matrix[node_index].resize(pattern_block_size);
            for (size_t patternIndex = 0; patternIndex < pattern_block_size; ++patternIndex)
            {
                // set the counts for this patter
                process_pattern_counts[patternIndex] = pattern_counts[patternIndex+pattern_block_start];

                charType &c = static_cast<charType &>( taxon.getCharacter(site_indices[indexOfSitePattern[patternIndex+pattern_block_start]]) );
                gap_matrix[node_index][patternIndex] = c.isGapState();

                if ( using_ambiguous_characters == true )
                {
                    // we use the actual state
                    ambiguous_char_matrix[node_index][patternIndex] = c.getState();
                }
                else if ( c.isGapState() == false )
                {
                    // we use the index of the state
                    char_matrix[node_index][patternIndex] = c.getStateIndex();
                    if ( c.getStateIndex() >= this->num_states )
                    {
                        throw RbException("Problem with state index in PhyloCTMC!");
                    }
                    
                }
                else
                {
                    // just to be safe
                    char_matrix[node_index][patternIndex] = -1;
                }

            }

        }

    }

    bool allow_ambiguous_as_invariant = true;

    // now copy back the pattern count vector
    pattern_counts = process_pattern_counts;

    // reset the vector if a site is invariant
    site_invariant.resize( pattern_block_size );
    invariant_site_index.clear();
    invariant_site_index.resize( pattern_block_size );
    size_t length = char_matrix.size();
        
    for (size_t i=0; i<pattern_block_size; ++i)
    {
        bool inv = true;
        size_t taxon_index = 0;

        while ( taxon_index<(length-1) && gap_matrix[taxon_index][i] == true  )
        {
            ++taxon_index;
        }

        if ( using_ambiguous_characters == true )
        {
            RbBitSet val = ambiguous_char_matrix[taxon_index][i];

            for (; taxon_index<length; ++taxon_index)
            {
                if ( gap_matrix[taxon_index][i] == false )
                {
                    val &= ambiguous_char_matrix[taxon_index][i];
                }

                if (   ( allow_ambiguous_as_invariant == true  &&  val.count() == 0 && gap_matrix[taxon_index][i] == false)
                    || ( allow_ambiguous_as_invariant == false && (val.count() == 0 || gap_matrix[taxon_index][i] == true ) ) )
                {
                    inv = false;
                    break;
                }
            }

            for ( size_t c = 0; c < this->num_states; c++ )
            {
                if ( val.test(c) )
                {
                    if ( c < 0 || c >= this->num_states )
                    {
                        throw RbException() << "Possible bug: Invar sites with ambiguous chars at index " << c << " out of bounds!";
                    }
                    invariant_site_index[i].push_back(c);
                }
            }
        }
        else
        {
            std::uint64_t c = char_matrix[taxon_index][i];
            
            if ( c < 0 || c >= this->num_states )
            {
                throw RbException() << "Possible bug: Invar sites with ambiguous chars at site " << i << " out of bounds! Site was " << (gap_matrix[taxon_index][i] ? "Gap" : "No Gap");
            }
            invariant_site_index[i].push_back(c);

            for (; taxon_index<length; ++taxon_index)
            {
                if (   ( allow_ambiguous_as_invariant == true  &&  c != char_matrix[taxon_index][i] && gap_matrix[taxon_index][i] == false)
                    || ( allow_ambiguous_as_invariant == false && (c != char_matrix[taxon_index][i] || gap_matrix[taxon_index][i] == true ) ) )
                {
                    inv = false;
                    break;
                }
            }
        }

        site_invariant[i] = inv;
        
    }

    // finally we resize the partial likelihood vectors to the new pattern counts
    resizeLikelihoodVectors();
    fillTipLikelihoods();
}


template<class charType>
double RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeLnProbability( void )
{
    const std::string& pls_approach = RbSettings::userSettings().getPartialLikelihoodStoring();
    if ( pls_approach == "branch" )
    {
        partial_likelihood_storing_approach = PARTIAL_LIKELIHOOD_STORING::BRANCH;
    }
    else if ( pls_approach == "node" )
    {
        partial_likelihood_storing_approach = PARTIAL_LIKELIHOOD_STORING::NODE;
    }
    else if ( pls_approach == "both" )
    {
        partial_likelihood_storing_approach = PARTIAL_LIKELIHOOD_STORING::BOTH;
    }
    
    // @todo: #thread
    // This part should be done on several threads if possible
    // That means we should probabily call this function as a job,
    // where a job is defined as computing the lnProbability for a subset of the data (block)
    // Sebastian: this call is very slow; a lot of work happens in nextCycle()
    
    // we need to check here if we still are listining to this tree for change events
    // the tree could have been replaced without telling us
    if ( tau->getValue().getTreeChangeEventHandler().isListening( this ) == false )
    {
        tau->getValue().getTreeChangeEventHandler().addListener( this );
        dirty_branches = std::vector<bool>(num_nodes, true);
        dirty_nodes    = std::vector<bool>(num_nodes, true);
        pmat_dirty_branches = std::vector<bool>(num_nodes, true);
    }

    checkInvariants();

    // update transition probability matrices
    this->updateTransitionProbabilityMatrices();
    // if we are not in MCMC mode, then we need to (temporarily) allocate memory
    if ( in_mcmc_mode == false )
    {
        if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BRANCH || partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BOTH )
        {
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
        }
        if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::NODE || partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BOTH )
        {
            partial_node_likelihoods = new double[2*active_node_likelihood_offset];
        }
        
#if defined( RB_BEAGLE )
        // We need to set up BEAGLE here if we are not in mcmc mode...
        charType tmp_char;
        if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
        {
            //-- If there are already BEAGLE instances, delete them.
            this->freeBeagleInstances();
            //-- Initialize fresh BEAGLE instances.
            this->initializeBeagleInstances();
        }
#endif /* RB_BEAGLE */
    }

    // compute the ln probability by recursively calling the probability calculation for each node
    const TopologyNode &root = tau->getValue().getRoot();

    // we start with the root and then traverse down the tree
    size_t root_index = root.getIndex();
    
    if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BRANCH )
    {
        // start by filling the likelihood vector for the children of the root
        if ( dirty_branches[root_index] ) // rooted trees have two children for the root
        {
            // start by filling the likelihood vector for the children of the root
            if ( root.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
            {
                const TopologyNode &left = root.getChild(0);
                size_t left_index = left.getIndex();
                const TopologyNode &right = root.getChild(1);
                size_t right_index = right.getIndex();

                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );

                computeRootLikelihood( root_index, left_index, right_index );
                
                //-- We only need to scale if we are not using BEAGLE
#if !defined ( RB_BEAGLE )
                // now call scaling
                scale(root_index, left_index, right_index);
#endif
            }
            else if ( root.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
            {
                const TopologyNode &left = root.getChild(0);
                size_t left_index = left.getIndex();
                const TopologyNode &right = root.getChild(1);
                size_t right_index = right.getIndex();
                const TopologyNode &middle = root.getChild(2);
                size_t middleIndex = middle.getIndex();

                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );
                fillLikelihoodVector( middle, middleIndex );

                computeRootLikelihood( root_index, left_index, right_index, middleIndex );
                
                //-- We only need to scale if we are not using BEAGLE
#if !defined ( RB_BEAGLE )
                // now call scaling
                scale(root_index, left_index, right_index, middleIndex);
#endif
            }
            else
            {
                throw RbException("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
            }

        }

        // sum the partials up
        this->lnProb = sumRootLikelihood();

    }
    else if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::NODE )
    {
        // only necessary if the root is actually dirty
        if ( dirty_nodes[root_index] == true )
        {
            // start by filling the likelihood vector for the children of the root
            if ( root.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
            {
                const TopologyNode &left = root.getChild(0);
                size_t left_index = left.getIndex();
                const TopologyNode &right = root.getChild(1);
                size_t right_index = right.getIndex();

                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );

                dirty_nodes[root_index] = false;
                computeRootLikelihoodNode( root_index, left_index, right_index );
            
                // now call scaling
                scale(root_index, left_index, right_index);

            }
            else if ( root.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
            {
                const TopologyNode &left = root.getChild(0);
                size_t left_index = left.getIndex();
                const TopologyNode &right = root.getChild(1);
                size_t right_index = right.getIndex();
                const TopologyNode &middle = root.getChild(2);
                size_t middleIndex = middle.getIndex();

                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );
                fillLikelihoodVector( middle, middleIndex );

                dirty_nodes[root_index] = false;
                computeRootLikelihoodNode( root_index, left_index, right_index, middleIndex );
            
                // now call scaling
                scale(root_index, left_index, right_index, middleIndex);
            }
            else
            {
                throw RbException("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
            }

        }

        // sum the partials up
        this->lnProb = sumRootLikelihood();
    }
    else if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BOTH )
    {
        // only necessary if the root is actually dirty
        if ( dirty_nodes[root_index] == true )
        {

            // start by filling the likelihood vector for the children of the root
            if ( root.getNumberOfChildren() == 2 ) // rooted trees have two children for the root
            {
                const TopologyNode &left = root.getChild(0);
                size_t left_index = left.getIndex();
                const TopologyNode &right = root.getChild(1);
                size_t right_index = right.getIndex();

                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );

                dirty_nodes[root_index] = false;
                computeRootLikelihoodBranchNode( root_index, left_index, right_index );
            
                // now call scaling
                scale(root_index, left_index, right_index);

            }
            else if ( root.getNumberOfChildren() == 3 ) // unrooted trees have three children for the root
            {
                const TopologyNode &left = root.getChild(0);
                size_t left_index = left.getIndex();
                const TopologyNode &right = root.getChild(1);
                size_t right_index = right.getIndex();
                const TopologyNode &middle = root.getChild(2);
                size_t middleIndex = middle.getIndex();

                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );
                fillLikelihoodVector( middle, middleIndex );

                dirty_nodes[root_index] = false;
                computeRootLikelihoodBranchNode( root_index, left_index, right_index, middleIndex );

                // now call scaling
                scale(root_index, left_index, right_index, middleIndex);
            }
            else
            {
                throw RbException("The root node has an unexpected number of children. Only 2 (for rooted trees) or 3 (for unrooted trees) are allowed.");
            }

        }

        // sum the partials up
        this->lnProb = sumRootLikelihood();
    }

    // if we are not in MCMC mode, then we need to (temporarily) free memory
    if ( in_mcmc_mode == false )
    {
#if defined( RB_BEAGLE )
        //-- Clean up after ourselves when not in MCMC mode.
        charType tmp_char;
        if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
        {
            this->freeBeagleInstances();
        }
#endif
        delete [] partial_branch_likelihoods;
        delete [] partial_node_likelihoods;
        partial_branch_likelihoods = NULL;
        partial_node_likelihoods   = NULL;
    }

    // set the ancestral states as stale
    has_ancestral_states = false;

    return this->lnProb;
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeMarginalNodeLikelihood( size_t node_index, size_t parentnode_index )
{

    // compute the transition probability matrix
    this->updateTransitionProbabilities( node_index );

    // get the pointers to the partial likelihoods and the marginal likelihoods
    const double*   p_node                  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;
    double*         p_node_marginal         = this->marginalLikelihoods + node_index*this->node_offset;
    const double*   p_parent_node_marginal  = this->marginalLikelihoods + parentnode_index*this->node_offset;

    // get pointers the likelihood for both subtrees
    const double*   p_mixture                   = p_node;
    double*         p_mixture_marginal          = p_node_marginal;
    const double*   p_parent_mixture_marginal   = p_parent_node_marginal;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // the transition probability matrix for this mixture category
        const double*    tp_begin                = this->transition_prob_matrices[mixture].theMatrix;

        // get pointers to the likelihood for this mixture category
        const double*   p_site_mixture                  = p_mixture;
        double*         p_site_mixture_marginal         = p_mixture_marginal;
        const double*   p_parent_site_mixture_marginal  = p_parent_mixture_marginal;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get the pointers to the likelihoods for this site and mixture category
            const double*   p_site_j                    = p_site_mixture;
            double*         p_site_marginal_j           = p_site_mixture_marginal;
            // iterate over all end states
            for (size_t j=0; j<num_states; ++j)
            {
                const double*   p_parent_site_marginal_k    = p_parent_site_mixture_marginal;
                double sum = 0;

                // iterator over all start states
                for (size_t k=0; k<num_states; ++k)
                {
                    // transition probability for k->j
                    const double tp_kj = *p_parent_site_marginal_k * tp_begin[ k * num_states + j ];

                    // add the probability of starting from this state
                    sum += *p_site_j * tp_kj;

                    // next parent state
                    ++p_parent_site_marginal_k;
                }
                *p_site_marginal_j = sum;

                // increment pointers
                ++p_site_j; ++p_site_marginal_j;
            }

            // increment the pointers to the next site
            p_site_mixture+=this->site_offset; p_site_mixture_marginal+=this->site_offset; p_parent_site_mixture_marginal+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixture_offset; p_mixture_marginal+=this->mixture_offset; p_parent_mixture_marginal+=this->mixture_offset;

    } // end-for over all mixtures (=rate categories)

}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeMarginalRootLikelihood( void )
{
    // get the root node
    const TopologyNode &root = tau->getValue().getRoot();

    // get root frequencies
    std::vector<std::vector<double> >   ff;
    getRootFrequencies(ff);

    // get the index of the root node
    size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods and the marginal likelihoods
    const double*   p_node           = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;
    double*         p_node_marginal  = this->marginalLikelihoods + node_index*this->node_offset;

    // get pointers the likelihood for both subtrees
    const double*   p_mixture           = p_node;
    double*         p_mixture_marginal  = p_node_marginal;

    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        // get root frequencies
        const std::vector<double>&          f           = ff[mixture % ff.size()];
        assert(f.size() == num_states);
        std::vector<double>::const_iterator f_end       = f.end();
        std::vector<double>::const_iterator f_begin     = f.begin();

        // get pointers to the likelihood for this mixture category
        const double*   p_site_mixture          = p_mixture;
        double*         p_site_mixture_marginal = p_mixture_marginal;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get the pointer to the stationary frequencies
            std::vector<double>::const_iterator f_j             = f_begin;
            // get the pointers to the likelihoods for this site and mixture category
            const double*   p_site_j            = p_site_mixture;
            double*         p_site_marginal_j   = p_site_mixture_marginal;
            // iterate over all starting states
            for (; f_j != f_end; ++f_j)
            {
                // add the probability of starting from this state
                *p_site_marginal_j = *p_site_j * *f_j;

                // increment pointers
                ++p_site_j; ++p_site_marginal_j;
            }

            // increment the pointers to the next site
            p_site_mixture+=this->site_offset; p_site_mixture_marginal+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixture_offset; p_mixture_marginal+=this->mixture_offset;

    } // end-for over all mixtures (=rate categories)

}


/**
 * Draw a vector of ancestral states from the marginal distribution (non-conditional of the other ancestral states).
 * Here we assume that the marginal likelihoods have been updated.
 */
template<class charType>
std::vector<charType> RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::drawAncestralStatesForNode(const TopologyNode &node)
{

    size_t node_index = node.getIndex();

    // get the marginal likelihoods
    std::vector< std::vector<double> >* marginals = sumMarginalLikelihoods(node_index);

    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector< charType > ancestralSeq = std::vector<charType>();

    for ( size_t i = 0; i < num_sites; ++i )
    {
		size_t pattern = i;
		// if the matrix is compressed use the pattern for this site
		if ( compressed == true )
        {
            pattern = site_pattern[i];
        }

        // create the character
        charType c = charType( num_states );
        c.setToFirstState();

        // sum the likelihoods for each character state
        const std::vector<double> siteMarginals = (*marginals)[pattern];
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

                    if (c.getStateIndex() + 1 >= c.getNumberOfStates())
                    {
                        stateIndex = 0;
                        c.setToFirstState();
                    }
                    else
                    {
                        c++;
                        stateIndex++;
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

    // we need to free the vector
    delete marginals;

    return ancestralSeq;
}


/**
 * Draw a vector of ancestral states from the joint-conditional distribution of states.
 */
template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::drawJointConditionalAncestralStates(std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates)
{

	// if we already have ancestral states, don't make new ones
    
    // MJL 181028: Disabling this flag to allow multiple monitors to work for same dnPhyloCTMC (e.g. ancestral states + stochastic mapping)
//	if ( has_ancestral_states == true )
//    {
//		return;
//    }
    
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get working variables
    std::vector<double> siteProbVector = getMixtureProbs();

    const TopologyNode &root = tau->getValue().getRoot();
    size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods and the marginal likelihoods
    double*         p_node  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;

    // get pointers the likelihood for both subtrees
    const double*   p_site           = p_node;

    // sample root states
    std::vector<double> p( this->num_site_mixtures*this->num_states, 0.0);

    // clear the container for sampling the site-rates
    sampled_site_mixtures.resize(this->num_sites);
    
    for (size_t i = 0; i < this->num_sites; ++i)
//    for (size_t i = pattern_block_start; i < this->pattern_block_end; ++i)
    {

        // create the character
        charType c = charType( template_state );

        // sum to sample
        double sum = 0.0;

        // if the matrix is compressed use the pattern for this site
        size_t pattern = i - pattern_block_start;
        if ( compressed == true )
        {
            pattern = site_pattern[i - pattern_block_start];
        }

        // get ptr to first mixture cat for site
        p_site          = p_node  + pattern * this->site_offset;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {

            // get pointers to the likelihood for this mixture category
            const double* p_site_mixture_j       = p_site;

            // iterate over all starting states
            for (size_t state = 0; state < this->num_states; ++state)
            {
                size_t k = this->num_states*mixture + state;
                p[k] = *p_site_mixture_j * siteProbVector[mixture];
                sum += p[k];

                // increment the pointers to the next state for (site,rate)
                p_site_mixture_j++;
            }

            // increment the pointers to the next mixture category for given site
            p_site       += this->mixture_offset;

        } // end-for over all mixtures (=rate categories)

        // sample char from p
        bool stop = false;
        double u = rng->uniform01() * sum;
        for (size_t mixture = 0; mixture < this->num_site_mixtures; mixture++)
        {
            c.setToFirstState();
            for (size_t state = 0; state < this->num_states; state++)
            {
                size_t k = this->num_states * mixture + state;
                u -= p[k];
                if (u < 0.0)
                {
                    startStates[root.getIndex()][i] = c;
                    sampled_site_mixtures[i] = mixture;
                    stop = true;
                    break;
                }
                if (c.getStateIndex() + 1 >= c.getNumberOfStates())
                {
                    c.setToFirstState();
                }
                else
                {
                    c++;
                }
            }
            if (stop) break;
        }

        endStates[node_index][i] = startStates[node_index][i];
    }

    // recurse
    std::vector<TopologyNode*> children = root.getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        // daughters identically inherit ancestral state
        startStates[ children[i]->getIndex() ] = endStates[ root.getIndex() ];

        // recurse towards tips
        if ( children[i]->isTip() == false )
        {
            recursivelyDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampled_site_mixtures);
        }
        else
        {
            tipDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampled_site_mixtures);
        }

    }

    // flag the ancestral states as sampled
    has_ancestral_states = true;

}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::drawSiteMixtureAllocations()
{

    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get working variables
    std::vector<double> siteProbVector = getMixtureProbs();

    const TopologyNode &root = tau->getValue().getRoot();
    size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods and the marginal likelihoods
    double*         p_node  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;

    // get pointers the likelihood for both subtrees
    const double*   p_site           = p_node;

    // sample root states
    std::vector<double> p( this->num_site_mixtures*this->num_states, 0.0);

    // clear the container for sampling the site-rates
    this->sampled_site_mixtures.resize(this->num_sites);

    for (size_t i = 0; i < this->num_sites; ++i)
//    for (size_t i = pattern_block_start; i < this->pattern_block_end; ++i)
    {

        // sum to sample
        double sum = 0.0;

        // if the matrix is compressed use the pattern for this site
        size_t pattern = i - pattern_block_start;
        if ( compressed == true )
        {
            pattern = site_pattern[i - pattern_block_start];
        }

        // get ptr to first mixture cat for site
        p_site          = p_node  + pattern * this->site_offset;

        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {

            // get pointers to the likelihood for this mixture category
            const double* p_site_mixture_j       = p_site;

            // iterate over all starting states
            for (size_t state = 0; state < this->num_states; ++state)
            {
                size_t k = this->num_states * mixture + state;
                p[k] = *p_site_mixture_j * siteProbVector[mixture];
                sum += p[k];

                // increment the pointers to the next state for (site,rate)
                p_site_mixture_j++;
            }

            // increment the pointers to the next mixture category for given site
            p_site       += this->mixture_offset;

        } // end-for over all mixtures (=rate categories)

        // sample char from p
        bool stop = false;
        double u = rng->uniform01() * sum;
        for (size_t mixture = 0; mixture < this->num_site_mixtures; mixture++)
        {
            for (size_t state = 0; state < this->num_states; state++)
            {
                size_t k = this->num_states * mixture + state;
                u -= p[k];
                if (u < 0.0)
                {
                	this->sampled_site_mixtures[i] = mixture;
                    stop = true;
                    break;
                }
            }
            if (stop) break;
        }

    }

}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::drawStochasticCharacterMap(std::vector<std::string>& character_histories, size_t site, bool use_simmap_default)
{

    bool success = false;
    size_t max_draws = 10;
    size_t n_draws = 0;

    while (!success && n_draws != max_draws) {
        // first draw joint ancestral states
        std::vector<std::vector<charType> > start_states(this->num_nodes, std::vector<charType>(this->num_sites, template_state));
        std::vector<std::vector<charType> > end_states(this->num_nodes, std::vector<charType>(this->num_sites, template_state));
        this->drawJointConditionalAncestralStates( start_states, end_states );

        // save the character history for the root
        const TopologyNode &root = this->tau->getValue().getRoot();
        size_t root_index = root.getIndex();
        std::string simmap_string = "{" + end_states[root_index][site].getStringValue() + "," + StringUtilities::toString( root.getBranchLength() ) + "}";
        character_histories[root_index] = simmap_string;

        // get the sampled site-matrix and site-rate indexes
        getSampledMixtureComponents(site, sampled_site_rate_component, sampled_site_matrix_component);

        // recurse towards tips
        const TopologyNode &right = root.getChild(0);
        const TopologyNode &left = root.getChild(1);
        success = recursivelyDrawStochasticCharacterMap(left,  character_histories, start_states, end_states, site, use_simmap_default);
        success &= recursivelyDrawStochasticCharacterMap(right, character_histories, start_states, end_states, site, use_simmap_default);

        if (n_draws != 0) {
            std::cout << "Warning: numerical instability in P(t)=exp(Qt) caused stochastic mapping to fail (attempt: " << n_draws << "/" << max_draws << ")\n";
        }
        n_draws++;
    }

    if (n_draws == max_draws) {
        throw RbException("Stochastic mapping failed due to numerical instability.");
    }

}

template<class charType>
bool RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::hasSiteRateMixture()
{
	bool ret = this->num_site_rates > 1;
	return ret;
}

template<class charType>
bool RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::hasSiteMatrixMixture()
{
	bool ret = false;

	if ( this->site_matrix_probs != NULL ) {
		ret = this->site_matrix_probs->getValue().size() > 1;
	}

	return ret;
}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getSampledMixtureComponents(size_t &site_index, size_t &rate_component, size_t &matrix_component )
{

	// get the mixture component (in vector form)
	size_t mixture_component_index = sampled_site_mixtures[site_index];

	matrix_component = 0;
    if (this->site_matrix_probs != NULL)
    {
    	matrix_component = mixture_component_index % num_matrices;
    }

    // determine the rate (row index)
    rate_component = 0;
    if (this->site_rates != NULL)
    {
    	rate_component = (mixture_component_index - matrix_component) / num_matrices;
    }


}

template<class charType>
bool RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::recursivelyDrawStochasticCharacterMap(const TopologyNode &node, std::vector<std::string>& character_histories, std::vector<std::vector<charType> >& start_states, std::vector<std::vector<charType> >& end_states, size_t site, bool use_simmap_default)
{

    bool success = false;

    // get the start and end states
    size_t node_index = node.getIndex();
    size_t start_state = start_states[node_index][site].getStateIndex();
    size_t end_state = start_state;

    // NOTE: ambiguous tip states are sampled along with internal node states
    end_state = end_states[node_index][site].getStateIndex();


    // set up vectors to hold the character transition events
    std::vector<size_t> transition_states;
    std::vector<double> transition_times;
    transition_states.push_back(start_state);

    // get the rate matrix for this branch (or site if using a mixture of matrices over sites)
    RateMatrix_JC jc(this->num_states); // BAD!!
    const RateGenerator *rate_matrix = &jc;
    if ( this->branch_heterogeneous_substitution_matrices == true )
    {
        if (this->heterogeneous_rate_matrices != NULL)
        {
            rate_matrix = &this->heterogeneous_rate_matrices->getValue()[node_index];
        }
        else if (this->homogeneous_rate_matrix != NULL)
        {
            rate_matrix = &this->homogeneous_rate_matrix->getValue();
        }
    }
    else
    {
        if (this->homogeneous_rate_matrix != NULL)
        {
            rate_matrix = &this->homogeneous_rate_matrix->getValue();
        }
        else if (this->site_matrix_probs != NULL)
        {
        	rate_matrix = &this->heterogeneous_rate_matrices->getValue()[this->sampled_site_matrix_component];
        }
    }

    // get the clock rate for the branch
    double clock_rate = 1.0;
    if ( this->branch_heterogeneous_clock_rates == true )
    {
        if (this->heterogeneous_clock_rates != NULL)
        {
            clock_rate = this->heterogeneous_clock_rates->getValue()[node_index];
        }
    }
    else
    {
        if (this->homogeneous_clock_rate != NULL)
        {
            clock_rate = this->homogeneous_clock_rate->getValue();
        }
    }

    // simulate stochastic map
    transition_states.push_back(end_state);
    if (mixture_model)
    {
        auto& model = mixture_model->getValue();
        int m = sampled_site_rate_component;
        success = const_cast<SiteMixtureModel&>(model).getComponent(m).simulateStochasticMapping(tau->getValue(), node_index, clock_rate, transition_states, transition_times);
    }
    else
    {
        auto [start_age, end_age] = getStartEndAge(node);

        // multiply by the clock-rate for the site
        if (this->site_rates != NULL)
        {
            // there is a mixture over site rates
            clock_rate *= this->site_rates->getValue()[this->sampled_site_rate_component];
        }

        success = const_cast<RateGenerator*>(rate_matrix)->simulateStochasticMapping(start_age, end_age, clock_rate, transition_states, transition_times);
    }

    // make SIMMAP string
    std::string simmap_string = "{";

    if (use_simmap_default == true)
    {
        for (size_t i = transition_times.size(); i > 0; i--)
        {
            simmap_string = simmap_string + StringUtilities::toString(transition_states[i - 1]) + "," + StringUtilities::toString(transition_times[i - 1]);
            if (i != 1)
            {
                simmap_string = simmap_string + ":";
            }
        }
    }
    else
    {
        for (size_t i = 0; i < transition_times.size(); i++)
        {
            if (i != 0)
            {
                simmap_string = simmap_string + ":";
            }
            simmap_string = simmap_string + StringUtilities::toString(transition_states[i]) + "," + StringUtilities::toString(transition_times[i]);
        }
    }
    simmap_string = simmap_string + "}";

    // save the character history for this branch
    character_histories[node_index] = simmap_string;

    // recurse towards tips
    if ( node.isTip() == false )
    {
        const TopologyNode &right = node.getChild(0);
        const TopologyNode &left = node.getChild(1);
        success &= recursivelyDrawStochasticCharacterMap(left, character_histories, start_states, end_states, site, use_simmap_default);
        success &= recursivelyDrawStochasticCharacterMap(right, character_histories, start_states, end_states, site, use_simmap_default);
    }

    return success;
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    if ( n == "siteLikelihoods" )
    {

        bool delete_partial_branch_likelihoods = false;

        // if we are not in MCMC mode, then we need to (temporarily) allocate memory
        if ( in_mcmc_mode == false )
        {
            delete_partial_branch_likelihoods = true;
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
            in_mcmc_mode = true;

            for (std::vector<bool>::iterator it = dirty_branches.begin(); it != dirty_branches.end(); ++it)
            {
                (*it) = true;
            }
            
            for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
            {
                (*it) = true;
            }
        }

        // make sure the likelihoods are updated
        const_cast<AbstractPhyloCTMCSiteHomogeneous<charType> *>( this )->computeLnProbability();

        // get the per site likelihood
        RbVector<double> tmp = RbVector<double>(num_patterns, 0.0);
        computeRootLikelihoods( tmp );

        // if we are not in MCMC mode, then we need to (temporarily) free memory
        if ( delete_partial_branch_likelihoods == true )
        {
            // free the partial likelihoods
            delete [] partial_branch_likelihoods;
            partial_branch_likelihoods = NULL;
            in_mcmc_mode = false;
        }

        // now match it back to the actual sites
        if ( this->compressed == true && num_sites > num_patterns )
        {
            rv = RbVector<double>(num_sites, 0.0);

            for (size_t i=0; i<num_sites; ++i)
            {
                size_t pattern_index = site_pattern[i];
                rv[i] = tmp[pattern_index] / pattern_counts[pattern_index];
            }
        }
        else
        {
            rv = tmp;
        }

    }
    else if ( n == "siteRates" )
    {

        bool delete_partial_branch_likelihoods = false;

        // if we are not in MCMC mode, then we need to (temporarily) allocate memory
        if ( in_mcmc_mode == false )
        {
            delete_partial_branch_likelihoods = true;
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
            in_mcmc_mode = true;

            for (std::vector<bool>::iterator it = dirty_branches.begin(); it != dirty_branches.end(); ++it)
            {
                (*it) = true;
            }
            for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
            {
                (*it) = true;
            }
        }

        // make sure the likelihoods are updated
        const_cast<AbstractPhyloCTMCSiteHomogeneous<charType> *>( this )->computeLnProbability();

        // get the site rates
        std::vector<double> r;
        if ( this->rate_variation_across_sites == true )
        {
            r = this->site_rates->getValue();
        }
        else
        {
            r.push_back(1.0);
        }

        size_t num_site_rates_withInv = num_site_rates;
        double prob_invariant = getPInv();
        if (prob_invariant > 0.0)
        {
            num_site_rates_withInv++;
            r.insert(r.begin(), 0.0);
        }

        // get the per site rate likelihood
        MatrixReal tmp = MatrixReal(num_patterns, num_site_rates_withInv, 0.0);
        computeRootLikelihoodsPerSiteRate( tmp );

        // if we are not in MCMC mode, then we need to (temporarily) free memory
        if ( delete_partial_branch_likelihoods == true )
        {
            // free the partial likelihoods
            delete [] partial_branch_likelihoods;
            partial_branch_likelihoods = NULL;
            in_mcmc_mode = false;
        }

        // now match it back to the actual sites
        MatrixReal siteRateConditionalProb = MatrixReal(num_sites, num_site_rates_withInv, 0.0);
        for (size_t i=0; i<num_sites; ++i)
        {
            size_t pattern_index = site_pattern[i];
            double siteLikelihoods = 0.0;

            for (size_t j=0; j<num_site_rates_withInv; ++j)
            {
                siteRateConditionalProb[i][j] = exp(tmp[pattern_index][j] / pattern_counts[pattern_index]);
                siteLikelihoods += siteRateConditionalProb[i][j];
            }
            for (size_t j=0; j<num_site_rates_withInv; ++j)
            {
                siteRateConditionalProb[i][j] = siteRateConditionalProb[i][j] / siteLikelihoods;
            }
        }

        rv = RbVector<double>(num_sites, 0.0);

        // now either sample the site rates according to the conditional posterior probability of each rate category
        // or compute them as the posterior mean
        std::string method = static_cast<const TypedDagNode<std::string>* >( args[0] )->getValue();

        if (method == "sampling")
        {
            RandomNumberGenerator* rng = GLOBAL_RNG;
            for (size_t i=0; i<num_sites; ++i)
            {
                size_t rateIndex = 0;
                double u = rng->uniform01();

                while ( true )
                {
                    u -= siteRateConditionalProb[i][rateIndex];

                    if ( u > 0.0 )
                    {
                        ++rateIndex;
                    }
                    else
                    {
                        break;
                    }
                }

                rv[i] = r[rateIndex];

            }

        }
        else if (method == "weightedAverage")
        {
            for (size_t i=0; i<num_sites; ++i)
            {
                for (size_t rateIndex=0; rateIndex < num_site_rates_withInv; ++rateIndex)
                {
                    rv[i] += r[rateIndex] * siteRateConditionalProb[i][rateIndex];
                }
            }
        }

    }
    else
    {
        throw RbException() << "The PhyloCTMC process does not have a member method called '" << n << "'.";
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, MatrixReal &rv) const
{

    if ( n == "siteRateLikelihoods" )
    {

        bool delete_partial_branch_likelihoods = false;

        // if we are not in MCMC mode, then we need to (temporarily) allocate memory
        if ( in_mcmc_mode == false )
        {
            delete_partial_branch_likelihoods = true;
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
            in_mcmc_mode = true;

            for (std::vector<bool>::iterator it = dirty_branches.begin(); it != dirty_branches.end(); ++it)
            {
                (*it) = true;
            }
            for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
            {
                (*it) = true;
            }
        }

        // make sure the likelihoods are updated
        const_cast<AbstractPhyloCTMCSiteHomogeneous<charType> *>( this )->computeLnProbability();

        // get the per site rate likelihood
        size_t num_site_rates_withInv = num_site_rates;
        double prob_invariant = getPInv();
        if (prob_invariant > 0.0)
        {
            num_site_rates_withInv++;
        }

        // get the per site rate likelihood
        MatrixReal tmp = MatrixReal(num_patterns, num_site_rates_withInv, 0.0);
        computeRootLikelihoodsPerSiteRate( tmp );

        // if we are not in MCMC mode, then we need to (temporarily) free memory
        if ( delete_partial_branch_likelihoods == true )
        {
            // free the partial likelihoods
            delete [] partial_branch_likelihoods;
            partial_branch_likelihoods = NULL;
            in_mcmc_mode = false;
        }

        // now match it back to the actual sites
        rv = MatrixReal(num_sites, num_site_rates_withInv, 0.0);
        for (size_t i=0; i<num_sites; ++i)
        {
            size_t pattern_index = site_pattern[i];
            for (size_t j=0; j<num_site_rates_withInv; ++j)
            {
                rv[i][j] = tmp[pattern_index][j] / pattern_counts[pattern_index];
            }
        }

    }
    else if ( n == "siteMixtureLikelihoods" )
    {

        bool delete_partial_branch_likelihoods = false;

        // if we are not in MCMC mode, then we need to (temporarily) allocate memory
        if ( in_mcmc_mode == false )
        {
            delete_partial_branch_likelihoods = true;
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];
            in_mcmc_mode = true;

            for (std::vector<bool>::iterator it = dirty_branches.begin(); it != dirty_branches.end(); ++it)
            {
                (*it) = true;
            }
            for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
            {
                (*it) = true;
            }
        }

        // make sure the likelihoods are updated
        const_cast<AbstractPhyloCTMCSiteHomogeneous<charType> *>( this )->computeLnProbability();

        // get the per site rate likelihood
        size_t num_site_mixture_withInv = num_site_mixtures;
        double prob_invariant = getPInv();
        if (prob_invariant > 0.0)
        {
            num_site_mixture_withInv += size_t(num_site_mixtures/num_site_rates);
        }

        // get the per site rate likelihood
        MatrixReal tmp = MatrixReal(num_patterns, num_site_mixture_withInv, 0.0);
        computeRootLikelihoodsPerSiteMixture( tmp );

        // if we are not in MCMC mode, then we need to (temporarily) free memory
        if ( delete_partial_branch_likelihoods == true )
        {
            // free the partial likelihoods
            delete [] partial_branch_likelihoods;
            partial_branch_likelihoods = NULL;
            in_mcmc_mode = false;
        }

        // now match it back to the actual sites
        rv = MatrixReal(num_sites, num_site_mixture_withInv, 0.0);
        for (size_t i=0; i<num_sites; ++i)
        {
            size_t pattern_index = site_pattern[i];
            for (size_t j=0; j<num_site_mixture_withInv; ++j)
            {
                rv[i][j] = tmp[pattern_index][j] / pattern_counts[pattern_index];
            }
        }

    }
    else
    {
        throw RbException() << "The PhyloCTMC process does not have a member method called '" << n << "'.";
    }

}




template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::recursivelyDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates, const std::vector<size_t>& sampledSiteRates)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get working variables
    size_t node_index = node.getIndex();
    size_t left = node.getChild(0).getIndex();
    size_t right = node.getChild(1).getIndex();
    //    size_t parentIndex = node.getParent().getIndex();

    // get transition probabilities
    this->updateTransitionProbabilities( node_index );

    // get the pointers to the partial likelihoods and the marginal likelihoods
    const double*   p_left  = this->partial_branch_likelihoods + this->active_branch_likelihood[left]  * this->active_branch_likelihood_offset + left  * this->node_offset;
    const double*   p_right = this->partial_branch_likelihoods + this->active_branch_likelihood[right] * this->active_branch_likelihood_offset + right * this->node_offset;

    // sample characters conditioned on start states, going to end states
    std::vector<double> p(this->num_states, 0.0);
    for (size_t i = 0; i < this->num_sites; i++)
    {
        size_t cat = sampledSiteRates[i];
        size_t k = startStates[node_index][i].getStateIndex();


        // sum to sample
        double sum = 0.0;

        // if the matrix is compressed use the pattern for this site
        size_t pattern = i;
        if ( compressed == true )
        {
            pattern = site_pattern[i];
        }

        // get ptr to first mixture cat for site
        //        p_site          = p_node  + cat * this->mixture_offset + pattern * this->site_offset;
        const double* p_left_site_mixture_j     = p_left  + cat * this->mixture_offset + pattern * this->site_offset;
        const double* p_right_site_mixture_j    = p_right + cat * this->mixture_offset + pattern * this->site_offset;

        // iterate over possible end states for each site given start state
        for (size_t j = 0; j < this->num_states; j++)
        {
            double tp_kj = this->transition_prob_matrices[cat][k][j];
            p[j] = tp_kj * *p_left_site_mixture_j * *p_right_site_mixture_j;
            sum += p[j];

            //            p_site_mixture_j++;
            p_left_site_mixture_j++;
            p_right_site_mixture_j++;
        }

        // sample char from p
        charType c = charType( template_state );
        double u = rng->uniform01() * sum;
        for (size_t state = 0; state < this->num_states; state++)
        {
            u -= p[state];
            if (u < 0.0)
            {
                endStates[node_index][i] = c;
                break;
            }
            if (c.getStateIndex() + 1 >= c.getNumberOfStates())
            {
                c.setToFirstState();
            }
            else
            {
                c++;
            }
        }
    }

    // recurse
    std::vector<TopologyNode*> children = node.getChildren();
    for (size_t i = 0; i < children.size(); i++)
    {
        // daughters identically inherit ancestral state
        startStates[ children[i]->getIndex() ] = endStates[ node.getIndex() ];

        // recurse towards tips
        if ( children[i]->isTip() == false )
        {
            recursivelyDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampledSiteRates);
        }
        else
        {
            tipDrawJointConditionalAncestralStates(*children[i], startStates, endStates, sampledSiteRates);
        }

    }

}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::tipDrawJointConditionalAncestralStates(const TopologyNode &node, std::vector<std::vector<charType> >& startStates, std::vector<std::vector<charType> >& endStates, const std::vector<size_t>& sampledSiteRates)
{

    // get working variables
    size_t node_index = node.getIndex();

    // get transition probabilities
    this->updateTransitionProbabilities( node_index );

    const AbstractHomologousDiscreteCharacterData& d = this->getValue();
    const HomologousDiscreteCharacterData<charType>& dd = static_cast<const HomologousDiscreteCharacterData<charType>& >( d );
    const DiscreteTaxonData<charType>& td = dd.getTaxonData( node.getName() );

    // ideally sample ambiguous tip states given the underlying process and ancestral state
    // for now, always sample the clamped character
    
    // create a vector with the correct site indices
    // some of the sites may have been excluded
    std::vector<size_t> site_indices = getIncludedSiteIndices();

    // sample characters conditioned on start states, going to end states
    std::vector<double> p(this->num_states, 0.0);
    for (size_t i = 0; i < this->num_sites; i++)
    {
        charType c = td.getCharacter(site_indices[i]);

        if ( c.isAmbiguous() == false )
        {
            // we know the tip state for unambiguous characters
            endStates[node_index][i] = c;
        }
        else
        {
            // we sample the tip state for ambiguous characters
            size_t cat = sampledSiteRates[i];
            size_t k = startStates[node_index][i].getStateIndex();

            // sum to sample
            double sum = 0.0;

            // if the matrix is compressed use the pattern for this site
            size_t pattern = i;
            if ( compressed == true )
            {
                pattern = site_pattern[i];
            }

            // get the ambiguous character's bitset for the tip taxon
            RbBitSet bs = RbBitSet(this->num_states);
            if ( c.isMissingState() )
                bs.set(); // set to all 1s.
            else
                bs = c.getState();

            // iterate over possible end states for each site given start state
            for (size_t j = 0; j < this->num_states; j++)
            {
                double tp_kj = this->transition_prob_matrices[cat][k][j] * bs[j];
                p[j] = tp_kj;
                sum += p[j];
            }

            // sample char from p
            charType c = charType( template_state );
            double u = GLOBAL_RNG->uniform01() * sum;
            for (size_t state = 0; state < this->num_states; state++)
            {
                u -= p[state];
                if (u < 0.0)
                {
                    endStates[node_index][i] = c;
                    break;
                }
                if (c.getStateIndex() + 1 >= c.getNumberOfStates())
                {
                    c.setToFirstState();
                }
                else
                {
                    c++;
                }
            }
        }
    }

    // no further recursion

}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::fillLikelihoodVector(const TopologyNode &node, size_t node_index)
{

    if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BRANCH )
    {

        // check for recomputation
        if ( dirty_branches[node_index] == true )
        {
            // mark as computed
            dirty_branches[node_index] = false;

            if ( node.isTip() == true )
            {
                // this is a tip node
                // compute the likelihood for the tip and we are done
                computeTipLikelihood(node, node_index);

                //-- We only need to rescale the tip likelihoods if we are not using BEAGLE
#if !defined ( RB_BEAGLE )
                
                // rescale likelihood vector
                scale(node_index);
#endif /* NOT RB_BEAGLE */
            }
            else
            {
                // this is an internal node
                const TopologyNode     &left        = node.getChild(0);
                size_t                  left_index  = left.getIndex();
                const TopologyNode     &right       = node.getChild(1);
                size_t                  right_index = right.getIndex();
                
                fillLikelihoodVector( left, left_index );
                fillLikelihoodVector( right, right_index );

                // now compute the likelihoods of this internal node
                computeInternalNodeLikelihoodBranchWise(node,node_index,left_index,right_index);
            
            
                //-- We only need to scale if we are not using BEAGLE
#if !defined ( RB_BEAGLE )
                // rescale likelihood vector
                scale(node_index,left_index,right_index);
#endif /* NOT RB_BEAGLE */

            }

        }
    }
    else if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::NODE )
    {
        // check for recomputation
        if ( dirty_nodes[node_index] == true && node.isTip() == false )
        {
            // mark as computed
            dirty_nodes[node_index] = false;

            // this is an internal node
            const TopologyNode     &left        = node.getChild(0);
            size_t                  left_index  = left.getIndex();
            const TopologyNode     &right       = node.getChild(1);
            size_t                  right_index = right.getIndex();
                
            fillLikelihoodVector( left, left_index );
            fillLikelihoodVector( right, right_index );

            // now compute the likelihoods of this internal node
            computeInternalNodeLikelihoodNodeWise(node,node_index,left_index,right_index);
            
            
            //-- We only need to scale only if we are not using BEAGLE
#if !defined ( RB_BEAGLE )
            // rescale likelihood vector
            scale(node_index,left_index,right_index);
#endif /* NOT RB_BEAGLE */

        }


    }
    else if ( partial_likelihood_storing_approach == PARTIAL_LIKELIHOOD_STORING::BOTH )
    {
        
        // check for recomputation
        if ( dirty_nodes[node_index] == true && node.isTip() == false )
        {
            // mark as computed
            dirty_nodes[node_index] = false;
            
            // this is an internal node
            const TopologyNode     &left        = node.getChild(0);
            size_t                  left_index  = left.getIndex();
            const TopologyNode     &right       = node.getChild(1);
            size_t                  right_index = right.getIndex();
                
            fillLikelihoodVector( left, left_index );
            fillLikelihoodVector( right, right_index );

            // now compute the likelihoods of this internal node
            computeInternalNodeLikelihoodBranchNodeWise(node,node_index,left_index,right_index);
                
            //-- We only need to scale only if we are not using BEAGLE
#if !defined ( RB_BEAGLE )
            // rescale likelihood vector
            scale(node_index,left_index,right_index);
#endif /* NOT RB_BEAGLE */

        }
    }

}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::fillTipLikelihoods( void )
{
    delete tip_likelihoods;
    
    // find the unique site patterns and compute their respective frequencies
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
    
    size_t number_of_nodes = nodes.size();
    
    tip_likelihoods = new double[num_tips*pattern_block_size*num_states];

    // reinitialize likelihood vectors
    for (size_t n = 0; n < number_of_nodes; ++n)
    {
        TopologyNode* node = nodes[n];
        
        if ( node->isTip() )
        {

            size_t data_tip_index = this->taxon_name_2_tip_index_map[ node->getName() ];

            const std::vector<bool> &gap_node = this->gap_matrix[data_tip_index];
            const std::vector<std::uint64_t> &char_node = this->char_matrix[data_tip_index];
            const std::vector<RbBitSet> &amb_char_node = this->ambiguous_char_matrix[data_tip_index];

            size_t char_data_node_index = this->value->indexOfTaxonWithName(node->getName());

            std::vector<size_t> site_indices;
            if ( this->using_weighted_characters == true )
                site_indices = this->getIncludedSiteIndices();

            double* p_node = tip_likelihoods + data_tip_index*(pattern_block_size*num_states);

            for (size_t site = 0; site < pattern_block_size; ++site)
            {
                
                // is this site a gap?
                if ( gap_node[site] )
                {
                    // since this is a gap we need to assume that the actual state could have been any state

                    // iterate over all initial states for the transitions
                    for (size_t c1 = 0; c1 < this->num_states; ++c1)
                    {
                        // store the likelihood
                        p_node[c1] = 1.0;
                    }
                }
                else // we have observed a character
                {

                    if ( this->using_ambiguous_characters == true && this->using_weighted_characters == false)
                    {
                        
                        // note, the observed state could be ambiguous!
                        const RbBitSet &val = amb_char_node[site];
                        
                        // iterate over all possible initial states
                        for (size_t c1 = 0; c1 < this->num_states; ++c1)
                        {
                            
                            // check whether we observed this state
                            if ( val.test(c1) == true )
                            {
                                // set the likelihood
                                p_node[c1] = 1.0;
                            }
                            else
                            {
                                // set the likelihood
                                p_node[c1] = 0.0;
                            }
                            
                        }
                    }
                    else if ( this->using_weighted_characters == true )
                    {
                        
                        // note, the observed state could be ambiguous!
                        size_t this_site_index = site_indices[site];
                        const std::vector< double >& weights = this->value->getCharacter(char_data_node_index, this_site_index).getWeights();

                        // iterate over all possible initial states
                        for (size_t c1 = 0; c1 < this->num_states; ++c1)
                        {
                            // set the weight
                            p_node[c1] = weights[c1];
                        }
                            
                    }
                    else // no ambiguous characters in use
                    {
                            
                        std::uint64_t org_val = char_node[site];
                        
                        // iterate over all possible initial states
                        for (size_t c1 = 0; c1 < this->num_states; ++c1)
                        {
                            
                            // check whether we observed this state
                            if ( org_val == c1 )
                            {
                                // set the likelihood
                                p_node[c1] = 1.0;
                            }
                            else
                            {
                                // set the likelihood
                                p_node[c1] = 0.0;
                            }
                            
                        }

                    } // end-for over all possible characters

                } // end-if a gap state
                
                p_node += num_states;
            } // end-for over all sites
            
        } // end-if this is a tip node

    } // end-for over all nodes
    
}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::fireTreeChangeEvent( const RevBayesCore::TopologyNode &n, const unsigned& m )
{

    // call a recursive flagging of all node above (closer to the root) and including this node
    if ( m == RevBayesCore::TreeChangeEventMessage::TOPOLOGY )
    {
        recursivelyFlagBranchDirty( n );
        recursivelyFlagNodeDirty( n );
    }
    else
    {
        recursivelyFlagBranchDirty( n );

        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }
    }
    
    if (m != RevBayesCore::TreeChangeEventMessage::TOPOLOGY)
    {
        recursivelyFlagBranchDirty( n );
        recursivelyFlagNodeDirty( n );
        flagNodeDirtyPmatrix(n.getIndex());
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::flagNodeDirtyPmatrix(size_t node_idx)
{
    pmat_dirty_branches[node_idx] = true;
    if ( pmat_changed_nodes[node_idx] == false )
    {
        active_pmatrices[node_idx] = (active_pmatrices[node_idx] == 0 ? 1 : 0);
        pmat_changed_nodes[node_idx] = true;
    }
    
}


template<class charType>
std::vector<size_t> RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getIncludedSiteIndices( void )
{
    return this->value->getIncludedSiteIndices();
}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getRootFrequencies( std::vector<std::vector<double> >& rf ) const
{
    if ( root_frequencies != NULL )
    {
        assert(not mixture_model);

        std::vector<double> f = root_frequencies->getValue();
        rf.push_back( f );
    }
    else if (mixture_model)
    {
        rf.clear();

        auto& mm = mixture_model->getValue();
        for(int i=0;i < mm.getNumberOfComponents(); i++)
            rf.push_back( mm.getComponent(i).getRootFrequencies() );
    }
    else if (heterogeneous_rate_matrices !=  NULL)
    {
        if ( this->branch_heterogeneous_substitution_matrices == true)
        {
            if (root_frequencies == NULL)
            {
                throw RbException("Using branch-heterogeneous rate matrices, but no root frequencies have been specified");
            }

            std::vector<double> f = root_frequencies->getValue();
            rf.push_back( f );
        }
        else
        {
            for (size_t matrix = 0; matrix < this->num_matrices; matrix++)
            {
                const RateMatrix *rm = dynamic_cast<const RateMatrix *>(&heterogeneous_rate_matrices->getValue()[matrix]);
                if ( rm != NULL )
                {
                    rf.push_back( rm->getStationaryFrequencies() );
                }
                else
                {
                    throw RbException("If you want to use RateGenerators that are not RateMatrices then you need to specify the root frequencies directly.");
                }
            }
        }
    }
    else if (homogeneous_rate_matrix != NULL)
    {
        const RateMatrix *rm = dynamic_cast<const RateMatrix *>(&homogeneous_rate_matrix->getValue());
        if ( rm != NULL )
        {
            rf.push_back( rm->getStationaryFrequencies() );
        }
        else
        {
            throw RbException("If you want to use RateGenerators that are not RateMatrices then you need to specify the root frequencies directly.");
        }

    }
    else
    {
        rf.push_back( std::vector<double>(num_states, 1.0/num_states) );
    }
}

template<class charType>
std::vector<double> RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getRootFrequencies( size_t mixture ) const
{
    if (mixture > this->num_site_mixtures)
    {
        throw(RbException("Site mixture index out of bounds"));
    }

    std::vector<std::vector<double> > rf;
    getRootFrequencies(rf);

    return rf[mixture % rf.size()];
}

template<class charType>
std::vector<double> RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getMixtureProbs( void ) const
{
    if (mixture_model)
        return mixture_model->getValue().componentProbs();

    std::vector<double> probs(num_site_mixtures, 1.0/num_site_mixtures);
    std::vector<double> rates_probs(num_site_rates, 1.0/num_site_rates);
    size_t num_site_matrices = num_site_mixtures/num_site_rates;
    std::vector<double> matrix_probs(num_site_matrices, 1.0/num_site_matrices);

    if ( site_rates_probs != NULL )
    {
        rates_probs = site_rates_probs->getValue();
    }

    if ( site_matrix_probs != NULL )
    {
        matrix_probs = site_matrix_probs->getValue();
    }

    for (size_t matrix = 0; matrix < num_site_matrices; ++matrix)
    {
        for (size_t j = 0; j < this->num_site_rates; ++j)
        {
            probs[j * num_site_matrices + matrix] = matrix_probs[matrix] * rates_probs[j];
        }
    }

    return probs;
}


template<class charType>
double RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getPInv( void ) const
{

    if ( p_inv != NULL )
    {
        return p_inv->getValue();
    }
    else
    {
        return 0.0;
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::keepSpecialization( const DagNode* affecter )
{

    // reset flags for likelihood computation
    touched = false;

    // we don't have a previous state anymore.
    this->storedLnProb = {};

    // reset all flags
    prev_dirty_nodes = {};
    prev_dirty_branches = {};
    prev_pmat_dirty_branches = {};

    for (std::vector<bool>::iterator it = this->changed_branches.begin(); it != this->changed_branches.end(); ++it)
    {
        (*it) = false;
    }

    for (std::vector<bool>::iterator it = this->changed_nodes.begin(); it != this->changed_nodes.end(); ++it)
    {
        (*it) = false;
    }
    
    for (std::vector<bool>::iterator it = this->pmat_changed_nodes.begin(); it != this->pmat_changed_nodes.end(); ++it)
    {
        (*it) = false;
    }


#if defined( RB_BEAGLE )
    // Reset eigensystems for BEAGLE
    for (std::vector<bool>::iterator it = this->touched_eigen_system.begin(); it != this->touched_eigen_system.end(); ++it)
    {
        (*it) = false;
    }
    
//    bool dest = this->active_branch_likelihood[this->tau->getValue().getRoot().getIndex()];
//    int destScalingIndex  = (int) 2*this->num_nodes+(dest ? 1 : 0);
//    int srcScalingIndex   = (int) 2*this->num_nodes+(dest ? 0 : 1);
//    beagleCopyScaleFactors(this->beagle_instance->getResourceID(), destScalingIndex, srcScalingIndex);

    
#endif /* RB_BEAGLE */
}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::recursivelyFlagBranchDirty( const RevBayesCore::TopologyNode &n )
{

    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();

    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( dirty_branches[index] == false )
    {
        // the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagBranchDirty( n.getParent() );
        }

        // set the flag
        dirty_branches[index] = true;

        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changed_branches[index] == false )
        {
            active_branch_likelihood[index] = (active_branch_likelihood[index] == 0 ? 1 : 0);
            changed_branches[index] = true;
        }

    }
}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::recursivelyFlagNodeDirty( const RevBayesCore::TopologyNode &n )
{

    // we need to flag this node and all ancestral nodes for recomputation
    size_t index = n.getIndex();

    // if this node is already dirty, the also all the ancestral nodes must have been flagged as dirty
    if ( dirty_nodes[index] == false )
    {
        // the root doesn't have an ancestor
        if ( n.isRoot() == false )
        {
            recursivelyFlagNodeDirty( n.getParent() );
        }

        // set the flag
        dirty_nodes   [index] = true;

        // if we previously haven't touched this node, then we need to change the active likelihood pointer
        if ( changed_nodes[index] == false )
        {
            active_node_likelihood[index] = (active_node_likelihood[index] == 0 ? 1 : 0);
            changed_nodes[index] = true;
        }

    }
    else
    {
#ifndef NDEBUG
	// Check the invariant that if a node is dirty, its parent is also dirty.
	auto n2 = &n;
	while(not n2->isRoot())
	{
	    n2 = &n2->getParent();
	    assert(dirty_nodes[n2->getIndex()]);
	}
#endif
    }

}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::recursiveMarginalLikelihoodComputation( size_t node_index )
{

    const TopologyNode &node = tau->getValue().getNode( node_index );

    for ( size_t i=0; i<node.getNumberOfChildren(); ++i )
    {
        const TopologyNode &child = node.getChild(i);

        if ( child.isTip() == false )
        {
            size_t childIndex = child.getIndex();
            computeMarginalNodeLikelihood( childIndex, node_index );
            recursiveMarginalLikelihoodComputation( childIndex );
        }

    }
}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::redrawValue( void )
{
    bool do_mask = this->dag_node != NULL && this->dag_node->isClamped() && gap_match_clamped;
    std::vector<std::vector<bool> > mask_gap        = std::vector<std::vector<bool> >(tau->getValue().getNumberOfTips(), std::vector<bool>());
    std::vector<std::vector<bool> > mask_missing    = std::vector<std::vector<bool> >(tau->getValue().getNumberOfTips(), std::vector<bool>());
    
    // we cannot use the stored gap matrix because it uses the pattern compression
    // therefore we create our own mask
    if ( do_mask == true )
    {
        this->value->fillMissingSitesMask(mask_gap, mask_missing);
    }

    // delete the old value first
    delete this->value;

    // create a new character data object
    this->value = new HomologousDiscreteCharacterData<charType>();

    // create a vector of taxon data
    std::vector< DiscreteTaxonData<charType> > taxa = std::vector< DiscreteTaxonData< charType > >( num_nodes, DiscreteTaxonData<charType>( Taxon("") ) );

    // first, simulate the per site rates
    RandomNumberGenerator* rng = GLOBAL_RNG;
    std::vector<size_t> perSiteMixtures = std::vector<size_t>(num_sites,0);
    std::vector<bool> inv = std::vector<bool>(num_sites,false);
    double prob_invariant = getPInv();
    for ( size_t i = 0; i < num_sites; ++i )
    {
        // draw if this site is invariant
        double u = rng->uniform01();
        if ( u < prob_invariant )
        {
            // this site is invariant
            inv[i] = true;

        }
        else if ( num_site_mixtures  > 1 )
        {
            // draw the rate for this site
            u = rng->uniform01();
            size_t mixtureIndex = 0;

            std::vector<double> mixtureProbs = getMixtureProbs();

            std::vector< double >::const_iterator freq = mixtureProbs.begin();
            while ( true )
            {
                u -= *freq;

                if ( u > 0.0 )
                {
                    ++mixtureIndex;
                    ++freq;
                }
                else
                {
                    break;
                }
            }

            perSiteMixtures[i] = mixtureIndex;
        }
        else
        {
            // there is only a single site rate so this is 1.0
            perSiteMixtures[i] = 0;

        }

    }

    std::vector<std::vector<double> > freqs;
    getRootFrequencies(freqs);
    // simulate the root sequence
    DiscreteTaxonData< charType > &root = taxa[ tau->getValue().getRoot().getIndex() ];
    for ( size_t i = 0; i < num_sites; ++i )
    {
        const std::vector< double > &stationary_freqs = freqs[perSiteMixtures[i] % freqs.size()];

        // create the character
        charType c = charType( num_states );
        c.setToFirstState();

        // draw the state
        double u = rng->uniform01();
        std::vector< double >::const_iterator freq = stationary_freqs.begin();
        while ( true )
        {
            u -= *freq;
            if ( u > 0.0 )
            {
                ++c;
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
    simulate( tau->getValue().getRoot(), taxa, inv, perSiteMixtures );

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
        this->value->applyMissingSitesMask(mask_gap, mask_missing);
    }

    // compress the data and initialize internal variables
    compress();

    for (std::vector<bool>::iterator it = dirty_branches.begin(); it != dirty_branches.end(); ++it)
    {
        (*it) = true;
    }
    for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
    {
        (*it) = true;
    }

    // flip the active likelihood pointers
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        if ( changed_nodes[index] == false )
        {
            active_node_likelihood[index] = (active_node_likelihood[index] == 0 ? 1 : 0);
            changed_nodes[index] = true;
        }
    }
    for (size_t index = 0; index < changed_branches.size(); ++index)
    {
        if ( changed_branches[index] == false )
        {
            active_branch_likelihood[index] = (active_branch_likelihood[index] == 0 ? 1 : 0);
            changed_branches[index] = true;
        }
    }
    
    for (std::vector<bool>::iterator it = pmat_dirty_branches.begin(); it != pmat_dirty_branches.end(); ++it)
    {
        (*it) = true;
    }
    
    for (size_t index = 0; index < pmat_changed_nodes.size(); ++index)
    {
        if ( pmat_changed_nodes[index] == false )
        {
            active_pmatrices[index] = (active_pmatrices[index] == 0 ? 1 : 0);
            pmat_changed_nodes[index] = true;
        }
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::reInitialized( void )
{
    // we need to recompress because the tree may have changed
    compress();

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::resizeLikelihoodVectors( void )
{
    if (mixture_model)
    {
        num_site_mixtures = mixture_model->getValue().getNumberOfComponents();
        // FIXME - not sure if this is right for covarion models.
        assert(num_states == mixture_model->getValue().getNumberOfStates());
    }
    else if (this->branch_heterogeneous_substitution_matrices == false)
    {
        this->num_site_mixtures = this->num_site_rates * this->num_matrices;
    }
    else
    {
        this->num_site_mixtures = this->num_site_rates;
    }

    // set the offsets for easier iteration through the likelihood vector
    site_offset                         =  num_states;
    mixture_offset                      =  pattern_block_size*site_offset;
    node_offset                         =  num_site_mixtures*mixture_offset;
    tip_offset                          =  pattern_block_size*site_offset;
    active_branch_likelihood_offset     =  num_nodes*node_offset;
    active_node_likelihood_offset       =  num_tips*node_offset;

    // only do this if we are in MCMC mode. This will safe memory
    if ( in_mcmc_mode == true )
    {
        // we resize the partial likelihood vectors to the new dimensions
        delete [] partial_branch_likelihoods;
        delete [] partial_node_likelihoods;

        if ( RbSettings::userSettings().getPartialLikelihoodStoring() == "branch" || RbSettings::userSettings().getPartialLikelihoodStoring() == "both" )
        {
            partial_branch_likelihoods = new double[2*active_branch_likelihood_offset];

            // reinitialize likelihood vectors
            for (size_t i = 0; i < 2*active_branch_likelihood_offset; i++)
            {
                partial_branch_likelihoods[i] = 0.0;
            }
        }

        if ( RbSettings::userSettings().getPartialLikelihoodStoring() == "node" || RbSettings::userSettings().getPartialLikelihoodStoring() == "both" )
        {
            partial_node_likelihoods = new double[2*active_node_likelihood_offset];

            // reinitialize likelihood vectors
            for (size_t i = 0; i < 2*active_node_likelihood_offset; i++)
            {
                partial_node_likelihoods[i] = 0.0;
            }
        }

        for(auto&& dirty_node: dirty_nodes)
        {
            dirty_node = true;
        }
        for(auto&& dirty_branch: dirty_branches)
        {
            dirty_branch = true;
        }
    }

    if ( useMarginalLikelihoods == true )
    {
        // we resize the partial likelihood vectors to the new dimensions
        delete [] marginalLikelihoods;

        marginalLikelihoods = new double[active_branch_likelihood_offset];

        // reinitialize likelihood vectors
        for (size_t i = 0; i < active_branch_likelihood_offset; i++)
        {
            marginalLikelihoods[i] = 0.0;
        }

    }

    per_node_site_mixture_log_scaling_factors   = std::vector<std::vector< std::vector< std::vector<double> > > >(2, std::vector<std::vector< std::vector<double> > >(num_nodes, std::vector< std::vector<double> >(num_site_mixtures, std::vector<double>(pattern_block_size, 0.0) ) ) );
    per_node_site_log_scaling_factors           = std::vector<std::vector< std::vector< double > > >(2, std::vector<std::vector< double > >(num_nodes, std::vector<double>(pattern_block_size, 0.0 ) ) );

    active_P_matrix_offset  =  num_nodes * num_site_mixtures;
    pmat_node_offset        =  num_site_mixtures;
    pmatrices               =  std::vector<TransitionProbabilityMatrix>(active_P_matrix_offset * 2, TransitionProbabilityMatrix(num_states));

    for (auto&& pmat_dirty_branch: pmat_dirty_branches)
    {
        pmat_dirty_branch = true;
    }

    transition_prob_matrices = std::vector<TransitionProbabilityMatrix>(num_site_mixtures, TransitionProbabilityMatrix(num_states) );
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::restoreSpecialization( const DagNode* affecter )
{
    // reset flags for likelihood computation
    touched = false;

    // reset the flags
    lnProb = *storedLnProb;
    storedLnProb = {};
    dirty_nodes = *prev_dirty_nodes;
    dirty_branches = *prev_dirty_branches;
    prev_dirty_nodes = {};
    prev_dirty_branches = {};

    // restore the active likelihoods vector
    for (size_t index = 0; index < changed_nodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changed_nodes[index] == true )
            active_node_likelihood[index] = (active_node_likelihood[index] == 0 ? 1 : 0);

        // set all flags to false
        changed_nodes[index] = false;
    }

    for (size_t index = 0; index < changed_branches.size(); ++index)
    {
        // we have to restore, that means if we have changed the active likelihood vector
        // then we need to revert this change
        if ( changed_branches[index] == true )
        {
            active_branch_likelihood[index] = (active_branch_likelihood[index] == 0 ? 1 : 0);
        }

        // set all flags to false
        changed_branches[index] = false;
    }

    // reset the flags
    pmat_dirty_branches = *prev_pmat_dirty_branches;
    prev_pmat_dirty_branches = {};

    // restore the active transition probability matrices vector
    for (size_t index = 0; index < pmat_changed_nodes.size(); ++index)
    {
        // we have to restore, that means if we have changed the active transition probability matrices vector
        // then we need to revert this change
        if ( pmat_changed_nodes[index] == true )
            active_pmatrices[index] = (active_pmatrices[index] == 0 ? 1 : 0);

        // set all flags to false
        pmat_changed_nodes[index] = false;
    }
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::scale( size_t node_index)
{
    
    bool test_underflow     = RbSettings::userSettings().getUseScaling() == true;
    bool test_this_node     = ((node_index+1) % RbSettings::userSettings().getScalingDensity() == 0);
    bool scale_threshold    = RbSettings::userSettings().getScalingMethod() == "threshold";
    bool scale_per_mixture  = RbSettings::userSettings().getScalingPerMixture();
    
    double* p_node = NULL;
    if ( partial_likelihood_storing_approach == BRANCH )
    {
        p_node = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;
    }
    else
    {
        p_node = this->partial_node_likelihoods + this->active_node_likelihood[node_index]*this->active_node_likelihood_offset + (node_index-this->num_tips)*this->node_offset;
    }
    
    if ( test_underflow == true && scale_per_mixture == false )
    {
        if ( test_this_node == true )
        {
            // iterate over all mixture categories
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // the max probability
                double max = 0.0;

                // compute the per site probabilities
                for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                {
                    // get the pointers to the likelihood for this mixture category
                    size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                    double* p_site_mixture = p_node + offset;

                    for ( size_t i=0; i<this->num_states; ++i)
                    {
                        if ( p_site_mixture[i] > max )
                        {
                            max = p_site_mixture[i];
                        }
                    }

                }
                
                // should we rescale by max-prob or by a fixed threshold
                if ( scale_threshold == false )
                {

                    // Don't divide by zero or NaN.
                    if (not (max > 0)) continue;

                    this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = -log(max);

                    // compute the per site probabilities
                    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                    {
                        // get the pointers to the likelihood for this mixture category
                        size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                        double* p_site_mixture = p_node + offset;

                        for ( size_t i=0; i<this->num_states; ++i)
                        {
                            p_site_mixture[i] /= max;
                        }
                    }
                    
                } // end-if scale by max prob
                else
                {
                    // only scale if the max-prob was smaller than a fixed threshold
                    if ( max < RbConstants::SCALING_THRESHOLD )
                    {
                        this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = 1;

                        // compute the per site probabilities
                        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                        {
                            // get the pointers to the likelihood for this mixture category
                            size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                            double* p_site_mixture = p_node + offset;

                            for ( size_t i=0; i<this->num_states; ++i)
                            {
                                p_site_mixture[i] /= RbConstants::SCALING_THRESHOLD;
                            }

                        }
                        
                    } // end-if the max prob was smaller than a fixed threshold
                    else
                    {
                        this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = 0;
                    }
                } // end-else scale by fixed threshold

            }
        } // end-if testing this node
        else
        {
            // iterate over all sites
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {
                this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = 0;
            }
        } // end-else not testing this node
    }
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::scale( size_t node_index, size_t left, size_t right )
{
    
    double* p_node = NULL;
    if ( partial_likelihood_storing_approach == BRANCH )
    {
        p_node = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;
    }
    else
    {
        p_node = this->partial_node_likelihoods + this->active_node_likelihood[node_index]*this->active_node_likelihood_offset + (node_index-this->num_tips)*this->node_offset;
    }
    
    bool test_underflow     = RbSettings::userSettings().getUseScaling() == true;
    bool test_this_node     = ((node_index+1) % RbSettings::userSettings().getScalingDensity() == 0);
    bool scale_threshold    = RbSettings::userSettings().getScalingMethod() == "threshold";
    bool scale_per_mixture  = RbSettings::userSettings().getScalingPerMixture();
    
    if ( test_underflow == true && scale_per_mixture == false )
    {
        
        // should we test this node for underflow?
        if ( test_this_node == true )
        {

            // iterate over all mixture categories
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // the max probability
                double max = 0.0;

                // compute the per site probabilities
                for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                {
                    // get the pointers to the likelihood for this mixture category
                    size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                    double*          p_site_mixture          = p_node + offset;

                    for ( size_t i=0; i<this->num_states; ++i)
                    {
                        if ( p_site_mixture[i] > max )
                        {
                            max = p_site_mixture[i];
                        }

                    }

                }

                // should we rescale by max-prob or by a fixed threshold
                if ( scale_threshold == false )
                {
                    // Don't divide by zero or NaN.
                    if (not (max > 0)) continue;

                    this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site] - log(max);

                    // compute the per site probabilities
                    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                    {
                        // get the pointers to the likelihood for this mixture category
                        size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                        double* p_site_mixture = p_node + offset;

                        for ( size_t i=0; i<this->num_states; ++i)
                        {
                            p_site_mixture[i] /= max;
                        }
                    }
                } // end-if scale by max prob
                else
                {
                    // only scale if the max-prob was smaller than a fixed threshold
                    if ( max < RbConstants::SCALING_THRESHOLD )
                    {
                        this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site] + 1;

                        // compute the per site probabilities
                        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                        {
                            // get the pointers to the likelihood for this mixture category
                            size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                            double* p_site_mixture = p_node + offset;

                            for ( size_t i=0; i<this->num_states; ++i)
                            {
                                p_site_mixture[i] /= RbConstants::SCALING_THRESHOLD;
                            }

                        }
                        
                    } // end-if the max prob was smaller than a fixed threshold
                    else
                    {
                        this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site];
                    }
                } // end-else scale by fixed threshold

            }

        } // end-if we should test this node for underflow
        else
        {
            // iterate over all mixture categories
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {
                this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site];
            }

        }
        
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::scale( size_t node_index, size_t left, size_t right, size_t middle )
{
    
    double* p_node = NULL;
    if ( partial_likelihood_storing_approach == BRANCH )
    {
        p_node = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index]*this->active_branch_likelihood_offset + node_index*this->node_offset;
    }
    else
    {
        p_node = this->partial_node_likelihoods + this->active_node_likelihood[node_index]*this->active_node_likelihood_offset + (node_index-this->num_tips)*this->node_offset;
    }

    
    bool test_underflow     = RbSettings::userSettings().getUseScaling() == true;
    bool test_this_node     = ((node_index+1) % RbSettings::userSettings().getScalingDensity() == 0);
    bool scale_threshold    = RbSettings::userSettings().getScalingMethod() == "threshold";
    bool scale_per_mixture  = RbSettings::userSettings().getScalingPerMixture();
    
    if ( test_underflow == true && scale_per_mixture == false )
    {
        
        // should we test this node for underflow?
        if ( test_this_node == true )
        {
            
            // iterate over all mixture categories
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {

                // the max probability
                double max = 0.0;

                // compute the per site probabilities
                for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                {
                    // get the pointers to the likelihood for this mixture category
                    size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                    double* p_site_mixture = p_node + offset;

                    for ( size_t i=0; i<this->num_states; ++i)
                    {
                        if ( p_site_mixture[i] > max )
                        {
                            max = p_site_mixture[i];
                        }
                    }

                }

                // should we rescale by max-prob or by a fixed threshold
                if ( scale_threshold == false )
                {
                    // Don't divide by zero or NaN.
                    if (not (max > 0)) continue;

                    this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[middle]][middle][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site] - log(max);

                    // compute the per site probabilities
                    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                    {
                        // get the pointers to the likelihood for this mixture category
                        size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                        double* p_site_mixture = p_node + offset;

                        for ( size_t i=0; i<this->num_states; ++i)
                        {
                            p_site_mixture[i] /= max;
                        }
                    }
                } // end-if scale by max prob
                else
                {
                    // only scale if the max-prob was smaller than a fixed threshold
                    if ( max < RbConstants::SCALING_THRESHOLD )
                    {
                        this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[middle]][middle][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site] + 1;

                        // compute the per site probabilities
                        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
                        {
                            // get the pointers to the likelihood for this mixture category
                            size_t offset = mixture*this->mixture_offset + site*this->site_offset;

                            double* p_site_mixture = p_node + offset;

                            for ( size_t i=0; i<this->num_states; ++i)
                            {
                                p_site_mixture[i] /= RbConstants::SCALING_THRESHOLD;
                            }

                        }
                        
                    } // end-if the max prob was smaller than a fixed threshold
                    else
                    {
                        this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[middle]][middle][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site];
                    }
                } // end-else scale by fixed threshold

            }

        } // end-if we should test this node for underflow
        else
        {
            // iterate over all mixture categories
            for (size_t site = 0; site < this->pattern_block_size ; ++site)
            {
                this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[middle]][middle][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[left]][left][site] + this->per_node_site_log_scaling_factors[this->active_branch_likelihood[right]][right][site];
            }

        }
        
    }

}



template <class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setActivePIDSpecialized(size_t a, size_t n)
{
    // we need to recompress the data
    this->compress();
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setValue(AbstractHomologousDiscreteCharacterData *v, bool force)
{
    if ( v->getMaxObservedStateIndex() > this->num_states - 1)
    {
        // We might use different sized matrices for different partitions depending on the observed number of states.
        std::stringstream ss;
        ss << "The number of observed states (" << v->getMaxObservedStateIndex() + 1 << ") is greater than the dimension of the Q matrix (" << this->num_states << ")" << std::endl;
        throw RbException(ss.str());
    }

    if (v->getDataType() != this->template_state.getDataType() )
    {
      // There is a mismatch between the data type of the data matrix and the data type of the CTMC.
      std::stringstream ss;
      ss << "The data type of the data matrix ("<< v->getDataType() <<") differs from that of the PhyloctMC object ("<< this->template_state.getDataType() <<")." << std::endl;
      throw RbException(ss.str());
    }

    // delegate to the parent class
    TypedDistribution< AbstractHomologousDiscreteCharacterData >::setValue(v, force);

    // reset the number of sites
    this->num_sites = v->getNumberOfIncludedCharacters();

    site_pattern.clear();
    site_pattern.resize(num_sites);

    // now compress the data and resize the likelihood vectors
    this->compress();
    
    // now we also set the template state
    template_state = charType( static_cast<const charType&>( this->value->getTaxonData(0).getCharacter(0) ) );
    template_state.setToFirstState();
    template_state.setGapState( false );
    template_state.setMissingState( false );
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::simulate( const TopologyNode &node, std::vector< DiscreteTaxonData< charType > > &taxa, const std::vector<bool> &invariant, const std::vector<size_t> &perSiteMixtures)
{

    // get the children of the node
    const std::vector<TopologyNode*>& children = node.getChildren();

    // get the sequence of this node
    size_t node_index = node.getIndex();
    const DiscreteTaxonData< charType > &parent = taxa[ node_index ];

    // simulate the sequence for each child
    RandomNumberGenerator* rng = GLOBAL_RNG;
    for (std::vector< TopologyNode* >::const_iterator it = children.begin(); it != children.end(); ++it)
    {
        const TopologyNode &child = *(*it);

        // update the transition probability matrix
        updateTransitionProbabilities( child.getIndex() );

        DiscreteTaxonData< charType > &taxon = taxa[ child.getIndex() ];
        for ( size_t i = 0; i < num_sites; ++i )
        {

            if ( invariant[i] == true )
            {

                // add the character to the sequence
                taxon.addCharacter( parent.getCharacter( i ) );
            }
            else
            {
                // get the ancestral character for this site
                std::uint64_t parentState = parent.getCharacter( i ).getStateIndex();

                double *freqs = transition_prob_matrices[ perSiteMixtures[i] ][ parentState ];

                // create the character
                charType c = charType( num_states );
                c.setToFirstState();
                // draw the state
                double u = rng->uniform01();
                size_t stateIndex = 0;
                while ( true )
                {
                    u -= *freqs;
                    ++stateIndex;

                    if ( u > 0.0 && stateIndex < this->num_states)
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

        }

        if ( child.isTip() )
        {
            taxon.setTaxon( child.getTaxon() );
        }
        else
        {
            // recursively simulate the sequences
            std::stringstream ss;
            ss << "Node" << child.getIndex();
            taxon.setTaxon( Taxon(ss.str()) );
            simulate( child, taxa, invariant, perSiteMixtures );
        }

    }

}


template<class charType>
const RevBayesCore::TypedDagNode<RevBayesCore::Tree>* RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::getTree()
{
    return tau;
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setClockRate(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( homogeneous_clock_rate != NULL )
    {
        this->removeParameter( homogeneous_clock_rate );
        homogeneous_clock_rate = NULL;
    }
    else // heterogeneousClockRate != NULL
    {
        this->removeParameter( heterogeneous_clock_rates );
        heterogeneous_clock_rates = NULL;
    }

    // set the value
    branch_heterogeneous_clock_rates = false;
    homogeneous_clock_rate = r;

    // add the new parameter
    this->addParameter( homogeneous_clock_rate );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setClockRate(const TypedDagNode< RbVector< double > > *r)
{

    // remove the old parameter first
    if ( homogeneous_clock_rate != NULL )
    {
        this->removeParameter( homogeneous_clock_rate );
        homogeneous_clock_rate = NULL;
    }
    else // heterogeneousClockRate != NULL
    {
        this->removeParameter( heterogeneous_clock_rates );
        heterogeneous_clock_rates = NULL;
    }

    // set the value
    branch_heterogeneous_clock_rates = true;
    heterogeneous_clock_rates = r;

    // add the new parameter
    this->addParameter( heterogeneous_clock_rates );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setMcmcMode(bool tf)
{
    // free old memory
    if ( in_mcmc_mode == true )
    {
#if defined( RB_BEAGLE )
        //-- If there is already a BEAGLE instance, delete it.
        charType tmp_char;
        if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
        {
            this->freeBeagleInstances();
        }
        else
        {
            delete [] partial_branch_likelihoods;
            delete [] partial_node_likelihoods;
            partial_branch_likelihoods = NULL;
            partial_node_likelihoods   = NULL;
        }
#else
        delete [] partial_branch_likelihoods;
        delete [] partial_node_likelihoods;
        partial_branch_likelihoods = NULL;
        partial_node_likelihoods   = NULL;
#endif
    }

    // set our internal flag
    in_mcmc_mode = tf;

    // resize likelihood buffers
    if ( in_mcmc_mode == true )
    {
#if defined( RB_BEAGLE )
        // Initialize a new BEAGLE instance for MCMC
        charType tmp_char;
        if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
        {
            this->initializeBeagleInstances();
        }
        else
        {
            resizeLikelihoodVectors();
            fillTipLikelihoods();
        }
#else
        resizeLikelihoodVectors();
        fillTipLikelihoods();
#endif
    }
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setPInv(const TypedDagNode< double > *r)
{

    // remove the old parameter first
    if ( p_inv != NULL )
    {
        this->removeParameter( p_inv );
        p_inv = NULL;
    }

    // set the value
    p_inv = r;

    // add the new parameter
    this->addParameter( p_inv );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setRateMatrix(const TypedDagNode< RateGenerator > *rm)
{

    // remove the old parameter first
    if ( homogeneous_rate_matrix != NULL )
    {
        this->removeParameter( homogeneous_rate_matrix );
        homogeneous_rate_matrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( heterogeneous_rate_matrices );
        heterogeneous_rate_matrices = NULL;
    }

    // set the value
    homogeneous_rate_matrix = rm;
    num_matrices = 1;

    this->resizeLikelihoodVectors();

    if (rm != NULL && rm->getValue().size() != num_states)
        throw RbException()<<"Rate generator dimensions (" << rm->getValue().size() << " do not match the number of character states (" << this->num_states << ")";

    // add the new parameter
    this->addParameter( homogeneous_rate_matrix );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setMixtureModel(const TypedDagNode< SiteMixtureModel> *mm)
{

    // remove the old parameter first
    if ( homogeneous_rate_matrix != NULL)
    {
        removeParameter( homogeneous_rate_matrix );
        homogeneous_rate_matrix = NULL;
    }

    if ( heterogeneous_rate_matrices != NULL)
    {
        removeParameter( heterogeneous_rate_matrices );
        heterogeneous_rate_matrices = NULL;
    }

    // set the value
    mixture_model = mm;
    num_matrices = 1;

    resizeLikelihoodVectors();

    if (mm && mm->getValue().getNumberOfStates() != this->num_states)
    {
        throw RbException()<<"Mixture model (" << mm->getValue().getNumberOfStates() << " do not match the number of character states (" << num_states << ")";
    }

    // add the new parameter
    addParameter( mixture_model );

    // redraw the current value
    if ( not dag_node or not dag_node->isClamped())
    {
        redrawValue();
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setRateMatrix(const TypedDagNode< RbVector< RateGenerator > > *rm)
{

    // remove the old parameter first
    if ( homogeneous_rate_matrix != NULL )
    {
        this->removeParameter( homogeneous_rate_matrix );
        homogeneous_rate_matrix = NULL;
    }
    else // heterogeneousRateMatrix != NULL
    {
        this->removeParameter( heterogeneous_rate_matrices );
        heterogeneous_rate_matrices = NULL;
    }

    // set the value
    heterogeneous_rate_matrices = rm;
    num_matrices = rm == NULL ? 1 : rm->getValue().size();

    this->resizeLikelihoodVectors();

    if (rm != NULL && rm->getValue()[0].size() != this->num_states)
    {
        std::stringstream ss;
        ss << "Rate generator dimensions (" << rm->getValue()[0].size() << " do not match the number of character states (" << this->num_states << ")";
        throw(RbException(ss.str()));
    }

    // add the new parameter
    this->addParameter( heterogeneous_rate_matrices );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setRootFrequencies(const TypedDagNode< Simplex > *f)
{

    // remove the old parameter first
    if ( root_frequencies != NULL )
    {
        this->removeParameter( root_frequencies );
        root_frequencies = NULL;
    }

    if ( f != NULL )
    {
        // set the value
        root_frequencies = f;
    }

    // add the new parameter
    this->addParameter( root_frequencies );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setSiteRates(const TypedDagNode< RbVector< double > > *r)
{

    // remove the old parameter first
    if ( site_rates != NULL )
    {
        this->removeParameter( site_rates );
        site_rates = NULL;
    }

    if ( r != NULL )
    {
        // set the value
        rate_variation_across_sites = true;
        site_rates = r;
        this->num_site_rates = r->getValue().size();
    }
    else
    {
        // set the value
        rate_variation_across_sites = false;
        site_rates = NULL;
        this->num_site_rates = 1;
    }

    this->resizeLikelihoodVectors();

    // add the new parameter
    this->addParameter( site_rates );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setSiteRatesProbs(const TypedDagNode< Simplex > *rp)
{

    // remove the old parameter first
    if ( site_rates_probs != NULL )
    {
        this->removeParameter( site_rates_probs );
        site_rates_probs = NULL;
    }

    if (rp != NULL)
    {
        // set the value
        site_rates_probs = rp;
    }

    // add the new parameter
    this->addParameter( site_rates_probs );

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setUseMarginalLikelihoods(bool tf)
{
    this->useMarginalLikelihoods = tf;
    this->resizeLikelihoodVectors();
}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::setUseSiteMatrices(bool use_sm, const TypedDagNode< Simplex > *s)
{

    if ( use_sm == false && s != NULL)
    {
        throw(RbException("Provided site matrix probs but not using site matrix mixture."));
    }

    // remove the old parameter first
    if ( site_matrix_probs != NULL )
    {
        this->removeParameter( site_matrix_probs );
        site_matrix_probs = NULL;
    }

    if ( use_sm == true )
    {
        // set the value
        site_matrix_probs = s;
    }

    // add the new parameter
    this->addParameter( site_matrix_probs );

    this->branch_heterogeneous_substitution_matrices = !use_sm;

    this->resizeLikelihoodVectors();

    // redraw the current value
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false )
    {
        this->redrawValue();
    }
}


template<class charType>
std::vector< std::vector<double> >* RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::sumMarginalLikelihoods( size_t node_index )
{
    std::vector< std::vector<double> >* per_mixture_Likelihoods = new std::vector< std::vector<double> >(this->pattern_block_size, std::vector<double>(num_states, 0.0) );

    std::vector<double> mixture_probs = getMixtureProbs();

    // get the pointers to the partial likelihoods and the marginal likelihoods
    double*         p_node_marginal         = this->marginalLikelihoods + node_index*this->node_offset;

    // get pointers the likelihood for both subtrees
    double*         p_mixture_marginal          = p_node_marginal;
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {

        // get pointers to the likelihood for this mixture category
        double*         p_site_mixture_marginal         = p_mixture_marginal;
        // iterate over all sites
        for (size_t site = 0; site < this->pattern_block_size; ++site)
        {
            // get the pointers to the likelihoods for this site and mixture category
            double*         p_site_marginal_j           = p_site_mixture_marginal;
            // iterate over all starting states
            for (size_t j=0; j<num_states; ++j)
            {
                // add the probability of being in this state
                (*per_mixture_Likelihoods)[site][j] += *p_site_marginal_j * mixture_probs[mixture];

                // increment pointers
                ++p_site_marginal_j;
            }

            // increment the pointers to the next site
            p_site_mixture_marginal+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture_marginal+=this->mixture_offset;

    } // end-for over all mixtures (=rate categories)

    return per_mixture_Likelihoods;
}




template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoods( std::vector<double> &rv ) const
{
    // get the root node
    const TopologyNode &root = tau->getValue().getRoot();

    // get the index of the root node
    size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods of the left and right subtree
    double*   p_node  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index] * this->active_branch_likelihood_offset  + node_index*this->node_offset;
    if ( RbSettings::userSettings().getPartialLikelihoodStoring() == "node" || RbSettings::userSettings().getPartialLikelihoodStoring() == "both" )
    {
        p_node  = this->partial_node_likelihoods + this->active_node_likelihood[node_index] * this->active_node_likelihood_offset  + (node_index-this->num_tips)*this->node_offset;
    }

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<double> per_mixture_Likelihoods = std::vector<double>(pattern_block_size,0.0);
    std::vector<double> per_mixture_scaling_factors = std::vector<double>(pattern_block_size,0.0);

    std::vector<double> site_mixture_probs = getMixtureProbs();

    // get pointer the likelihood
    double*   p_mixture     = p_node;
    
    bool use_scaling        = RbSettings::userSettings().getUseScaling();
    bool scale_threshold    = RbSettings::userSettings().getScalingMethod() == "threshold";
    bool scale_per_mixture  = RbSettings::userSettings().getScalingPerMixture();


    // iterate over all sites
    for (size_t site = 0; site < pattern_block_size; ++site)
    {
        
        if ( use_scaling == true )
        {
            if ( scale_per_mixture == true )
            {
                // first find the max scaling
                double max_scaling = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site];
                // iterate over all mixture categories
                for (size_t mixture = 1; mixture < this->num_site_mixtures; ++mixture)
                {
                    double tmp_sf = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site];
                    if ( scale_threshold == true )
                    {
                        max_scaling = (max_scaling < tmp_sf ? max_scaling : tmp_sf);
                    }
                    else
                    {
                        max_scaling = (max_scaling < tmp_sf ? max_scaling : tmp_sf);
                    }
                }
                per_mixture_scaling_factors[site] = max_scaling;
            }
            else
            {
                per_mixture_scaling_factors[site] = this->per_node_site_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][site];
            }
        }
        
        // iterate over all mixture categories
        for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
        {

            // get pointers to the likelihood for this mixture category
            double*   p_site_mixture     = p_mixture + (mixture*this->mixture_offset) + (site*this->site_offset);

            // temporary variable storing the likelihood
            double tmp = 0.0;
            // get the pointers to the likelihoods for this site and mixture category
            double* p_site_j   = p_site_mixture;
            // iterate over all starting states
            for (size_t i=0; i<num_states; ++i)
            {
                // add the probability of starting from this state
                tmp += *p_site_j;

                // increment pointers
                ++p_site_j;
            }
            
            if ( use_scaling == true )
            {
                if ( scale_per_mixture == true )
                {
                    // add the likelihood for this mixture category
                    double tmp_sf = this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][mixture][site];
                    if ( scale_threshold == false )
                    {
                        double test_nan = exp(per_mixture_scaling_factors[site] - tmp_sf) * site_mixture_probs[mixture] * tmp;
                        if ( RbMath::isFinite(test_nan) == false )
                        {
                            std::cerr << "exp(" << per_mixture_scaling_factors[site] << "-" << tmp_sf<< ") * site_mixture_probs[mixture] *" << tmp << std::endl;
                        }
                        per_mixture_Likelihoods[site] += exp(per_mixture_scaling_factors[site] - tmp_sf) * site_mixture_probs[mixture] * tmp;
                    }
                    else
                    {
                        per_mixture_Likelihoods[site] += pow(2, 256*(per_mixture_scaling_factors[site] - tmp_sf)) * site_mixture_probs[mixture] * tmp;
                    }
                }
                else
                {
                    per_mixture_Likelihoods[site] += site_mixture_probs[mixture] * tmp;
                }
            }
            else
            {
                per_mixture_Likelihoods[site] += site_mixture_probs[mixture] * tmp;
            }
            // increment the pointers to the next site
//            p_site_mixture+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
//        p_mixture+=this->mixture_offset;

    } // end-for over all mixtures

    double prob_invariant = getPInv();
    
    double oneMinusPInv = 1.0 - prob_invariant;
    std::vector< size_t >::const_iterator patterns = this->pattern_counts.begin();
    if ( prob_invariant > 0.0 )
    {
        // get the mean root frequency vector
        std::vector<double> f;
        if (this->branch_heterogeneous_substitution_matrices == true)
        {
            f = this->getRootFrequencies(0);
        }
        else
        {
            std::vector<std::vector<double> > ff;
            getRootFrequencies(ff);

            std::vector<double> matrix_probs(num_matrices, 1.0/num_matrices);

            if (site_matrix_probs != NULL)
            {
                matrix_probs = site_matrix_probs->getValue();
            }

            f = std::vector<double>(ff[0].size(), 0.0);

            for (size_t matrix = 0; matrix < ff.size(); matrix++)
            {
                // get the root frequencies
                const std::vector<double> &fm = ff[matrix];

                for (size_t i = 0; i < fm.size(); i++)
                {
                    f[i] += fm[i] * matrix_probs[matrix];
                }
            }
        }

        for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
        {

           
            if ( use_scaling == true )
            {
                
                // get the scaling factor for this site
                double this_sf = per_mixture_scaling_factors[site];
                
                if ( this->site_invariant[site] == true )
                {
                    double ftotal = 0.0;
                    for ( size_t c = 0; c < this->invariant_site_index[site].size(); c++ )
                    {
                        ftotal += f[this->invariant_site_index[site][c]];
                    }

                    if ( scale_threshold == false )
                    {
                        rv[site] = log( prob_invariant * ftotal + oneMinusPInv * per_mixture_Likelihoods[site] * exp(-this_sf) ) * *patterns;
                    }
                    else
                    {
                        rv[site] = log( prob_invariant * ftotal + oneMinusPInv * per_mixture_Likelihoods[site] * pow(RbConstants::SCALING_THRESHOLD, this_sf) ) * *patterns;
                    }
                }
                else
                {
                    rv[site] = log( oneMinusPInv * per_mixture_Likelihoods[site] ) * *patterns;
                    if ( scale_threshold == false )
                    {
                        rv[site] -= this_sf * *patterns;
                    }
                    else
                    {
                        rv[site] -= (RbConstants::LN2 * 256 * this_sf) * *patterns;
                    }

                }

            }
            else // no scaling
            {

                if ( this->site_invariant[site] == true )
                {
                    double ftotal = 0.0;
                    for ( size_t c = 0; c < this->invariant_site_index[site].size(); c++ )
                    {
                        ftotal += f[this->invariant_site_index[site][c]];
                    }
                    
                    rv[site] = log( prob_invariant * ftotal  + oneMinusPInv * per_mixture_Likelihoods[site] ) * *patterns;
                }
                else
                {
                    rv[site] = log( oneMinusPInv * per_mixture_Likelihoods[site] ) * *patterns;
                }

            }

        }
        
    }
    else
    {

        for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
        {
            rv[site] = log( per_mixture_Likelihoods[site] ) * *patterns;

            if ( use_scaling == true )
            {
                // get the scaling factor for this site
                double this_sf = per_mixture_scaling_factors[site];
                
                if ( scale_threshold == false )
                {
                    rv[site] -= this_sf * *patterns;
                }
                else
                {
                    rv[site] -= (RbConstants::LN2 * 256 * this_sf) * *patterns;
                }
            }

        }

    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoodsPerSiteMixture( MatrixReal &rv ) const
{
    // get the root node
    const TopologyNode &root = tau->getValue().getRoot();

    // get the index of the root node
    size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods of the left and right subtree
    double*   p_node  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index] * this->active_branch_likelihood_offset  + node_index*this->node_offset;

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<std::vector<double> > per_site_mixture_Likelihoods = std::vector<std::vector<double> >(pattern_block_size, std::vector<double>(num_site_mixtures, 0.0));

    std::vector<double> site_mixture_probs = getMixtureProbs();
    
    bool use_scaling     = RbSettings::userSettings().getUseScaling();
    bool scale_threshold = RbSettings::userSettings().getScalingMethod() == "threshold";

    // get pointer the likelihood
    double*   p_mixture     = p_node;
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {

        // get pointers to the likelihood for this mixture category
        double*   p_site_mixture     = p_mixture;
        // iterate over all sites

        for (size_t site = 0; site < pattern_block_size; ++site)
        {
            // temporary variable storing the likelihood
            double tmp = 0.0;
            // get the pointers to the likelihoods for this site and mixture category
            double* p_site_j   = p_site_mixture;
            // iterate over all starting states
            for (size_t i=0; i<num_states; ++i)
            {
                // add the probability of starting from this state
                tmp += *p_site_j;

                // increment pointers
                ++p_site_j;
            }
            // add the likelihood for this mixture category
            per_site_mixture_Likelihoods[site][mixture] += tmp * site_mixture_probs[mixture];

            // increment the pointers to the next site
            p_site_mixture+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixture_offset;

    } // end-for over all mixtures

    double prob_invariant = getPInv();
    double oneMinusPInv = 1.0 - prob_invariant;
    std::vector< size_t >::const_iterator patterns = this->pattern_counts.begin();
    if ( prob_invariant > 0.0 )
    {
        // get the root frequency vector(s)
        std::vector<std::vector<double> > ff;
        std::vector<double> f;
        if (this->branch_heterogeneous_substitution_matrices == true)
        {
            f = this->getRootFrequencies(0);
            ff.push_back(f);
        }
        else
        {
            getRootFrequencies(ff);
        }

        size_t num_site_matrices = num_site_mixtures / num_site_rates;
        std::vector<double> matrix_probs(num_site_matrices, 1.0/num_site_matrices);

        if (site_matrix_probs != NULL)
        {
            matrix_probs = site_matrix_probs->getValue();
        }

        for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
        {
            for (size_t matrix = 0; matrix < num_site_matrices; ++matrix)
            {
                // the first rate category is the invariant
                if ( this->site_invariant[site] == true )
                {
                    double ftotal = 0.0;
                    for ( size_t c = 0; c < this->invariant_site_index[site].size(); c++ )
                    {
                        ftotal += f[this->invariant_site_index[site][c]];
                    }

                    rv[site][matrix] = log( prob_invariant * ftotal * matrix_probs[matrix] ) * *patterns;
                }
                else
                {
                    rv[site][matrix] = RbConstants::Double::neginf;
                }

                // the remaining variant rate categories
                for (size_t site_rate_index = 1; site_rate_index < num_site_rates + 1; ++site_rate_index)
                {
                    rv[site][site_rate_index * num_site_matrices + matrix] = log( oneMinusPInv * per_site_mixture_Likelihoods[site][site_rate_index * num_site_matrices + matrix] ) * *patterns;

                    if ( use_scaling == true )
                    {
                        if ( scale_threshold == false )
                        {
                            rv[site][site_rate_index * num_site_matrices + matrix] -= this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site] * *patterns;
                        }
                        else
                        {
                            rv[site][site_rate_index * num_site_matrices + matrix] -= (RbConstants::LN2 * 256 * this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site]) * *patterns;
                        }
                            
                    }

                }

            }

        }

    }
    else
    {

        for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
        {
            for (size_t mixture = 0; mixture < num_site_mixtures; ++mixture)
            {
                rv[site][mixture] = log( per_site_mixture_Likelihoods[site][mixture] ) * *patterns;

                if ( use_scaling == true )
                {
                    if ( scale_threshold == false )
                    {
                        rv[site][mixture] -= this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site] * *patterns;
                    }
                    else
                    {
                        rv[site][mixture] -= (RbConstants::LN2 * 256 * this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site]) * *patterns;
                    }
                }
            }

        }

    }

}


template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::computeRootLikelihoodsPerSiteRate( MatrixReal &rv ) const
{
    // get the root node
    const TopologyNode &root = tau->getValue().getRoot();

    // get the index of the root node
    size_t node_index = root.getIndex();

    // get the pointers to the partial likelihoods of the left and right subtree
    double*   p_node  = this->partial_branch_likelihoods + this->active_branch_likelihood[node_index] * this->active_branch_likelihood_offset  + node_index*this->node_offset;

    size_t num_site_matrices = num_site_mixtures/num_site_rates;

    // create a vector for the per mixture likelihoods
    // we need this vector to sum over the different mixture likelihoods
    std::vector<std::vector<double> > per_site_rate_Likelihoods = std::vector<std::vector<double> >(pattern_block_size, std::vector<double>(num_site_rates, 0.0));

    std::vector<double> site_mixture_probs = getMixtureProbs();
    
    
    bool use_scaling     = RbSettings::userSettings().getUseScaling();
    bool scale_threshold = RbSettings::userSettings().getScalingMethod() == "threshold";
    
    // get pointer the likelihood
    double*   p_mixture     = p_node;
    // iterate over all mixture categories
    for (size_t mixture = 0; mixture < this->num_site_mixtures; ++mixture)
    {
        size_t site_rate_index = mixture / num_site_matrices;

        // get pointers to the likelihood for this mixture category
        double*   p_site_mixture     = p_mixture;
        // iterate over all sites

        for (size_t site = 0; site < pattern_block_size; ++site)
        {
            // temporary variable storing the likelihood
            double tmp = 0.0;
            // get the pointers to the likelihoods for this site and mixture category
            double* p_site_j   = p_site_mixture;
            // iterate over all starting states
            for (size_t i=0; i<num_states; ++i)
            {
                // add the probability of starting from this state
                tmp += *p_site_j;

                // increment pointers
                ++p_site_j;
            }
            // add the likelihood for this mixture category
            per_site_rate_Likelihoods[site][site_rate_index] += tmp * site_mixture_probs[mixture];

            // increment the pointers to the next site
            p_site_mixture+=this->site_offset;

        } // end-for over all sites (=patterns)

        // increment the pointers to the next mixture category
        p_mixture+=this->mixture_offset;

    } // end-for over all mixtures (=rate categories)

    double prob_invariant = getPInv();
    double oneMinusPInv = 1.0 - prob_invariant;
    std::vector< size_t >::const_iterator patterns = this->pattern_counts.begin();
    if ( prob_invariant > 0.0 )
    {
        // get the mean root frequency vector
        std::vector<double> f;
        if (this->branch_heterogeneous_substitution_matrices == true)
        {
            f = this->getRootFrequencies(0);
        }
        else
        {
            std::vector<std::vector<double> > ff;
            getRootFrequencies(ff);

            std::vector<double> matrix_probs(num_matrices, 1.0/num_matrices);

            if (site_matrix_probs != NULL)
            {
                matrix_probs = site_matrix_probs->getValue();
            }

            f = std::vector<double>(ff[0].size(), 0.0);

            for (size_t matrix = 0; matrix < ff.size(); matrix++)
            {
                // get the root frequencies
                const std::vector<double> &fm = ff[matrix];

                for (size_t i = 0; i < fm.size(); i++)
                {
                    f[i] += fm[i] * matrix_probs[matrix];
                }
            }
        }

        size_t num_site_rates_withInv = num_site_rates + 1;

        for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
        {
            // the first rate category is the invariant
            if ( this->site_invariant[site] == true )
            {
                double ftotal = 0.0;
                for ( size_t c = 0; c < this->invariant_site_index[site].size(); c++ )
                {
                    ftotal += f[this->invariant_site_index[site][c]];
                }

                rv[site][0] = log( prob_invariant * ftotal ) * *patterns;
            }
            else
            {
                rv[site][0] = RbConstants::Double::neginf;
            }

            // the remaining variant rate categories
            for (size_t site_rate_index = 1; site_rate_index < num_site_rates_withInv; ++site_rate_index)
            {
                rv[site][site_rate_index] = log( oneMinusPInv * per_site_rate_Likelihoods[site][site_rate_index - 1] ) * *patterns;

                if ( use_scaling == true )
                {
                    if ( scale_threshold == false )
                    {
                        rv[site][site_rate_index] -= this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site] * *patterns;
                    }
                    else
                    {
                        rv[site][site_rate_index] -= (RbConstants::LN2 * 256 * this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site]) * *patterns;
                    }
                }

            }

        }

    }
    else
    {

        for (size_t site = 0; site < pattern_block_size; ++site, ++patterns)
        {
            for (size_t site_rate_index = 0; site_rate_index < num_site_rates; ++site_rate_index)
            {
                rv[site][site_rate_index] = log( per_site_rate_Likelihoods[site][site_rate_index] ) * *patterns;

                if ( use_scaling == true )
                {
                    if ( scale_threshold == false )
                    {
                        rv[site][site_rate_index] -= this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site] * *patterns;
                    }
                    else
                    {
                        rv[site][site_rate_index] -= (RbConstants::LN2 * 256 * this->per_node_site_mixture_log_scaling_factors[this->active_branch_likelihood[node_index]][node_index][0][site]) * *patterns;
                    }

                }

            }

        }

    }

}



template<class charType>
double RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::sumRootLikelihood( void )
{
    // TODO: Decide here where to get the likelihood: either from BEAGLE directly or as site likelihoods

    std::vector<double> site_likelihoods = std::vector<double>(pattern_block_size,0.0);
    computeRootLikelihoods( site_likelihoods );

    double sum_partial_probs = 0.0;

    for (size_t site = 0; site < pattern_block_size; ++site)
    {
        sum_partial_probs += site_likelihoods[site];
    }

#ifdef RB_MPI

    // we only need to send message if there is more than one process
    if ( num_processes > 1 )
    {

        // send the likelihood from the helpers to the master
        if ( process_active == false )
        {
            // send from the workers the log-likelihood to the master
            MPI_Send(&sum_partial_probs, 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD);
        }

        // receive the likelihoods from the helpers
        if ( process_active == true )
        {
            for (size_t i=active_PID+1; i<active_PID+num_processes; ++i)
            {
                double tmp = 0;
                MPI_Status status;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, int(i), 0, MPI_COMM_WORLD, &status);
                sum_partial_probs += tmp;
            }
        }

        // now send back the combined likelihood to the helpers
        if ( process_active == true )
        {
            for (size_t i=active_PID+1; i<active_PID+num_processes; ++i)
            {
                MPI_Send(&sum_partial_probs, 1, MPI_DOUBLE, int(i), 0, MPI_COMM_WORLD);
            }
        }
        else
        {
            MPI_Status status;
            MPI_Recv(&sum_partial_probs, 1, MPI_DOUBLE, active_PID, 0, MPI_COMM_WORLD, &status);
        }

    }

#endif

    return sum_partial_probs;
}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::swap_taxon_name_2_tip_index(std::string tip1, std::string tip2)
{
    size_t index1 = taxon_name_2_tip_index_map[tip1];
    size_t index2 = taxon_name_2_tip_index_map[tip2];
    taxon_name_2_tip_index_map[tip1] = index2;
    taxon_name_2_tip_index_map[tip2] = index1;
}

/** Swap a parameter of the distribution */
template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if (oldP == homogeneous_clock_rate)
    {
        homogeneous_clock_rate = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == heterogeneous_clock_rates)
    {
        heterogeneous_clock_rates = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else if (oldP == mixture_model)
    {
        mixture_model = static_cast<const TypedDagNode< SiteMixtureModel >* >( newP );
    }
    else if (oldP == homogeneous_rate_matrix)
    {
        homogeneous_rate_matrix = static_cast<const TypedDagNode< RateGenerator >* >( newP );
    }
    else if (oldP == heterogeneous_rate_matrices)
    {
        heterogeneous_rate_matrices = static_cast<const TypedDagNode< RbVector< RateGenerator > >* >( newP );
    }
    else if (oldP == root_frequencies)
    {
        root_frequencies = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == site_matrix_probs)
    {
        site_matrix_probs = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == site_rates)
    {
        site_rates = static_cast<const TypedDagNode< RbVector< double > >* >( newP );
    }
    else if (oldP == site_rates_probs)
    {
        site_rates_probs = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == p_inv)
    {
        p_inv = static_cast<const TypedDagNode< double >* >( newP );
    }
    else if (oldP == tau)
    {
        tau->getValue().getTreeChangeEventHandler().removeListener( this );

        tau = static_cast<const TypedDagNode<Tree>* >( newP );

        tau->getValue().getTreeChangeEventHandler().addListener( this );

        num_nodes = tau->getValue().getNumberOfNodes();
    }

}

template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::touchSpecialization( const DagNode* affecter, bool touch_all )
{
    if ( touched == false )
    {
        touched = true;
        this->storedLnProb = this->lnProb;

        prev_dirty_nodes = dirty_nodes;
        prev_dirty_branches = dirty_branches;
        prev_pmat_dirty_branches = pmat_dirty_branches;
    }


    // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    if ( affecter == heterogeneous_clock_rates )
    {
        const std::set<size_t> &indices = heterogeneous_clock_rates->getTouchedElementIndices();

        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 || indices.size() == this->tau->getValue().getNodes().size() )
        {
            // just flag everyting for recomputation
            touch_all = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
                this->recursivelyFlagBranchDirty( *nodes[*it] );
                this->flagNodeDirtyPmatrix( *it );
            }
        }
    }
    else if ( affecter == heterogeneous_rate_matrices && branch_heterogeneous_substitution_matrices == true)
    {
        const std::set<size_t> &indices = heterogeneous_rate_matrices->getTouchedElementIndices();

        // maybe all of them have been touched or the flags haven't been set properly
        if ( indices.size() == 0 )
        {
            // just flag everyting for recomputation
            touch_all = true;
        }
        else
        {
            const std::vector<TopologyNode *> &nodes = this->tau->getValue().getNodes();
            // flag recomputation only for the nodes
            for (std::set<size_t>::iterator it = indices.begin(); it != indices.end(); ++it)
            {
                this->recursivelyFlagNodeDirty( *nodes[*it] );
                this->recursivelyFlagBranchDirty( *nodes[*it] );
                this->flagNodeDirtyPmatrix( *it );
            }
        }
    }
    else if ( affecter == root_frequencies )
    {
        const TopologyNode &root = this->tau->getValue().getRoot();
        this->recursivelyFlagNodeDirty( root );
    }
    else if ( affecter == p_inv )
    {
        touch_all = true;
    }
    else if ( affecter == site_rates_probs || affecter == site_matrix_probs )
    {
	// This doesn't affect the cached conditional likelihoods (so don't touch all of them).
	// But it does affect the final likelihood (so we need to recompute that).
    }
    else if ( affecter != tau ) // if the topology wasn't the culprit for the touch, then we just flag everything as dirty
    {
        touch_all = true;
    }

    if ( touch_all == true )
    {

        for (std::vector<bool>::iterator it = dirty_branches.begin(); it != dirty_branches.end(); ++it)
        {
            (*it) = true;
        }
        for (std::vector<bool>::iterator it = dirty_nodes.begin(); it != dirty_nodes.end(); ++it)
        {
            (*it) = true;
        }

        // flip the active likelihood pointers
        for (size_t index = 0; index < changed_nodes.size(); ++index)
        {
            if ( changed_nodes[index] == false )
            {
                active_node_likelihood[index] = (active_node_likelihood[index] == 0 ? 1 : 0);
                changed_nodes[index] = true;
            }
        }
        for (size_t index = 0; index < changed_branches.size(); ++index)
        {
            if ( changed_branches[index] == false )
            {
                active_branch_likelihood[index] = (active_branch_likelihood[index] == 0 ? 1 : 0);
                changed_branches[index] = true;
            }
        }
        
        for (std::vector<bool>::iterator it = pmat_dirty_branches.begin(); it != pmat_dirty_branches.end(); ++it)
        {
            (*it) = true;
        }
        
        // flip the active transition probability matrices pointers
        for (size_t index = 0; index < pmat_changed_nodes.size(); ++index)
        {
            if ( pmat_changed_nodes[index] == false )
            {
                active_pmatrices[index] = (active_pmatrices[index] == 0 ? 1 : 0);
                pmat_changed_nodes[index] = true;
            }
        }
    }

}



template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateMarginalNodeLikelihoods( void )
{

    // calculate the root marginal likelihood, then start the recursive call down the tree
    this->computeMarginalRootLikelihood();

    // update the marginal likelihoods by a recursive downpass
    this->recursiveMarginalLikelihoodComputation( tau->getValue().getRoot().getIndex() );

}



/*
 * Update the transition probability matrices for the branch attached to the given node index.
 */
template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateTransitionProbabilities(size_t node_idx)
{
    const TopologyNode* node = tau->getValue().getNodes()[node_idx];

    if ( node->isRoot() == true )
    {
        throw RbException("dnPhyloCTMC called updateTransitionProbabilities for the root node\n");
    }

    // 1. Get the clock rate for the branch
    double rate = 1.0;
    if ( this->branch_heterogeneous_clock_rates == true )
    {
        rate = this->heterogeneous_clock_rates->getValue()[node_idx];
    }
    else if (homogeneous_clock_rate != NULL)
    {
        rate = this->homogeneous_clock_rate->getValue();
    }

    // 2. Handle the mixture model object case.
    if (mixture_model)
    {
        transition_prob_matrices = mixture_model->getValue().calculateTransitionProbabilities(tau->getValue(), node_idx, rate);
        return;
    }

    auto [start_age, end_age] = getStartEndAge(*node);

    // we rescale the rate by the inverse of the proportion of invariant sites
    rate /= ( 1.0 - getPInv() );

    // first, get the rate matrix for this branch
    RateMatrix_JC jc(this->num_states);

    if (this->branch_heterogeneous_substitution_matrices == false )
    {
        // loop now over all per-site rate matrices (could also be only a single one, as by default)
        for (size_t matrix = 0; matrix < this->num_matrices; ++matrix)
        {
            const RateGenerator *rm = nullptr;

            // get the i-th rate matrix
            if ( this->heterogeneous_rate_matrices != NULL )
            {
                rm = &this->heterogeneous_rate_matrices->getValue()[matrix];
            }
            else if ( this->homogeneous_rate_matrix != NULL )
            {
                rm = &this->homogeneous_rate_matrix->getValue();
            }
            else
            {
                rm = &jc;
            }

	    // The rm can change behind our back if the user redefines it.
	    if (rm->size() != num_states)
		throw RbException()<<"Rate generator with "<<rm->size()<<" states does not match data with "<<num_states<<" states";

            // now also get the site specific rates
            for (size_t j = 0; j < this->num_site_rates; ++j)
            {
                double r = 1.0;
                if ( this->rate_variation_across_sites == true )
                {
                    r = this->site_rates->getValue()[j];
                }

                rm->calculateTransitionProbabilities( start_age, end_age,  rate * r, this->transition_prob_matrices[j*this->num_matrices + matrix] );
            }
        }
    }
    else
    {
        const RateGenerator *rm = nullptr;

        if ( this->heterogeneous_rate_matrices != NULL )
        {
            rm = &this->heterogeneous_rate_matrices->getValue()[node_idx];
        }
        else if ( this->homogeneous_rate_matrix != NULL )
        {
            rm = &this->homogeneous_rate_matrix->getValue();
        }
        else
        {
            rm = &jc; // BAD!!
        }

	// The rm can change behind our back if the user redefines it.
	if (rm->size() != num_states)
	    throw RbException()<<"Rate generator with "<<rm->size()<<" states does not match data with "<<num_states<<" states";

        for (size_t j = 0; j < this->num_site_rates; ++j)
        {
            double r = 1.0;
            if ( this->rate_variation_across_sites == true )
            {
                r = this->site_rates->getValue()[j];
            }

            rm->calculateTransitionProbabilities( start_age, end_age,  rate * r, this->transition_prob_matrices[j] );
        }
    }
}


/*
 * Update the transition probability matrices for the branch attached to the given node index.
 */
template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateTransitionProbabilityMatrix(size_t node_idx)
{
    const TopologyNode* node = tau->getValue().getNodes()[node_idx];
    
    if (node->isRoot()) throw RbException("dnPhyloCTMC called updateTransitionProbabilityMatrix for the root node\n");
    
    size_t pmat_offset = this->active_pmatrices[node_idx] * this->active_P_matrix_offset + node_idx * this->pmat_node_offset;

    // 1. Get the clock rate for the branch
    double rate = 1.0;
    if ( this->branch_heterogeneous_clock_rates == true )
    {
        rate = this->heterogeneous_clock_rates->getValue()[node_idx];
    }
    else if (homogeneous_clock_rate != NULL)
    {
        rate = this->homogeneous_clock_rate->getValue();
    }

    // 2. Handle the mixture model object case.
    if (mixture_model)
    {
        auto& model = mixture_model->getValue();
        for(int m=0; m < model.getNumberOfComponents(); m++)
            pmatrices[pmat_offset + m] = model.getComponent(m).calculateTransitionProbabilities(tau->getValue(), node_idx, rate);
        return;
    }

    // 3. Handle non-mixture-models.
    auto [start_age, end_age] = getStartEndAge(*node);

    // first, get the rate matrix for this branch
    RateMatrix_JC jc(this->num_states); // BAD!!!
    const RateGenerator *rm = &jc;

    // we rescale the rate by the inverse of the proportion of invariant sites
    rate /= ( 1.0 - getPInv() );

    if (this->branch_heterogeneous_substitution_matrices == false )
    {
        for (size_t matrix = 0; matrix < this->num_matrices; ++matrix)
        {
            if ( this->heterogeneous_rate_matrices != NULL )
            {
                rm = &this->heterogeneous_rate_matrices->getValue()[matrix];
            }
            else if ( this->homogeneous_rate_matrix != NULL )
            {
                rm = &this->homogeneous_rate_matrix->getValue();
            }

	    // The rm can change behind our back if the user redefines it.
	    if (rm->size() != num_states)
		throw RbException()<<"Rate generator with "<<rm->size()<<" states does not match data with "<<num_states<<" states";

            for (size_t j = 0; j < this->num_site_rates; ++j)
            {
                double r = 1.0;
                if ( this->rate_variation_across_sites == true )
                {
                    r = this->site_rates->getValue()[j];
                }
                
                rm->calculateTransitionProbabilities( start_age, end_age,  rate * r, this->pmatrices[pmat_offset + j * this->num_matrices + matrix] );
            }
        }
    }
    else
    {
        if ( this->heterogeneous_rate_matrices != NULL )
        {
            rm = &this->heterogeneous_rate_matrices->getValue()[node_idx];
        }
        else if ( this->homogeneous_rate_matrix != NULL )
        {
            rm = &this->homogeneous_rate_matrix->getValue();
        }

	// The rm can change behind our back if the user redefines it.
	if (rm->size() != num_states)
	    throw RbException()<<"Rate generator with "<<rm->size()<<" states does not match data with "<<num_states<<" states";

	for (size_t j = 0; j < this->num_site_rates; ++j)
        {
            double r = 1.0;
            if ( this->rate_variation_across_sites == true )
            {
                r = this->site_rates->getValue()[j];
            }
            
            rm->calculateTransitionProbabilities( start_age, end_age,  rate * r, this->pmatrices[pmat_offset + j] );
        }
    }
}


/*
 * Update the transition probability matrices for each branch that is marked dirty.
 */
template<class charType>
void RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateTransitionProbabilityMatrices( void )
{
    
    std::vector<TopologyNode*> nodes = tau->getValue().getNodes();
        
    for (std::vector<TopologyNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        size_t node_index = (*it)->getIndex();
        if (pmat_dirty_branches[node_index] == true)
        {
            if ((*it)->isRoot() == false)
            {
                updateTransitionProbabilityMatrix(node_index);
            }

            // mark as computed
            pmat_dirty_branches[node_index] = false;
        }
    }
    
}




/* Collection of BEAGLE helper methods (when beagle-lib is available). */
#if defined( RB_BEAGLE )

template<class charType>
void
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::initializeBeagleInstances( void )
{
    // Return and do nothing if we are not in BEAGLE
    charType tmp_char;
    if ( RbSettings::userSettings().getUseBeagle() == false || ( tmp_char.getDataType() != "Protein" && tmp_char.getDataType() != "DNA" ) )
    {
        return;
    }

    // For now we do not allow for partitioned analyses. Maybe in the future
    //size_t num_models = 1;

#if defined ( RB_BEAGLE_INFO )
    RBOUT("BEAGLE Instances: " + std::to_string(num_models) + "\n");
#endif 

    this->num_mixtures           = (this->homogeneous_rate_matrix)
                                 ? 1 : this->heterogeneous_rate_matrices->getValue().size();
    
    bool   b_use_scaling         = RbSettings::userSettings().getUseScaling();

    int    b_tipCount            = int( this->num_tips );
    int    b_partialsBufferCount = 2 * this->num_nodes
                                 + ( this->using_ambiguous_characters
                                   ? this->num_tips
                                   : 0
                                   );
    int    b_compactBufferCount  = this->num_tips
                                 - ( this->using_ambiguous_characters
                                   ? this->num_tips
                                   : 0 );
    int    b_stateCount          = this->num_states;
    int    b_patternCount        = this->pattern_block_size;
    int    b_eigenBufferCount    = this->num_mixtures * 2;
    int    b_matrixBufferCount   = this->num_mixtures * this->num_nodes * 2;
    int    b_categoryCount       = this->num_site_rates +
                                   ( this->getPInv() > std::numeric_limits<double>::epsilon()
                                     ? 1 : 0 );
    // @todo: need to fix this, don't need to allocate as many buffers (Sebastian)
    int    b_scaleBufferCount    = b_use_scaling ? (2 + 2 * this->num_nodes + ( this->using_ambiguous_characters ? this->tau->getValue().getNumberOfTips() : 0)) : 0;
    BeagleInstance *b_instance   = new BeagleInstance();

    // Create the BEAGLE instance
    b_instance->createBEAGLE( b_tipCount,
                              b_partialsBufferCount,
                              b_compactBufferCount,
                              b_stateCount,
                              b_patternCount,
                              b_eigenBufferCount,
                              b_matrixBufferCount,
                              b_categoryCount,
                              b_scaleBufferCount );

    // free old beagle instance
    if ( beagle_instance != NULL )
    {
        freeBeagleInstances();
    }
    
    // And set the global instance
    beagle_instance = b_instance;
        

    // Initialize tips for models
    this->initializeBeagleTips();
    
    if ( RbSettings::userSettings().getUseScaling() == true )
    {
        for (size_t i=0; i<num_nodes; ++i)
        {
            beagleResetScaleFactors(beagle_instance->getResourceID(), i);
            beagleResetScaleFactors(beagle_instance->getResourceID(), i+num_nodes);
        }
        beagleResetScaleFactors(beagle_instance->getResourceID(), 2*num_nodes);
        beagleResetScaleFactors(beagle_instance->getResourceID(), 2*num_nodes+1);
    }
}


template<class charType>
void
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::freeBeagleInstances( void )
{
    charType tmp_char;
    if ( RbSettings::userSettings().getUseBeagle() == true && ( tmp_char.getDataType() == "Protein" || tmp_char.getDataType() == "DNA" ) )
    {
        if ( beagle_instance != NULL )
        {
            beagle_instance->freeBEAGLE();
            delete beagle_instance;
            beagle_instance = NULL;
        }
        
    }
}


template<class charType>
void
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::initializeBeagleTips( void )
{
    // TODO - This method does not work with patitioned analyses!!
    
    int     b_tipIndex   = 0;
    int*    b_inStates   = NULL;
    double* b_inPartials = NULL;
    int     b_ret_code   = 0;

    if ( this->using_ambiguous_characters == true ) {
        b_inPartials = new double[this->pattern_block_size * this->num_states];
    } else {
        b_inStates = new int[this->pattern_block_size];
    }

    //-- Iterate over the all leaf nodes in the tree.
    std::vector<TopologyNode*> nodes = this->tau->getValue().getNodes();
    for ( std::vector<TopologyNode*>::iterator it = nodes.begin() ; it != nodes.end() ; ++it ) {
        if ( (*it)->isTip() ) {
            b_tipIndex                                  = (*it)->getIndex();
            size_t data_tip_index                       = this->taxon_name_2_tip_index_map[(*it)->getName()];
            const std::vector<bool> &gap_node           = this->gap_matrix[data_tip_index];
            const std::vector<std::uint64_t> &char_node = this->char_matrix[data_tip_index];
            const std::vector<RbBitSet> &amb_char_node  = this->ambiguous_char_matrix[data_tip_index];

            // iterate over all sites
            for ( size_t b_pattern = 0; b_pattern < this->pattern_block_size; ++b_pattern ) {
                // is this site a gap?
                if ( gap_node[b_pattern] ) {
                    if ( this->using_ambiguous_characters == true ) {
                        for ( size_t c = 0; c < this->num_states; ++c ) {
                            b_inPartials[b_pattern * this->num_states + c] = 1.0;
                        }
                    }
                    else {
                        b_inStates[b_pattern] = (int) this->num_states;
                    }
                } else {
                    if ( this->using_ambiguous_characters == true ) {
                        for ( size_t c = 0; c < this->num_states; ++c ) {
                            b_inPartials[b_pattern * this->num_states + c] =
                                (amb_char_node[b_pattern][c] ? 1.0 : 0.0);
                        }
                    } else {
                        b_inStates[b_pattern] = (int) char_node[b_pattern];
                    }
                }
            }

            if ( this->using_ambiguous_characters == true ) {
                b_ret_code = beagleSetTipPartials( beagle_instance->getResourceID(),
                                                   b_tipIndex,
                                                   b_inPartials );
                //-- Check to see if we could set the tips
                if (b_ret_code != 0) {
                    throw RbException("Could not set tip partials for model" +
                                      BeagleUtilities::printErrorCode(b_ret_code));
                }
            } else {
                b_ret_code = beagleSetTipStates( this->beagle_instance->getResourceID(),
                                                 b_tipIndex,
                                                 b_inStates );
                //-- Check to see if we could set the tips
                if (b_ret_code != 0) {
                    throw RbException("Could not set tip partials for model" +
                                      BeagleUtilities::printErrorCode(b_ret_code));
                }
            }
            
        }
    }

    this->b_inPatternWeights = std::vector<double>(this->pattern_block_size);
    for ( size_t b_pattern = 0; b_pattern < this->pattern_block_size; ++b_pattern ) {
        this->b_inPatternWeights[b_pattern] = (double) this->pattern_counts[b_pattern];
    }

    //-- Set pattern weights 
    b_ret_code = beagleSetPatternWeights( this->beagle_instance->getResourceID(),
                                          &this->b_inPatternWeights[0] );
    //-- Check to see if we could set the pattern weights
    if (b_ret_code != 0) {
        throw RbException("Could not set pattern weights for model " +
                          BeagleUtilities::printErrorCode(b_ret_code));
    }

    delete[] b_inStates;
    delete[] b_inPartials;
}


template<class charType>
double
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::calculateBranchLength
( const TopologyNode &node
, size_t node_index
)
{
    double branch_len;
    double rate;

    if ( this->branch_heterogeneous_clock_rates == true ) {
        rate = this->heterogeneous_clock_rates->getValue()[node_index];
    } else if ( this->homogeneous_clock_rate != NULL) {
        rate = this->homogeneous_clock_rate->getValue();
    } else {
        rate = 1.0;
    }

    rate /= ( 1.0 - getPInv() );

    branch_len = rate * node.getBranchLength();
    if ( branch_len < 0 ) {
      throw RbException("Error : Negative branch length!");
    }

    return branch_len;
}


template<class charType>
void
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateBeagleEigensystems
( void )
{
    /*
     * This method updates the eigensystems for a BEAGLE instance.
     * Note: This method is not set up to handle partitioned datasets.
     * Work still needs to be done for this.
     */

    //-- Substitution model indexing.
    size_t                               b_model_idx;
    size_t                               b_num_models;

    //-- Discrete rate matrices for models.
    const AbstractRateMatrix*            arm_ptr;
    RbVector<MatrixReal>                 rate_matrices;

    //-- Stationary frequencies for models.
    std::vector<std::vector<double>>     model_pi_vectors;
    size_t                               b_stateFrequenciesIndex;
    std::vector<double>                  b_inStateFrequencies;

    //-- Mixture weights
    std::vector<double>                  b_patternWeights;

    //-- BEAGLE return code for error checking.
    int                                  b_ret_code;

    //-- Eigensystem temporary structures.
#if defined ( RB_USE_EIGEN3 )
    //-- New  Eigen library way... but requires 2 copies...
    Eigen::EigenSolver<Eigen::MatrixXd>  eigen_system;
    Eigen::MatrixXd                      rate_matrix(this->num_states, this->num_states);
    Eigen::VectorXd                      eigenvalues;
    Eigen::MatrixXd                      eigenvectors;
    Eigen::MatrixXd                      inv_eigenvectors;
#endif /* RB_BEAGLE_DEBUG */

    //-- Get the 'root frequencies' (stationary distributions) for all models.
    this->getRootFrequencies(model_pi_vectors);
    
    //-- Get the rate matrices of all models.
    if ( this->homogeneous_rate_matrix ) {
        arm_ptr = dynamic_cast<const AbstractRateMatrix*>(&this->homogeneous_rate_matrix->getValue());
        rate_matrices.push_back(arm_ptr->getRateMatrix());
    } else {
        for ( size_t i = 0; i < this->heterogeneous_rate_matrices->getValue().size(); ++i ) {
            arm_ptr = dynamic_cast<const AbstractRateMatrix*>(&this->heterogeneous_rate_matrices->getValue()[i]);
            rate_matrices.push_back(arm_ptr->getRateMatrix());
        }
    }

    for ( size_t i = 0; i < rate_matrices.size(); ++i ) {
        //-- TODO : Maybe add checks to only update if eigensystem changes (use touched bitmap)
        b_model_idx             = i;
        b_stateFrequenciesIndex = i;
        b_inStateFrequencies    = model_pi_vectors[i];

        b_ret_code =
            beagleSetStateFrequencies( this->beagle_instance->getResourceID(),
                                       b_stateFrequenciesIndex,
                                       &b_inStateFrequencies[0] );
        if ( b_ret_code != 0 ) {
            throw RbException( "Could not set state frequencies for model '" +
                               std::to_string(i) + "'. " +
                               BeagleUtilities::printErrorCode(b_ret_code));
        }
        
        //-- Set eigensystem for model.
#if defined ( RB_USE_EIGEN3 )
        //-- Copy into Eigen structure
        for (size_t j = 0; j < this->num_states; ++j ) {
            for (size_t k = 0; k < this->num_states; ++k ) {
                rate_matrix(j,k) = (rate_matrices[i])[j][k]; 
            }
        }
        eigen_system.compute(rate_matrix, true); //-- compute eigen values and vectors

        // TODO - There should be a way to map contiguous memory to a Eigen datatype.
        // Need to figure out how to extract from RbVector... if possible.
        //Eigen::Map<Eigen::MatrixXd> rate_matrix(rate_matrices[i].elements.data(), this->num_states, this->num_states);
        //eigen_system.compute(rate_matrix, true);
        
        eigenvalues      = eigen_system.eigenvalues().real();
        eigenvectors     = eigen_system.eigenvectors().real().transpose();
        inv_eigenvectors = eigen_system.eigenvectors().inverse().real().transpose();
        
        //-- Should at some point figure out if there is a way just to use eigen matrix data pointer instead of copy.
        std::vector<double> b_flat_eigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.rows() * eigenvalues.cols());
        std::vector<double> b_flat_eigenvectors(eigenvectors.data(), eigenvectors.data() + eigenvectors.rows() * eigenvectors.cols());
        std::vector<double> b_flat_inv_eigenvectors(inv_eigenvectors.data(), inv_eigenvectors.data() + inv_eigenvectors.rows() * inv_eigenvectors.cols());

        // Set eigensystem for model in BEAGLE
        b_ret_code =
            beagleSetEigenDecomposition( this->beagle_instance->getResourceID()
                                       , b_model_idx
                                       , &b_flat_eigenvectors[0]
                                       , &b_flat_inv_eigenvectors[0]
                                       , &b_flat_eigenvalues[0]
                                       );
        if ( b_ret_code != 0 ) {
            throw RbException( "Could not set eigen decomposition for model '"
                               + std::to_string(i) + "'. "
                               + BeagleUtilities::printErrorCode(b_ret_code));
        }
    
#endif /* RB_USE_EIGEN3 */
    } // End models for loop
}


template<class charType>
void
RevBayesCore::AbstractPhyloCTMCSiteHomogeneous<charType>::updateBeagleSiteRates
( void )
{
    //-- Multiple sets of site rates not supported, so we only use index '0'.
    int b_categoryWeightsIndex = 0;

    //-- Return code for beagle operations.
    int b_ret_code = 0;

    //-- Get invariable sites
    double pInv = this->getPInv();

    //-- Set the site rate probabilities (corresponds with beagle 'category weights').
    std::vector<double> rates_probs(this->num_site_rates, 1.0 / this->num_site_rates);
    if ( this->site_rates_probs != NULL ) {
        rates_probs = this->site_rates_probs->getValue();
    }
    
    //-- Set the category weights and rates.
    if ( this->num_site_rates > 1 ) {
        //-- Set the site rate probabilities (corresponds with beagle 'category weights').
        std::vector<double> rates_probs(this->num_site_rates, 1.0 / this->num_site_rates);
        this->b_inCategoryWeights = rates_probs; 

        //-- Set the site rate values (corresponds with beagle 'category rates').
        this->b_inCategoryRates = dynamic_cast<const RbVector<double>&>(this->site_rates->getValue());
    } else {
        this->b_inCategoryWeights = { 1.0 };
        this->b_inCategoryRates   = { 1.0 };
    }

    //-- If we are using invariable sites we need to add it to the model
    //   and adjust the category weights and rates accordingly so that we
    //   do not violate the model.
    if ( pInv > 0.0 ) {
        //-- We need the category weights to sum to 1.0.
        for ( size_t i = 0; i < this->b_inCategoryWeights.size(); ++i ) {
            this->b_inCategoryWeights[i] -= (pInv / this->b_inCategoryWeights.size());
        }

        //-- And push the invariable sites weight and rate, respectively.
        this->b_inCategoryWeights.push_back(pInv);
        this->b_inCategoryRates.push_back(0.0);
    }

    //-- Set BEAGLE category rates for rates mixture
    b_ret_code =
        beagleSetCategoryRatesWithIndex( this->beagle_instance->getResourceID(),
                                         b_categoryWeightsIndex,
                                         &this->b_inCategoryRates[0] );
    if ( b_ret_code != 0 ) {
        throw RbException( "Could not set category rates for ASRV model " +
                           std::to_string(b_categoryWeightsIndex) + " : " +
                           BeagleUtilities::printErrorCode(b_ret_code) );
    }
    
    //-- Set BEAGLE category weights for mixture
    b_ret_code =
        beagleSetCategoryWeights( this->beagle_instance->getResourceID(),
                                  b_categoryWeightsIndex,
                                  &this->b_inCategoryWeights[0] );
    if (b_ret_code != 0) {
        throw RbException( "Could not set category weights for ASRV model " +
                           std::to_string(b_categoryWeightsIndex) + " : " +
                           BeagleUtilities::printErrorCode(b_ret_code));
    }
}

#endif /* END RB_BEAGLE */

#endif /* END AbstractPhyloCTMCSiteHomogeneous_H */
