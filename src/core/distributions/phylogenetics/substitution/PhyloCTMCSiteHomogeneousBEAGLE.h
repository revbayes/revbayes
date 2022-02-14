/**
 * @file
 * 
 * This file contains BEAGLE specific methods that implement the
 * virtual functions from AbstractPhyloCTMCSiteHomogeneous.
 *
 * There are still some things that we need to finish:
 *
 *     - Using transition matrices with ASRV does not return correct likelihoods. This is probably
 *       because there is a transition matrix per site rate. Need to allocate more memory in BEAGLE
 *       and adjust accordingly.
 *
 *     - The allocation and deallocation of 'partialLikelihoods' should be removed when using BEAGLE
 *       in 'AbstractPhyloCTMCSiteHomogeneous.h'. However, it appears that something is still required
 *       inside of 'resizeLikelihoodvectors' or 'scale'. Determine what needs to be refactored.
 * 
 *     - BEAGLE instances should be in a vector. Then we can use with MPI and mixture models.
 *
 *     - Configure for usage with MPI.
 *
 *     - Configure for usage with mixture models.
 *
 */

#ifndef PhyloCTMCSiteHomogeneousBEAGLE_H
#define PhyloCTMCSiteHomogeneousBEAGLE_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"


//#define RB_BEAGLE_DEBUG
//#define RB_BEAGLE_DEBUG_TIP


namespace RevBayesCore
{

    template<class charType>
    class PhyloCTMCSiteHomogeneousBEAGLE : public AbstractPhyloCTMCSiteHomogeneous<charType>
    {

        public:

            //----====  Constructors  ====----

            //-- Default constructor
            PhyloCTMCSiteHomogeneousBEAGLE ( const TypedDagNode<Tree>* t
                                           , size_t nChars
                                           , bool c
                                           , size_t nSites
                                           , bool amb
                                           , bool internal
                                           , bool gapmatch
                                           );

            //-- Clone constructor
            PhyloCTMCSiteHomogeneousBEAGLE* clone ( void ) const;


            //-- Destructor
            virtual ~PhyloCTMCSiteHomogeneousBEAGLE ( void );


        protected:

            //----====  Protected Methods  ====----

            //-- Return the computed likelihood.
            virtual double sumRootLikelihood             ( void );

            //-- BEAGLE compute lnLikelihood of a rooted tree.
            virtual void   computeRootLikelihood         ( size_t root
                                                         , size_t l
                                                         , size_t r
                                                         );

            //-- BEAGLE compute lnLikelihood of an unrooted tree.
            virtual void   computeRootLikelihood         ( size_t root
                                                         , size_t l
                                                         , size_t r
                                                         , size_t m
                                                         );

            //-- Collect a BEAGLE operation for an internal node into the computation queue.
            virtual void   computeInternalNodeLikelihood ( const TopologyNode &n
                                                         , size_t nIdx
                                                         , size_t l
                                                         , size_t r
                                                         );

            //-- Collect a BEAGLE operation for an internal node into the computation queue.
            virtual void   computeInternalNodeLikelihood ( const TopologyNode &n
                                                         , size_t nIdx
                                                         , size_t l
                                                         , size_t r
                                                         , size_t m
                                                         );

            //-- Collect a BEAGLE operation for a leaf node into the computation queue.
            virtual void   computeTipLikelihood          ( const TopologyNode &node
                                                         , size_t nIdx
                                                         );

        private:

            //----====  Private Variables  ====----

            //-- Accumulate BEAGLE lnLikelihood across all models.
            double ln_beagle_probability;
        
    };

} //-- End namespace


#include "ConstantNode.h"
#include "DiscreteCharacterState.h"
#include "RandomNumberFactory.h"

#include <cmath>
#include <cstring>



template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::PhyloCTMCSiteHomogeneousBEAGLE
  ( const TypedDagNode<Tree>* t
  , size_t nChars
  , bool c
  , size_t nSites
  , bool amb
  , bool internal
  , bool gapmatch
  ) : AbstractPhyloCTMCSiteHomogeneous<charType> ( t
                                                 , nChars
                                                 , 1
                                                 , c
                                                 , nSites
                                                 , amb
                                                 , internal
                                                 , gapmatch
                                                 )
{ }



template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::~PhyloCTMCSiteHomogeneousBEAGLE ( void )
{
    // We don't delete the parameters, because they might be used somewhere else too.
    // The model needs to do that!

    // When we are done clean up BEAGLE instance
    this->freeBeagleInstance();
}



template<class charType>
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>*
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::clone ( void ) const
{
    return new PhyloCTMCSiteHomogeneousBEAGLE<charType>(*this);
}


template<class charType>
double
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::sumRootLikelihood (void )
{
    //-- We have already computed the ln_beagle_probability, so just return it
    return this->ln_beagle_probability;
}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeRootLikelihood
  ( size_t root
  , size_t left
  , size_t right
  )
{
#if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << std::endl << "Using rooted lnLikelihood calculation" << std::endl;
    RBOUT(ss.str());
#endif /* RB_BEAGLE_DEBUG */

    size_t b_model_idx;
    size_t num_taxa  = (this->num_nodes + 1) / 2;
    
    size_t root_idx  = root;
    if ( RbSettings::userSettings().getUseBeagleLikelihoodStoring() == true )
    {
        root_idx  = root + this->num_nodes * this->activeLikelihood[root];
    }
    
    size_t left_idx  = left;
    size_t right_idx = right;
    if ( RbSettings::userSettings().getUseBeagleLikelihoodStoring() == true )
    {
        left_idx  = left   + this->num_nodes * this->activeLikelihood[left];
        right_idx = right  + this->num_nodes * this->activeLikelihood[right];
    }

    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials  = (left  < num_taxa) ? left  : left_idx;
    size_t right_partials = (right < num_taxa) ? right : right_idx;

    //-- Push the last operation onto the queue
    BeagleOperation b_operation =
        { .destinationPartials    = (int) root_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };
    //this->b_ops.push_back(b_operation);  //-- TESTING! -- not in original

    //-- BEAGLE model parameters.
    int     b_parentBufferIndices     = (int) root_idx;
    int     b_childBufferIndices      = (int) NULL;
    int     b_probabilityIndices      = (int) NULL;
    int*    b_firstDerivativeIndices  = NULL;
    int*    b_secondDerivativeIndices = NULL;
    int     b_categoryWeightsIndices  = 0; //(int) model;  //0;
    int     b_stateFrequenciesIndices = 0; //(int) model;  //0;
    int     b_cumulativeScaleIndices  = BEAGLE_OP_NONE;
    int     b_count                   = 1;
    double  b_outSumLogLikelihood     = 0; //     = NULL; //0;
    double* b_outSumFirstDerivative   = NULL;
    double* b_outSumSecondDerivative  = NULL;

    //-- Return codes for BEAGLE operations.
    int b_code_update_transitions;
    int b_code_update_partials;
    int b_code_calc_edges;

    //-- Update rates across sites 
    this->updateBeagleSiteRates();

#if defined ( RB_BEAGLE_DEBUG )
    ss << "updated site rates" << std::endl;
#endif /* RB_BEAGLE_DEBUG */

#if defined ( RB_BEAGLE_EIGEN )
    this->updateBeagleEigensystem();


    for ( size_t i = 0; i < this->num_site_mixtures; ++i )
    {
        //b_model_idx = this->active_eigen_system[i];
        b_model_idx = i + 0* this->active_eigen_system[i] * this->num_site_mixtures;

        //-- Update all transition matrices for model i.
        b_code_update_transitions =
            beagleUpdateTransitionMatrices( this->beagle_instance->getResourceID()
                                          , b_model_idx
                                          , &this->b_node_indices[0]
                                          , NULL
                                          , NULL
                                          , &this->b_branch_lengths[0]
                                          , this->b_branch_lengths.size()
                                          );
        if ( b_code_update_transitions != 0 )
        {
            throw RbException( "Could not update transition matrix for model '"
                             + std::to_string(i) + "'. "
                             + BeagleUtilities::printErrorCode(b_code_update_transitions));
        }
    }

    //-- Calculate and update all partial likelihood buffers
    b_code_update_partials = beagleUpdatePartials( this->beagle_instance->getResourceID()
                                                 , &this->b_ops[0]
                                                 , this->b_ops.size()
                                                 , BEAGLE_OP_NONE
                                                 );
#else
    b_code_update_partials = beagleUpdatePartials( this->beagle_instance->getResourceID()
                                                 , &b_operation
                                                 , 1
                                                 , BEAGLE_OP_NONE
                                                 );
#endif /* RB_BEAGLE_EIGEN */

#if defined ( RB_BEAGLE_DEBUG )
    ss << "updated partials" << std::endl;
#endif /* RB_BEAGLE_DEBUG */
    
	if ( b_code_update_partials != 0 )
	{
        throw RbException( "Could not update partials for models '"
	  		             + BeagleUtilities::printErrorCode(b_code_update_partials));
	}

    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();

    //-- Calclulate the lnLikelihood of the model
    b_code_calc_edges =
        beagleCalculateEdgeLogLikelihoods( this->beagle_instance->getResourceID()
                                         , &b_parentBufferIndices
                                         , &b_childBufferIndices
                                         , &b_probabilityIndices
                                         , b_firstDerivativeIndices
                                         , b_secondDerivativeIndices

                                         , &b_categoryWeightsIndices
                                         , &b_stateFrequenciesIndices

				                         , &b_cumulativeScaleIndices
                                         , b_count
                                         , &b_outSumLogLikelihood
                                         , b_outSumFirstDerivative
                                         , b_outSumSecondDerivative
                                         );
    if ( b_code_calc_edges != 0 )
    {
	    throw RbException("Could not calculate edge log likelihood for models '"
                          + BeagleUtilities::printErrorCode(b_code_calc_edges));
    }

    this->ln_beagle_probability = b_outSumLogLikelihood;

#if defined ( RB_BEAGLE_DEBUG )
    RBOUT(ss.str());
#endif /* RB_BEAGLE_DEBUG */
}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeRootLikelihood
  ( size_t root
  , size_t left
  , size_t right
  , size_t middle
  )
{
#if defined ( RB_BEAGLE_DEBUG )
    std::stringstream ss;
    ss << std::endl << "Using unrooted lnLikelihood calculation" << std::endl;
#endif /* RB_BEAGLE_DEBUG */

    size_t b_model_idx;

    //-- Get the number of taxa in the tree
    size_t num_taxa = (this->num_nodes + 2) / 2;

    //-- Determine the number of unique eigensystems we will have.
    size_t b_num_models = this->num_site_mixtures - this->num_site_rates;
    if ( b_num_models < 1 ) { b_num_models = 1; }
    
    //-- Calculate the node indices accounting for active/inactive offests.
    size_t root_idx  = root   + this->num_nodes * this->activeLikelihood[root];
    size_t mid_idx   = middle + this->num_nodes * this->activeLikelihood[middle];
    size_t left_idx  = left   + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right  + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once in BEAGLE, so we cant have active/inactive offests.
    size_t mid_partials   = (middle < num_taxa) ? middle : mid_idx;
    size_t left_partials  = (left   < num_taxa) ? left   : left_idx;
    size_t right_partials = (right  < num_taxa) ? right  : right_idx;

#if defined ( RB_USE_EIGEN3 )
    //-- Set the ASRV category for each site. Since we do not allow for partitions, this is always 0.
    std::vector<int> categoryIndicesASRV(this->pattern_block_size, 0);
    //-- And Update the respective ASRV BEAGLE buffers.
    this->updateBeagleSiteRates();
#else
    std::vector<int> categoryIndicesASRV = NULL;
#endif
    
    //-- Configure BEAGLE model parameters.
    int      b_parentBufferIndices     = (int) root_idx;
    int      b_childBufferIndices      = (int) mid_partials;
    int      b_probabilityIndices      = (int) mid_idx;
    int *    b_firstDerivativeIndices  = NULL;
    int *    b_secondDerivativeIndices = NULL;
    int *    b_categoryWeightsIndices  = &categoryIndicesASRV[0];
    int      b_stateFrequenciesIndices = 0; //(int) model;  //0;
    int      b_cumulativeScaleIndices  = BEAGLE_OP_NONE;
    int      b_count                   = 1;
    double   b_outSumLogLikelihood     = -1 * std::numeric_limits<double>::max();
    double * b_outSumFirstDerivative   = NULL;
    double * b_outSumSecondDerivative  = NULL;

    //-- Create BEAGLE operation.
    BeagleOperation b_operation =
        { .destinationPartials    = (int) root_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    //-- Push operation and root index onto respective vectors to prepare for likelihood calculation.
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(root_idx);

    //-- Return codes for BEAGLE operations.
    int b_code_update_transitions;
    int b_code_update_partials;
    int b_code_calc_edges;

    
#if defined ( RB_USE_EIGEN3 )
    //-- Update Eigensystem BEAGLE buffers
    this->updateBeagleEigensystem();

    for ( size_t i = 0; i < b_num_models; ++i )
    {
        b_model_idx = i + this->active_eigen_system[i] * this->num_site_mixtures;

        //-- Update all transition matrices for model i.
        b_code_update_transitions =
            beagleUpdateTransitionMatrices( this->beagle_instance->getResourceID()
                                          , b_model_idx
                                          , &this->b_node_indices[0]
                                          , NULL
                                          , NULL
                                          , &this->b_branch_lengths[0]
                                          , this->b_branch_lengths.size()
                                          );
        if ( b_code_update_transitions != 0 )
        {
            throw RbException( "Could not update transition matrix for model '"
                             + std::to_string(i) + "'. "
                             + BeagleUtilities::printErrorCode(b_code_update_transitions));
        }
    }

    //-- Calculate and update all partial likelihood buffers
    b_code_update_partials = beagleUpdatePartials( this->beagle_instance->getResourceID()
                                                 , &this->b_ops[0]
                                                 , this->b_ops.size()
                                                 , BEAGLE_OP_NONE
                                                 );

#else
    // TODO - remove this from set eigensystem as only homogeneous models supported
    std::vector<std::vector<double>> model_pi_vectors;
    this->getRootFrequencies(model_pi_vectors);
    int b_stateFrequenciesIndex = 0;
    std::vector<double> b_inStateFrequencies = model_pi_vectors[0];

    beagleSetStateFrequencies( this->beagle_instance->getResourceID()
                             , b_stateFrequenciesIndex
                             , &b_inStateFrequencies[0]
                             );

    b_code_update_partials = beagleUpdatePartials( this->beagle_instance->getResourceID()
                                                 , &b_operation
                                                 , 1
                                                 , BEAGLE_OP_NONE
                                                 );
#endif //-- RB_BEAGLE_EIGEN

    if ( b_code_update_partials != 0 )
	{
        throw RbException( "Could not update partials for models '"
	  		             + BeagleUtilities::printErrorCode(b_code_update_partials));
	}

    //-- Reset the beagle operations queues 
    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();

    //-- Calclulate the lnLikelihood of the model
    b_code_calc_edges =
        beagleCalculateEdgeLogLikelihoods( this->beagle_instance->getResourceID()
                                         , &b_parentBufferIndices
                                         , &b_childBufferIndices
                                         , &b_probabilityIndices
                                         , b_firstDerivativeIndices
                                         , b_secondDerivativeIndices

                                         , b_categoryWeightsIndices
                                         , &b_stateFrequenciesIndices

				                         , &b_cumulativeScaleIndices
                                         , b_count
                                         , &b_outSumLogLikelihood
                                         , b_outSumFirstDerivative
                                         , b_outSumSecondDerivative
                                         );
    if ( b_code_calc_edges != 0 )
    {
	    throw RbException("Could not calculate edge log likelihood for models '"
                          + BeagleUtilities::printErrorCode(b_code_calc_edges));
    }

    this->ln_beagle_probability = b_outSumLogLikelihood;
}



template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeInternalNodeLikelihood
  ( const TopologyNode &node
  , size_t node_index
  , size_t left
  , size_t right
  )
{
    //-- Calculate the node indices accounting for active/inactive offests.
    size_t b_node_idx  = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t b_left_idx  = left       + this->num_nodes * this->activeLikelihood[left];
    size_t b_right_idx = right      + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t b_left_partials  = (node.getChild(0).isTip()) ? left  : b_left_idx;
    size_t b_right_partials = (node.getChild(1).isTip()) ? right : b_right_idx;

    // Compute the branch length
    double b_branch_length = this->calculateBranchLength(node, node_index);

    // Construct the BEAGLE operation that will be pushed onto the compute queue.
    BeagleOperation b_operation =
        { .destinationPartials    = (int) b_node_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) b_left_partials
        , .child1TransitionMatrix = (int) b_left_idx
        , .child2Partials         = (int) b_right_partials
        , .child2TransitionMatrix = (int) b_right_idx
        };

    //-- Push operations, nodes, and branches into respective vectors
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(b_node_idx);
    this->b_branch_lengths.push_back(b_branch_length);

#if !defined ( RB_USE_EIGEN3 )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    // Compute the transition probability matrix in revbayes
    this->updateTransitionProbabilities( node_index );

    // Set the transition probability matrix in BEAGLE
    const double * b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance->getResourceID()
                             , (int) b_node_idx
                             , b_tp_begin
                             , (double) 1.0
                             );

    // Update partial buffers for the node
    beagleUpdatePartials( this->beagle_instance->getResourceID()
                        , &b_operation
                        , 1
                        , BEAGLE_OP_NONE
                        );
#endif //-- !RB_BEAGLE_EIGEN 
}


//TODO : This should probably never exist.... Why is this here
template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeInternalNodeLikelihood
  ( const TopologyNode &node
  , size_t node_index
  , size_t left
  , size_t right
  , size_t middle
  )
{
    size_t node_idx   = node_index;
    size_t left_idx   = left;
    size_t right_idx  = right;
    size_t middle_idx = middle;
    if ( RbSettings::userSettings().getUseBeagleLikelihoodStoring() == true )
    {
        node_idx   = node_index + this->num_nodes * this->activeLikelihood[node_index];
        left_idx   = left       + this->num_nodes * this->activeLikelihood[left];
        right_idx  = right      + this->num_nodes * this->activeLikelihood[right];
        middle_idx = middle     + this->num_nodes * this->activeLikelihood[middle];
    }
    
    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials   = (node.getChild(0).isTip()) ? left   : left_idx;
    size_t right_partials  = (node.getChild(1).isTip()) ? right  : right_idx;
    size_t middle_partials = (node.getChild(2).isTip()) ? middle : middle_idx;

    double branch_length = this->calculateBranchLength(node, node_index);

    //-- TODO : Check which operation for middle
    BeagleOperation b_operation =
        { .destinationPartials    = (int) node_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    //-- push operations, nodes, and branches into respective vectors
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(node_idx);
    this->b_branch_lengths.push_back(branch_length);

#if !defined ( RB_BEAGLE_EIGEN )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    this->updateTransitionProbabilities( node_index );

    const double * b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance->getResourceID()
                               //, (int) node_idx
                             , (int) node_index
                             , b_tp_begin
                             , (double) 1.0
                             );

    beagleUpdatePartials( this->beagle_instance->getResourceID()
                        , &b_operation
                        , 1
                        , BEAGLE_OP_NONE
                        );
#endif //-- RB_BEAGLE_EIGEN 

}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeTipLikelihood
  ( const TopologyNode &node
  , size_t node_index
  )
{
    size_t b_node_idx      = node_index + this->num_nodes * this->activeLikelihood[node_index];
    double b_branch_length = this->calculateBranchLength(node, node_index);

    this->b_node_indices.push_back(b_node_idx);
    this->b_branch_lengths.push_back(b_branch_length);

#if !defined ( RB_USE_EIGEN3 )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    // Compute the transition probability matrix in revbayes
    this->updateTransitionProbabilities( node_index );

    // Set the transition probability matrix in BEAGLE
    const double * b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance->getResourceID()
                             , (int) b_node_idx
                             , b_tp_begin
                             , (double) 1.0
                             );
#endif //-- !RB_BEAGLE_EIGEN 
}


#endif

