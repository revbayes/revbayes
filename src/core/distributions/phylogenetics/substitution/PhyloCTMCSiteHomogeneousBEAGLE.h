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
 *       and adjust accordingly. For now, only use ASRV when compiling with Eigen3 support.
 *
 *     - Configure for usage with MPI.
 *
 *     - Configure for usage with mixture models.
 *
 */

#ifndef PhyloCTMCSiteHomogeneousBEAGLE_H
#define PhyloCTMCSiteHomogeneousBEAGLE_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"


#define RB_BEAGLE_DEBUG
#define RB_BEAGLE_DEBUG_TIP
#define RB_BEAGLE_INFO
//#undef RB_BEAGLE_DEBUG


namespace RevBayesCore
{

    template<class charType>
    class PhyloCTMCSiteHomogeneousBEAGLE : public AbstractPhyloCTMCSiteHomogeneous<charType>
    {

        public:

            //----====  Constructors  ====----

            //-- Default constructor
            PhyloCTMCSiteHomogeneousBEAGLE ( const TypedDagNode<Tree>* t, size_t nChars, bool c, size_t nSites, bool amb, bool internal, bool gapmatch );

            //-- Destructor
            virtual ~PhyloCTMCSiteHomogeneousBEAGLE ( void );
            
        
            //-- Clone constructor
            PhyloCTMCSiteHomogeneousBEAGLE*         clone ( void ) const;


        protected:

            virtual double                          sumRootLikelihood( void );                                                                          //!< Return the computed likelihood.
            virtual void                            computeRootLikelihood( size_t root, size_t l, size_t r );                                           //!< BEAGLE compute lnLikelihood of a rooted tree.
            virtual void                            computeRootLikelihood( size_t root, size_t l, size_t r, size_t m);                                  //!< BEAGLE compute lnLikelihood of an unrooted tree.
            virtual void                            computeInternalNodeLikelihood ( const TopologyNode &n, size_t nIdx, size_t l, size_t r );           //!< Collect a BEAGLE operation for an internal node into the computation queue.
            virtual void                            computeInternalNodeLikelihood ( const TopologyNode &n, size_t nIdx, size_t l, size_t r, size_t m);  //!< Collect a BEAGLE operation for an internal node into the computation queue.
            virtual void                            computeTipLikelihood( const TopologyNode &node, size_t nIdx);                                       //!< Collect a BEAGLE operation for a leaf node into the computation queue.

        private:

            std::vector<BeagleOperation>            b_ops;
            std::vector<int>                        b_node_indices;
            std::vector<int>                        b_scale_indices;
            std::vector<double>                     b_branch_lengths;

            double                                  ln_beagle_probability;                  //!< Accumulate BEAGLE lnLikelihood across all models.
        
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
    this->ln_beagle_probability = this->beagle_instance->getStoredLnLikelihood();

    //-- We have already computed the ln_beagle_probability, so just return it
    return this->ln_beagle_probability;
}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeRootLikelihood( size_t root, size_t left, size_t right )
{

    //-- Return codes for BEAGLE operations.
    size_t b_ret_code;

    //-- Get the number of taxa in the tree
    size_t num_taxa = (this->num_nodes + 2) / 2;

    //-- Determine the number of unique eigensystems we will have.
    size_t b_num_models = this->num_site_mixtures - this->num_site_rates;
    if ( b_num_models < 1 ) { b_num_models = 1; }
    
    
    size_t root_idx  = root  + this->num_nodes * this->active_likelihood[root];
    size_t left_idx  = left  + this->num_nodes * this->active_likelihood[left];
    size_t right_idx = right + this->num_nodes * this->active_likelihood[right];

    //-- Tips are actually just stored once, so we dont need offsets.
    size_t left_partials  = (left  < num_taxa) ? left  : left_idx;
    size_t right_partials = (right < num_taxa) ? right : right_idx;

#if defined ( RB_USE_EIGEN3 )
    //-- Set the ASRV category for each site. Since we do not allow for partitions, this is always 0.
    std::vector<int> categoryIndicesASRV(this->pattern_block_size, 0);
    //-- And update the respective ASRV BEAGLE buffers.
    this->updateBeagleSiteRates();
#else
    std::vector<int> categoryIndicesASRV;
#endif
    
    int scaler_index_read  = BEAGLE_OP_NONE;
    int scaler_index_write = (int) root_idx;
    
    //-- Push the last operation onto the queue
    BeagleOperation b_operation =
        { .destinationPartials    = (int) root_idx
        , .destinationScaleWrite  = scaler_index_write
        , .destinationScaleRead   = scaler_index_read
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    //-- Push operation and root index onto respective vectors to prepare for likelihood calculation.
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back((int)root_idx);
    this->b_scale_indices.push_back((int)root_idx);

#if defined ( RB_USE_EIGEN3 )
    //-- Update Eigensystem BEAGLE buffers
    this->updateBeagleEigensystems();  // TODO should be in abstract class

    //-- TODO - not working with mixture models yet!!
    if ( this->num_mixtures > 1 )
    {
        throw RbException( "Mixture models not supported when using BEAGLE! Aborting.'");
        // This is a start to get mixture models working. Right now the `b_category_indices`
        // do not seem to be set correctly. Maybe need to look closer at how to set the
        // `category weights` for BEAGLE that is for the model `categories` and not for
        // the `rates`... Documentation not clear...
        std::vector<int> b_model_indices;
        std::vector<int> b_category_indices;
        for (int i = 0; i < this->num_mixtures; ++i)
        {
            b_model_indices.push_back(i);
            b_category_indices.push_back(i);
        }
        b_ret_code = beagleUpdateTransitionMatricesWithMultipleModels(
            this->beagle_instance->getResourceID(),
            &b_model_indices[0],
            &b_category_indices[0],
            &this->b_node_indices[0],
            NULL,
            NULL,
            &this->b_branch_lengths[0],
            this->b_branch_lengths.size()
            );
    }
    else
    {
        b_ret_code = beagleUpdateTransitionMatrices(
            this->beagle_instance->getResourceID(),
            0,
            &this->b_node_indices[0],
            NULL,
            NULL,
            &this->b_branch_lengths[0],
            this->b_branch_lengths.size()
            );
    }
    
    if (b_ret_code != 0)
    {
        throw RbException( "Could not update transition matrix. "
                            + BeagleUtilities::printErrorCode(b_ret_code));
    }

    //-- Calculate and update all partial likelihood buffers
    b_ret_code = beagleUpdatePartials( this->beagle_instance->getResourceID(),
                                       &this->b_ops[0],
                                        this->b_ops.size(),
                                        BEAGLE_OP_NONE );
        
#else
    std::vector<std::vector<double>> model_pi_vectors;
    this->getRootFrequencies(model_pi_vectors);
    int b_stateFrequenciesIndex = 0;
    
    b_stateFrequenciesIndex = 0;
    std::vector<double> b_inStateFrequencies = model_pi_vectors[0];
    
    beagleSetStateFrequencies( this->beagle_instance->getResourceID(),
                               b_stateFrequenciesIndex,
                               &b_inStateFrequencies[0] );
    
    b_ret_code = beagleUpdatePartials( this->beagle_instance->getResourceID(),
                                       &b_operation,
                                       1,
                                       BEAGLE_OP_NONE );
#endif /* RB_USE_EIGEN3 */
    
    if ( b_ret_code != 0 )
    {
        throw RbException( "Could not update partials for models '" +
			   BeagleUtilities::printErrorCode((int)b_ret_code));
    }

    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();
//    this->b_scale_indices.clear();
    

    //-- BEAGLE model parameters.
    int     b_parentBufferIndices     = (int) root_idx;
    int*    b_categoryWeightsIndices  = &categoryIndicesASRV[0];
    int     b_stateFrequenciesIndices = 0; //(int) model;  //0;
    int     b_cumulativeScaleIndices  = (int) 2*this->num_nodes+this->active_likelihood[root];
    int     b_count                   = 1;
    double  b_outSumLogLikelihood     = std::numeric_limits<double>::min();
    
    beagleResetScaleFactors(this->beagle_instance->getResourceID(), b_cumulativeScaleIndices);
    beagleAccumulateScaleFactors(this->beagle_instance->getResourceID(), &b_scale_indices[0], b_scale_indices.size(),
                                 b_cumulativeScaleIndices);
    
    this->b_scale_indices.clear();

    b_ret_code = beagleCalculateRootLogLikelihoods( this->beagle_instance->getResourceID(),
                                                   &b_parentBufferIndices,
                                                   b_categoryWeightsIndices,
                                                   &b_stateFrequenciesIndices,
                                                   &b_cumulativeScaleIndices,
                                                   b_count,
                                                   &b_outSumLogLikelihood );
    
    if (b_ret_code != 0)
    {
        throw RbException("Could not calculate edge log likelihood for models '" +
                BeagleUtilities::printErrorCode(b_ret_code));
    }

    //this->ln_beagle_probability = b_outSumLogLikelihood;
    this->beagle_instance->setStoredLnLikelihood(b_outSumLogLikelihood);
    
}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeRootLikelihood( size_t root, size_t left, size_t right, size_t middle )
{
    size_t b_model_idx;

    //-- Return codes for BEAGLE operations.
    size_t b_ret_code;

    //-- Get the number of taxa in the tree
    size_t num_taxa = (this->num_nodes + 2) / 2;

    //-- Determine the number of unique eigensystems we will have.
    size_t b_num_models = this->num_site_mixtures - this->num_site_rates;
    if ( b_num_models < 1 ) { b_num_models = 1; }
    
    //-- Calculate the node indices accounting for active/inactive offests.
    size_t root_idx  = root   + this->num_nodes * this->active_likelihood[root];
    size_t mid_idx   = middle + this->num_nodes * this->active_likelihood[middle];
    size_t left_idx  = left   + this->num_nodes * this->active_likelihood[left];
    size_t right_idx = right  + this->num_nodes * this->active_likelihood[right];

    //-- Tips are actually just stored once in BEAGLE, so we cant have active/inactive offests.
    size_t mid_partials   = (middle < num_taxa) ? middle : mid_idx;
    size_t left_partials  = (left   < num_taxa) ? left   : left_idx;
    size_t right_partials = (right  < num_taxa) ? right  : right_idx;

#if defined ( RB_USE_EIGEN3 )
    //-- Set the ASRV category for each site. Since we do not allow for partitions, this is always 0.
    std::vector<int> categoryIndicesASRV(this->pattern_block_size, 0);
    //-- And update the respective ASRV BEAGLE buffers.
    this->updateBeagleSiteRates();
#else
    std::vector<int> categoryIndicesASRV;
#endif

//    int scaler_index_read  = BEAGLE_OP_NONE;
//    int scaler_index_write = BEAGLE_OP_NONE;
    
    int scaler_index_read  = BEAGLE_OP_NONE;
    int scaler_index_write = (int) root_idx;
    
    //-- Create BEAGLE operation.
    BeagleOperation b_operation =
        { .destinationPartials    = (int) root_idx
        , .destinationScaleWrite  = scaler_index_write
        , .destinationScaleRead   = scaler_index_read
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    //-- Push operation and root index onto respective vectors to prepare for likelihood calculation.
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back((int)root_idx);
    this->b_scale_indices.push_back((int)root_idx);
    
#if defined ( RB_USE_EIGEN3 )
    //-- Update Eigensystem BEAGLE buffers
    this->updateBeagleEigensystems();  // TODO should be in abstract class

    b_model_idx = 0;

    //-- TODO - not working with mixture models yet!!
    if ( this->num_mixtures > 1 )
    {
        throw RbException( "Mixture models not supported when using BEAGLE! Aborting.'");
        // This is a start to get mixture models working. Right now the `b_category_indices`
        // do not seem to be set correctly. Maybe need to look closer at how to set the
        // `category weights` for BEAGLE that is for the model `categories` and not for
        // the `rates`... Documentation not clear...
        std::vector<int> b_model_indices;
        std::vector<int> b_category_indices;
        for (int i = 0; i < this->num_mixtures; ++i)
        {
            b_model_indices.push_back(i);
            b_category_indices.push_back(i);
        }
        b_ret_code = beagleUpdateTransitionMatricesWithMultipleModels(
                                                                      this->beagle_instance->getResourceID(),
                                                                      &b_model_indices[0],
                                                                      &b_category_indices[0],
                                                                      &this->b_node_indices[0],
                                                                      NULL,
                                                                      NULL,
                                                                      &this->b_branch_lengths[0],
                                                                      this->b_branch_lengths.size()
                                                                      );
    }
    else
    {
            b_ret_code = beagleUpdateTransitionMatrices(
                                                        this->beagle_instance->getResourceID(),
                                                        0,
                                                        &this->b_node_indices[0],
                                                        NULL,
                                                        NULL,
                                                        &this->b_branch_lengths[0],
                                                        this->b_branch_lengths.size()
                                                        );
    }
    
    if (b_ret_code != 0)
    {
            throw RbException( "Could not update transition matrix. "
                               + BeagleUtilities::printErrorCode(b_ret_code));
    }

    //-- Calculate and update all partial likelihood buffers
    int      b_cumulativeScaleIndices  = (int) 2*this->num_nodes+1;
    beagleResetScaleFactors(this->beagle_instance->getResourceID(), b_cumulativeScaleIndices);

    b_ret_code = beagleUpdatePartials( this->beagle_instance->getResourceID(),
                                       &this->b_ops[0],
                                       this->b_ops.size(),
                                       b_cumulativeScaleIndices );
#else
    std::vector< std::vector<double> > model_pi_vectors;
    this->getRootFrequencies(model_pi_vectors);
    int b_stateFrequenciesIndex = 0;
    
    b_stateFrequenciesIndex = 0;
    std::vector<double> b_inStateFrequencies = model_pi_vectors[0];
    
    beagleSetStateFrequencies( this->beagle_instance->getResourceID(),
                               b_stateFrequenciesIndex,
                               &b_inStateFrequencies[0] );
    
    b_ret_code = beagleUpdatePartials( this->beagle_instance->getResourceID(),
                                       &b_operation,
                                       1,
                                       BEAGLE_OP_NONE );
#endif //-- RB_USE_EIGEN3

    if ( b_ret_code != 0 )
    {
        throw RbException("Could not update partials for models '" +
                          BeagleUtilities::printErrorCode(b_ret_code));
    }

    //-- Reset the beagle operations queues 
    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();
    
    //-- Configure BEAGLE model parameters.
    int      b_parentBufferIndices     = (int) root_idx;
    int      b_childBufferIndices      = (int) mid_partials;
    int      b_probabilityIndices      = (int) mid_idx;
    int *    b_firstDerivativeIndices  = NULL;
    int *    b_secondDerivativeIndices = NULL;
    int *    b_categoryWeightsIndices  = &categoryIndicesASRV[0];
    int      b_stateFrequenciesIndices = 0; //(int) model;  //0;
//    int      b_cumulativeScaleIndices  = BEAGLE_OP_NONE;
//    int      b_cumulativeScaleIndices  = (int) 2*this->num_nodes+this->active_likelihood[root];
//    int      b_cumulativeScaleIndices  = (int) 2*this->num_nodes+1;
    int      b_count                   = 1;
    double   b_outSumLogLikelihood     = std::numeric_limits<double>::min();
    double * b_outSumFirstDerivative   = NULL;
    double * b_outSumSecondDerivative  = NULL;
        
//    beagleResetScaleFactors(this->beagle_instance->getResourceID(), b_cumulativeScaleIndices);
//    beagleAccumulateScaleFactors(this->beagle_instance->getResourceID(), &b_scale_indices[0], b_scale_indices.size(),
//                                 b_cumulativeScaleIndices);
    
    this->b_scale_indices.clear();
    
    //-- Calclulate the lnLikelihood of the model
    b_ret_code = beagleCalculateEdgeLogLikelihoods(
                                                   this->beagle_instance->getResourceID(),
                                                   &b_parentBufferIndices,
                                                   &b_childBufferIndices,
                                                   &b_probabilityIndices,
                                                   b_firstDerivativeIndices,
                                                   b_secondDerivativeIndices,
                                                   b_categoryWeightsIndices,
                                                   &b_stateFrequenciesIndices,
                                                   &b_cumulativeScaleIndices,
                                                   b_count,
                                                   &b_outSumLogLikelihood,
                                                   b_outSumFirstDerivative,
                                                   b_outSumSecondDerivative );
    if (b_ret_code != 0)
    {
        throw RbException("Could not calculate edge log likelihood for models '" +
                            BeagleUtilities::printErrorCode(b_ret_code));
    }
    
    this->beagle_instance->setStoredLnLikelihood(b_outSumLogLikelihood);
    
}



template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeInternalNodeLikelihood( const TopologyNode &node, size_t node_index, size_t left, size_t right )
{
    //-- Calculate the node indices accounting for active/inactive offests.
    size_t b_node_idx  = node_index + this->num_nodes * this->active_likelihood[node_index];
    size_t b_left_idx  = left       + this->num_nodes * this->active_likelihood[left];
    size_t b_right_idx = right      + this->num_nodes * this->active_likelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t b_left_partials  = (node.getChild(0).isTip()) ? left  : b_left_idx;
    size_t b_right_partials = (node.getChild(1).isTip()) ? right : b_right_idx;

    // Compute the branch length
    double b_branch_length = this->calculateBranchLength(node, node_index);
    
//    int scaler_index_read  = BEAGLE_OP_NONE;
//    int scaler_index_write = BEAGLE_OP_NONE;
    
    int scaler_index_read  = BEAGLE_OP_NONE;
    int scaler_index_write = (int) b_node_idx;

    // Construct the BEAGLE operation that will be pushed onto the compute queue.
    BeagleOperation b_operation =
        { .destinationPartials    = (int) b_node_idx
        , .destinationScaleWrite  = scaler_index_write
        , .destinationScaleRead   = scaler_index_read
        , .child1Partials         = (int) b_left_partials
        , .child1TransitionMatrix = (int) b_left_idx
        , .child2Partials         = (int) b_right_partials
        , .child2TransitionMatrix = (int) b_right_idx
        };

    //-- Push operations, nodes, and branches into respective vectors
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back((int)b_node_idx);
    this->b_scale_indices.push_back((int)b_node_idx);
    this->b_branch_lengths.push_back(b_branch_length);

#if !defined ( RB_USE_EIGEN3 )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    // Compute the transition probability matrix in revbayes
    this->updateTransitionProbabilities( node_index );

    // Set the transition probability matrices in BEAGLE
    const double *b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance->getResourceID(),
                               (int) b_node_idx,
                               b_tp_begin,
                               1.0 );

    // Update partial buffers for the node
    beagleUpdatePartials( this->beagle_instance->getResourceID(),
                          &b_operation,
                          1,
                          BEAGLE_OP_NONE );
#endif //-- !RB_USE_EIGEN3
}


//TODO : This should probably never exist.... Why is this here
template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeInternalNodeLikelihood( const TopologyNode &node, size_t node_index, size_t left, size_t right, size_t middle )
{
    size_t node_idx   = node_index + this->num_nodes * this->active_likelihood[node_index];
    size_t left_idx   = left       + this->num_nodes * this->active_likelihood[left];
    size_t right_idx  = right      + this->num_nodes * this->active_likelihood[right];
    size_t middle_idx = middle     + this->num_nodes * this->active_likelihood[middle];
    
    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials   = (node.getChild(0).isTip()) ? left   : left_idx;
    size_t right_partials  = (node.getChild(1).isTip()) ? right  : right_idx;
    size_t middle_partials = (node.getChild(2).isTip()) ? middle : middle_idx;

    double branch_length = this->calculateBranchLength(node, node_index);

//    int scaler_index_read  = BEAGLE_OP_NONE;
//    int scaler_index_write = BEAGLE_OP_NONE;
    
    int scaler_index_read  = BEAGLE_OP_NONE;
    int scaler_index_write = (int) node_idx;
    
    //-- TODO : Check which operation for middle
    BeagleOperation b_operation =
        { .destinationPartials    = (int) node_idx
        , .destinationScaleWrite  = scaler_index_write
        , .destinationScaleRead   = scaler_index_read
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    //-- push operations, nodes, and branches into respective vectors
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back((int)node_idx);
    this->b_scale_indices.push_back((int)node_idx);
    this->b_branch_lengths.push_back(branch_length);

#if !defined ( RB_USE_EIGEN3 )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    // Compute and update the transition probability matrices in revbayes
    this->updateTransitionProbabilities( node_index );

    // Update all model transition matrices and partials
	const double *b_tp_begin = this->transition_prob_matrices[0].theMatrix;

	beagleSetTransitionMatrix( this->beagle_instance->getResourceID(),
                              (int) node_index,
                              b_tp_begin,
                              1.0 );

	beagleUpdatePartials( this->beagle_instance->getResourceID(),
                          &b_operation,
                          1,
                          BEAGLE_OP_NONE );
#endif //-- RB_USE_EIGEN3

}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeTipLikelihood( const TopologyNode &node, size_t node_index )
{
    size_t b_node_idx      = node_index + this->num_nodes * this->active_likelihood[node_index];
    double b_branch_length = this->calculateBranchLength(node, node_index);

    this->b_node_indices.push_back((int)b_node_idx);
    this->b_scale_indices.push_back((int)b_node_idx);
    this->b_branch_lengths.push_back(b_branch_length);

#if !defined ( RB_USE_EIGEN3 )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    // Compute and update the transition probability matrices in revbayes
    this->updateTransitionProbabilities( node_index );

    // Update transition matrix
	// Set the transition probability matrix in BEAGLE
	const double * b_tp_begin = this->transition_prob_matrices[0].theMatrix;

	beagleSetTransitionMatrix( this->beagle_instance->getResourceID(),
                               (int) b_node_idx,
                               b_tp_begin,
                               1.0 );
#endif //-- !RB_USE_EIGEN3
}


#endif

