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
#undef RB_BEAGLE_DEBUG


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
            std::vector<BeagleOperation>   b_ops;
            std::vector<int>               b_node_indices;
            std::vector<double>            b_branch_lengths;

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
    this->ln_beagle_probability = 0.0;

    for ( auto b_instance : this->beagle_instances ) {
        this->ln_beagle_probability += b_instance->getStoredLnLikelihood();
    }

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

    size_t b_model   = 0;
    size_t num_taxa  = (this->num_nodes + 1) / 2;
    
    size_t root_idx  = root  + this->num_nodes * this->activeLikelihood[root];
    size_t left_idx  = left  + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right + this->num_nodes * this->activeLikelihood[right];

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
    int b_ret_code;

    //-- Update rates across sites 
    this->updateBeagleSiteRates(); // TODO should be in abstract class

#if defined ( RB_BEAGLE_DEBUG )
    ss << "updated site rates" << std::endl;
#endif /* RB_BEAGLE_DEBUG */

#if defined ( RB_BEAGLE_EIGEN )
    this->updateBeagleEigensystems(); // TODO should be in abstract class

    for ( size_t i = 0; i < this->beagle_instances.size(); ++i ) {
        //b_model = this->active_eigen_system[i];
        //b_model = i + 0 * this->active_eigen_system[i] * this->num_site_mixtures;
	b_model = i;

        //-- Update transition matrix for model i.
        b_ret_code = beagleUpdateTransitionMatrices( this->beagle_instances[i]->getResourceID(),
						     b_model,
						     &this->b_node_indices[0],
						     NULL,
						     NULL,
						     &this->b_branch_lengths[0],
						     this->b_branch_lengths.size() );
        if ( b_ret_code != 0 ) {
            throw RbException("Could not update transition matrix for model '" +
                              std::to_string(i) + "'. " +
                              BeagleUtilities::printErrorCode(b_ret_code));
        }

        //-- Calculate and update all partial likelihood buffers
        b_ret_code = beagleUpdatePartials( this->beagle_instances[i]->getResourceID(),
					   &this->b_ops[0],
					   this->b_ops.size(),
					   BEAGLE_OP_NONE );
    }
#else
    for ( size_t i = 0; i < this->beagle_instances.size(); ++i ) {
	b_ret_code = beagleUpdatePartials( this->beagle_instances[i]->getResourceID(),
					   &b_operation,
					   1,
					   BEAGLE_OP_NONE );
    };
#endif /* RB_BEAGLE_EIGEN */

#if defined ( RB_BEAGLE_DEBUG )
    ss << "updated partials" << std::endl;
#endif /* RB_BEAGLE_DEBUG */
    
    if ( b_ret_code != 0 ) {
        throw RbException( "Could not update partials for models '" +
			   BeagleUtilities::printErrorCode(b_ret_code));
    }

    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();

    for ( size_t i = 0; i < this->beagle_instances.size(); ++i ) {
        //-- Calclulate the lnLikelihood of the model
        b_ret_code = beagleCalculateEdgeLogLikelihoods( this->beagle_instances[i]->getResourceID(),
							&b_parentBufferIndices,
							&b_childBufferIndices,
							&b_probabilityIndices,
							b_firstDerivativeIndices,
							b_secondDerivativeIndices,
							&b_categoryWeightsIndices,
							&b_stateFrequenciesIndices,
							&b_cumulativeScaleIndices,
							b_count,
							&b_outSumLogLikelihood,
							b_outSumFirstDerivative,
							b_outSumSecondDerivative );
        if (b_ret_code != 0) {
	    throw RbException("Could not calculate edge log likelihood for models '" +
			      BeagleUtilities::printErrorCode(b_ret_code));
        }

        //this->ln_beagle_probability = b_outSumLogLikelihood;
	this->beagle_instances[i]->setStoredLnLikelihood(b_outSumLogLikelihood);
    }

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

    //-- Return codes for BEAGLE operations.
    size_t b_ret_code;

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
    //-- And update the respective ASRV BEAGLE buffers.
    this->updateBeagleSiteRates();
#else
    std::vector<int> categoryIndicesASRV;
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
    double   b_outSumLogLikelihood     = std::numeric_limits<double>::min();
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

    
#if defined ( RB_USE_EIGEN3 )
    //-- Update Eigensystem BEAGLE buffers
    this->updateBeagleEigensystems();  // TODO should be in abstract class

    for ( size_t i = 0; i < this->beagle_instances.size(); ++i )
    {
        //b_model_idx = i + this->active_eigen_system[i] * this->num_site_mixtures;
        b_model_idx = i;

        //-- Update all transition matrices for model i.
        b_ret_code = beagleUpdateTransitionMatrices( this->beagle_instances[i]->getResourceID(),
						     b_model_idx,
						     &this->b_node_indices[0],
						     NULL,
						     NULL,
						     &this->b_branch_lengths[0],
						     this->b_branch_lengths.size() );
        if ( b_ret_code != 0 )
        {
            throw RbException( "Could not update transition matrix for model '"
                             + std::to_string(i) + "'. "
                             + BeagleUtilities::printErrorCode(b_ret_code));
        }

        //-- Calculate and update all partial likelihood buffers
        b_ret_code = beagleUpdatePartials( this->beagle_instances[i]->getResourceID(),
					   &this->b_ops[0],
					   this->b_ops.size(),
					   BEAGLE_OP_NONE );
    }
#else
    std::vector<std::vector<double>> model_pi_vectors;
    this->getRootFrequencies(model_pi_vectors);
    int b_stateFrequenciesIndex = 0;
    
    for ( size_t i = 0; i < this->beagle_instances.size(); ++i )
    {
	b_stateFrequenciesIndex = i;
        std::vector<double> b_inStateFrequencies = model_pi_vectors[i];
    
        beagleSetStateFrequencies( this->beagle_instances[i]->getResourceID()
                                 , b_stateFrequenciesIndex
                                 , &b_inStateFrequencies[0]
                                 );
    
        b_ret_code = beagleUpdatePartials( this->beagle_instances[i]->getResourceID(),
    				       &b_operation,
    				       1,
    				       BEAGLE_OP_NONE );
    }
#endif //-- RB_BEAGLE_EIGEN

    if ( b_ret_code != 0 ) {
        throw RbException("Could not update partials for models '" +
			  BeagleUtilities::printErrorCode(b_ret_code));
    }

    //-- Reset the beagle operations queues 
    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();

    for ( size_t i = 0; i < this->beagle_instances.size(); ++i )
    {
        //-- Calclulate the lnLikelihood of the model
        b_ret_code =
            beagleCalculateEdgeLogLikelihoods(
                this->beagle_instances[i]->getResourceID(),
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
        if (b_ret_code != 0) {
            throw RbException("Could not calculate edge log likelihood for models '" +
                              BeagleUtilities::printErrorCode(b_ret_code));
        }
        
        this->beagle_instances[i]->setStoredLnLikelihood(b_outSumLogLikelihood);
    }

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

    // Set the transition probability matrices for each model in BEAGLE
    for ( size_t i = 0; i < this->beagle_instances.size(); ++i )
    {
      const double *b_tp_begin = this->transition_prob_matrices[i].theMatrix;
      beagleSetTransitionMatrix( this->beagle_instances[i]->getResourceID(),
                                (int) b_node_idx,
				 b_tp_begin,
				 (double) 1.0 );

      // Update partial buffers for the node
      beagleUpdatePartials( this->beagle_instances[i]->getResourceID(),
			    &b_operation,
			    1,
			    BEAGLE_OP_NONE );
    }
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
    size_t node_idx   = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t left_idx   = left       + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx  = right      + this->num_nodes * this->activeLikelihood[right];
    size_t middle_idx = middle     + this->num_nodes * this->activeLikelihood[middle];
    
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

    // Compute and update the transition probability matrices in revbayes
    this->updateTransitionProbabilities( node_index );

    // Update all model transition matrices and partials
    for ( size_t i = 0; i < this->beagle_instances.size(); ++i ) {
	const double *b_tp_begin = this->transition_prob_matrices[i].theMatrix;

	beagleSetTransitionMatrix( this->beagle_instances[i]->getResourceID(),
				   (int) node_index,
				   b_tp_begin,
				   (double) 1.0 );

	beagleUpdatePartials( this->beagle_instances[i]->getResourceID(),
			      &b_operation,
			      1,
			      BEAGLE_OP_NONE );
    }
#endif //-- RB_BEAGLE_EIGEN 

}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeTipLikelihood
( const TopologyNode &node,
  size_t node_index
)
{
    size_t b_node_idx      = node_index + this->num_nodes * this->activeLikelihood[node_index];
    double b_branch_length = this->calculateBranchLength(node, node_index);

    this->b_node_indices.push_back(b_node_idx);
    this->b_branch_lengths.push_back(b_branch_length);

#if !defined ( RB_USE_EIGEN3 )
    //-- If we are not using the eigensystem, we will need to update and set the
    //   transition probability matrices.

    // Compute and update the transition probability matrices in revbayes
    this->updateTransitionProbabilities( node_index );

    // Update all model transition matrices
    for ( size_t i = 0; i < this->beagle_instances.size(); ++i )
    {
	// Set the transition probability matrix in BEAGLE
	const double * b_tp_begin = this->transition_prob_matrices[i].theMatrix;

	beagleSetTransitionMatrix( this->beagle_instances[i]->getResourceID(),
				   (int) b_node_idx,
				   b_tp_begin,
				   (double) 1.0 );
    }
#endif //-- !RB_BEAGLE_EIGEN 
}


#endif

