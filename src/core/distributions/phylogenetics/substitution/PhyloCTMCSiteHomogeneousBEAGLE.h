#ifndef PhyloCTMCSiteHomogeneousBEAGLE_H
#define PhyloCTMCSiteHomogeneousBEAGLE_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"


#define RB_BEAGLE_DEBUG
//#define RB_BEAGLE_DEBUG_EIGEN
//#define RB_BEAGLE_DEBUG_TIP
//#define RB_BEAGLE_DEBUG_BRANCH

// #define RB_BEAGLE_EIGEN

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

    size_t root_idx  = root + this->num_nodes * this->activeLikelihood[root];
    //this->b_node_indices.push_back(root_idx); //-- TESTING! -- not in original

    size_t left_idx  = left   + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right  + this->num_nodes * this->activeLikelihood[right];

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
            beagleUpdateTransitionMatrices( this->beagle_instance
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
                             + this->parseBeagleReturnCode(b_code_update_transitions));
        }
    }

    //-- Calculate and update all partial likelihood buffers
    b_code_update_partials = beagleUpdatePartials( this->beagle_instance
                                                 , &this->b_ops[0]
                                                 , this->b_ops.size()
                                                 , BEAGLE_OP_NONE
                                                 );
#else
    b_code_update_partials = beagleUpdatePartials( this->beagle_instance
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
	  		             + this->parseBeagleReturnCode(b_code_update_partials));
	}

    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();

    //-- Calclulate the lnLikelihood of the model
    b_code_calc_edges =
        beagleCalculateEdgeLogLikelihoods( this->beagle_instance
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
                          + this->parseBeagleReturnCode(b_code_calc_edges));
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

    //-- TODO - maybe active indices are out of phase??
    //for ( size_t i = 0; i < this->activeLikelihood.size(); ++i ) {
    //    this->activeLikelihood[i] = (this->activeLikelihood[i] == 0 ? 1 : 0);
    //}
    
    //-- Calculate the node indices accounting for active/inactive offests.
    size_t root_idx  = root   + this->num_nodes * this->activeLikelihood[root];
    size_t mid_idx   = middle + this->num_nodes * this->activeLikelihood[middle];
    size_t left_idx  = left   + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right  + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once in BEAGLE, so we cant have active/inactive offests.
    size_t mid_partials   = (middle < num_taxa) ? middle : mid_idx;
    size_t left_partials  = (left   < num_taxa) ? left   : left_idx;
    size_t right_partials = (right  < num_taxa) ? right  : right_idx;

    //-- Set the ASRV category for each site. Since we do not allow for partitions, this is always 0.
    std::vector<int> categoryIndicesASRV(this->pattern_block_size, 0);
    //-- And Update the respective ASRV BEAGLE buffers.
    this->updateBeagleSiteRates();
    
    //-- Configure BEAGLE model parameters.
    int      b_parentBufferIndices     = (int) root_idx;
    int      b_childBufferIndices      = (int) mid_partials;
    int      b_probabilityIndices      = (int) mid_idx;
    int *    b_firstDerivativeIndices  = NULL;
    int *    b_secondDerivativeIndices = NULL;
    int *    b_categoryWeightsIndices = &categoryIndicesASRV[0];
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
        , .child1TransitionMatrix = (int) left_partials
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_partials
        };

    //-- Push operation and root index onto respective vectors to prepare for likelihood calculation.
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(root_idx);

    //-- Return codes for BEAGLE operations.
    int b_code_update_transitions;
    int b_code_update_partials;
    int b_code_calc_edges;

    //-- TODO - uncomment when seteigensystem is fixed to only homogeneous models
    //std::vector<std::vector<double>> model_pi_vectors;
    //this->getRootFrequencies(model_pi_vectors);
    //int b_stateFrequenciesIndex = 0;
    //std::vector<double> b_inStateFrequencies = model_pi_vectors[0];

    //beagleSetStateFrequencies( this->beagle_instance
    //                         , b_stateFrequenciesIndex
    //                         , &b_inStateFrequencies[0]
    //                         );

    
#if defined ( RB_BEAGLE_EIGEN )
    //-- Update Eigensystem BEAGLE buffers
    this->updateBeagleEigensystem();


    for ( size_t i = 0; i < b_num_models; ++i )
    {
        //b_model_idx = 0; //this->active_eigen_system[i];
        b_model_idx = i + this->active_eigen_system[i] * this->num_site_mixtures;

        //-- Update all transition matrices for model i.
        b_code_update_transitions =
            beagleUpdateTransitionMatrices( this->beagle_instance
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
                             + this->parseBeagleReturnCode(b_code_update_transitions));
        }
    }

#if defined ( RB_BEAGLE_DEBUG )
    RBOUT("\n");
    RBOUT("b_node_indices len: " + std::to_string(this->b_node_indices.size()));
    RBOUT("b_ops len: " + std::to_string(this->b_ops.size()));
    RBOUT("\n");
#endif //-- RB_BEAGLE_DEBUG
    //-- Calculate and update all partial likelihood buffers
    b_code_update_partials = beagleUpdatePartials( this->beagle_instance
                                                 , &this->b_ops[0]
                                                 , this->b_ops.size()
                                                 , BEAGLE_OP_NONE
                                                 );


#if defined ( RB_BEAGLE_DEBUG )
    //-- Debug transitions and partial buffers (verbose! only use with small number of sites and taxa)
    std::stringstream debug_ss;
    double * t = (double *) malloc(16 * sizeof(double));
    for ( size_t i = 0; i < this->b_node_indices.size() - 1; ++i ) { //-- we dont care about root in unrooted case
        if ( this->b_node_indices[i] >= num_taxa ) { 
            beagleGetTransitionMatrix(this->beagle_instance, this->b_node_indices[i], t);
            debug_ss << "Transition for node " + std::to_string(this->b_node_indices[i]) << ":\t";
            for ( size_t j = 0; j < this->num_chars * this->num_chars; ++j )
            {
                if ( (j % this->num_chars) == 0) { debug_ss << "\n\t\t"; };
                debug_ss << std::fixed << std::setw(8) << std::setprecision(4) << t[j] << " ";
            }
            debug_ss << std::endl;
        }
    }
    free(t);
    debug_ss << std::endl;

    double * partial = (double *) malloc(this->num_chars*this->num_sites*sizeof(double));
    for ( size_t i = 0; i < this->b_node_indices.size(); ++i ) {
        if (this->b_node_indices[i] >= num_taxa) {
            beagleGetPartials(this->beagle_instance, this->b_node_indices[i], NULL, partial);
            debug_ss << "Partial buffer for node " + std::to_string(this->b_node_indices[i]) << ":\t";
            for (size_t j = 0; j < this->num_chars*this->num_sites; ++j ) {
                if ( (j % this->num_chars) == 0) { debug_ss << "\n\t\t"; };
                debug_ss << std::fixed << std::setw(8) << std::setprecision(12)
                         << partial[j] << " ";
            }
            debug_ss << std::endl;
        }
    }

//    for ( size_t i = 0; i < this->b_node_indices.size(); ++i ) {
//        if (this->b_node_indices[i] >= num_taxa) {
//            beagleGetPartials(this->beagle_instance, this->b_node_indices[i], NULL, partial);
//            double lnL_total = 0;
//            //for ( size_t j = 0; j < this->num_patterns; ++j ) {
//            for ( size_t j = 0; j < this->num_sites; ++j ) {
//                double site_L = 0;
//                for ( size_t k = 0; k < this->num_chars; ++k ) {
//                    site_L += this->pattern_counts[j] * model_pi_vectors[0][k] * partial[k+j*this->num_chars];     
//                }
//                if ( std::isnan(site_L) || site_L <= 0 || site_L > 1 ) {
//                    debug_ss << "Error: site " << std::to_string(j) << " val = " << std::to_string(site_L) << std::endl;
//
//                    for ( size_t k = 0; k < this->num_chars; ++k ) {
//                        debug_ss << "\tpat_count " << std::to_string(this->pattern_counts[j])
//                                 << "\t" << std::to_string(model_pi_vectors[0][k])
//                                 << "\t" << std::to_string(partial[k+j*this->num_chars]) << std::endl;     
//                    }
//                }
//                lnL_total += std::log(site_L);
//            }
//            debug_ss << "lnL node " << std::to_string(this->b_node_indices[i]) << ": " << std::to_string(lnL_total) << std::endl;
//
//            double * scale = (double*)malloc(this->num_sites*sizeof(double));
//            //beagleGetScaleFactors(this->beagle_instance, this->b_node_indices[i], scale);
//            beagleGetScaleFactors(this->beagle_instance, this->b_node_indices[i] - (1+this->activeLikelihood[this->b_node_indices[i]])*num_taxa, scale);
//            debug_ss <<  "scale: ";
//            for ( size_t i = 0; i < this->num_sites; ++i ) {
//                debug_ss << " " << std::to_string(scale[i]);
//            }
//            free(scale);
//            debug_ss << std::endl;
//
//        }
//    }

    free(partial);

    std::vector<std::vector<double>> model_pi_vectors;
    this->getRootFrequencies(model_pi_vectors);
    debug_ss << std::endl << "Current Site Patterns..." << std::endl;
    debug_ss << "\tCounts : ";
    for ( auto x : this->b_inPatternWeights ) {
        debug_ss << std::to_string(x) << " ";
    }
    debug_ss << std::endl << std::endl;
    debug_ss << "\tWeights : ";
    for ( auto x : this->pattern_counts ) {
        debug_ss << std::to_string(x) << " ";
    }
    debug_ss << std::endl << std::endl;

    RBOUT(debug_ss.str());
#endif //-- RB_BEAGLE_DEBUG
    
#else
    // TODO - remove this from set eigensystem as only homogeneous models supported
    std::vector<std::vector<double>> model_pi_vectors;
    this->getRootFrequencies(model_pi_vectors);
    int b_stateFrequenciesIndex = 0;
    std::vector<double> b_inStateFrequencies = model_pi_vectors[0];

    beagleSetStateFrequencies( this->beagle_instance
                             , b_stateFrequenciesIndex
                             , &b_inStateFrequencies[0]
                             );

    b_code_update_partials = beagleUpdatePartials( this->beagle_instance
                                                 , &b_operation
                                                 , 1
                                                 , BEAGLE_OP_NONE
                                                 );
#endif //-- RB_BEAGLE_EIGEN

    if ( b_code_update_partials != 0 )
	{
        throw RbException( "Could not update partials for models '"
	  		             + this->parseBeagleReturnCode(b_code_update_partials));
	}

#if defined ( RB_BEAGLE_DEBUG )
    ss << std::endl << "Branch Lengths:" << std::endl << "\t";
    for (size_t i = 0; i < this->b_branch_lengths.size(); ++i) {
        ss << std::to_string(this->b_branch_lengths[i]) << " ";
    }
    ss << std::endl;

    ss << std::endl << "Node Indices:" << std::endl << "\t";
    for (size_t i = 0; i < this->b_node_indices.size(); ++i) {
        ss << std::to_string(this->b_node_indices[i]) << " ";
    }
    ss << std::endl;
#endif //-- RB_BEAGLE_DEBUG

    //-- Reset the beagle operations queues 
    this->b_ops.clear();
    this->b_branch_lengths.clear();
    this->b_node_indices.clear();

    //-- Calclulate the lnLikelihood of the model
    b_code_calc_edges =
        beagleCalculateEdgeLogLikelihoods( this->beagle_instance
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
    //-- TODO print val
    if ( b_code_calc_edges != 0 )
    {
	    throw RbException("Could not calculate edge log likelihood for models '"
                          + this->parseBeagleReturnCode(b_code_calc_edges));
    }

    this->ln_beagle_probability = b_outSumLogLikelihood;
    
#if defined ( RB_BEAGLE_DEBUG )
    if ( b_outSumLogLikelihood != -std::numeric_limits<double>::infinity() ) {
        ss << std::endl;
        //ss << "Calculated lnLikelihood: " << std::to_string(b_outSumLogLikelihood);
    }
    else {
        ss << std::endl;
        ss << "Calculated lnLikelihood: -inf";
    }
    ss << std::endl;
    ss << "End of BEAGLE calculation for instance " << std::to_string(this->beagle_instance);
    ss << std::endl;
    RBOUT(ss.str());
#endif /* RB_BEAGLE_DEBUG */

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
    // compute the transition probability matrix
    this->updateTransitionProbabilities( node_index );

    size_t node_idx  = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t left_idx  = left       + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right      + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials  = (node.getChild(0).isTip()) ? left  : left_idx;
    size_t right_partials = (node.getChild(1).isTip()) ? right : right_idx;

    // compute the branch length
    double branch_length = this->calculateBranchLength(node, node_index);

    BeagleOperation b_operation =
        { .destinationPartials    = (int) node_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_partials
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_partials
        };

#if !defined ( RB_BEAGLE_EIGEN )
    const double * b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance
                             , (int) node_idx
                             , b_tp_begin
                             , (double) 1.0
                             );

    beagleUpdatePartials( this->beagle_instance
                        , &b_operation
                        , 1
                        , BEAGLE_OP_NONE
                        );
#endif //-- RB_BEAGLE_EIGEN 

    //-- push operations, nodes, and branches into respective vectors
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(node_idx);
    this->b_branch_lengths.push_back(branch_length);
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
    // compute the transition probabilities
    this->updateTransitionProbabilities( node_index );

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

#if !defined ( RB_BEAGLE_EIGEN )
    const double * b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance
                             , (int) node_index
                             , b_tp_begin
                             , (double) 1.0
                             );

    beagleUpdatePartials( this->beagle_instance
                        , &b_operation
                        , 1
                        , BEAGLE_OP_NONE
                        );
#endif //-- RB_BEAGLE_EIGEN 

    //-- push operations, nodes, and branches into respective vectors
    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(node_idx);
    this->b_branch_lengths.push_back(branch_length);
}


template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeTipLikelihood
  ( const TopologyNode &node
  , size_t node_index
  )
{
    // compute the transition probabilities
    this->updateTransitionProbabilities( node_index );

    //size_t node_idx = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t node_idx = node_index;

    double branch_length = this->calculateBranchLength(node, node_index);

#if !defined ( RB_BEAGLE_EIGEN )
    const double* b_tp_begin = this->transition_prob_matrices[0].theMatrix;
    beagleSetTransitionMatrix( this->beagle_instance
                             , node_idx
                             , b_tp_begin
                             , (double) 1.0
                             );
#endif //-- RB_BEAGLE_EIGEN 

    this->b_branch_lengths.push_back(branch_length);
    this->b_node_indices.push_back(node_idx);
}


#endif

