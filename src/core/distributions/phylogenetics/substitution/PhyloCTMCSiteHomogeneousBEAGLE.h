#ifndef PhyloCTMCSiteHomogeneousBEAGLE_H
#define PhyloCTMCSiteHomogeneousBEAGLE_H

#include "AbstractPhyloCTMCSiteHomogeneous.h"


#define RB_BEAGLE_DEBUG
//#define RB_BEAGLE_DEBUG_EIGEN
//#define RB_BEAGLE_DEBUG_TIP
//#define RB_BEAGLE_DEBUG_BRANCH


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


            //----====  Private Methods  ====----
            //-- Calculate the tree branch length between adjacent tree nodes..
            double      calculateBranchLength     ( const TopologyNode &node
                                                  , size_t node_index
                                                  );

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

//--}


//-----------------------------------------------------------------------[ Private Helper Methods ]
//--{1

template<class charType>
double
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::calculateBranchLength ( const TopologyNode &node, size_t node_index )
{
    double branch_len;
    double rate ;

    if ( this->branch_heterogeneous_clock_rates == true )
    {
        rate = this->heterogeneous_clock_rates->getValue()[node_index];
    }
    else if ( this->homogeneous_clock_rate != NULL)
    {
        rate = this->homogeneous_clock_rate->getValue();
    }
    else
    {
        rate = 1.0;
    }

    //-- TODO: check if this should be invariable site...
    //rate /= this->homogeneous_clock_rate->getValue();

    branch_len = rate * node.getBranchLength();
    if ( branch_len < 0 )
    {
      throw RbException("Error : Negative branch length!");
    }

#if defined ( RB_BEAGLE_DEBUG_BRANCH )
    std::stringstream ss;
    ss << "Branch length: " << std::to_string(branch_len) << "\n";
    RBOUT(ss.str());
#endif /* RB_BEAGLE_DEBUG_BRANCH */

    return branch_len;
}

//--}


//----------------------------------------------------------------------[ Likelihood Calculations ]
//--{1

template<class charType>
double
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::sumRootLikelihood (void )
{

    #if defined ( RB_BEAGLE_DEBUG )
        std::stringstream ss;
        ss << std::endl << "SumRootlikelihood = "
           << std::to_string(this->ln_beagle_probability)
           << std::endl << std::endl;
        RBOUT(ss.str());
    #endif /* RB_BEAGLE_DEBUG */

    return this->ln_beagle_probability;
}



template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeRootLikelihood ( size_t root, size_t left, size_t right )
{
    //-- TODO : Calculate the lnLikelihood for rooted trees. Should be able to copy the code for the
    //          unroooted case, but must change the child indexes/partials buffers.
}



template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeRootLikelihood ( size_t root
                                                                              , size_t left
                                                                              , size_t right
                                                                              , size_t middle
                                                                              )
{
    #if defined ( RB_BEAGLE_DEBUG )
        std::stringstream ss;
        ss << "using beagle instance " << std::to_string(this->beagle_instance) << std::endl;
    #endif /* RB_BEAGLE_DEBUG */

    //-- Reset the stored probability.
    this->ln_beagle_probability = 0.0;

    size_t b_model_idx;
    size_t num_taxa  = (this->num_nodes + 2) / 2;

    size_t root_idx  = root + this->num_nodes * this->activeLikelihood[root];
    //this->b_node_indices.push_back(root_idx); //-- TESTING! -- not in original

    size_t mid_idx   = middle + this->num_nodes * this->activeLikelihood[middle];
    size_t left_idx  = left   + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right  + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t mid_partials   = (middle < num_taxa) ? middle : mid_idx;
    size_t left_partials  = (left   < num_taxa) ? left   : left_idx;
    size_t right_partials = (right  < num_taxa) ? right  : right_idx;

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
    int     b_childBufferIndices      = (int) mid_partials;
    int     b_probabilityIndices      = (int) mid_idx;
    int*    b_firstDerivativeIndices  = NULL;
    int*    b_secondDerivativeIndices = NULL;
    int     b_categoryWeightsIndices  = 0; //(int) model;  //0;
    int     b_stateFrequenciesIndices = 0; //(int) model;  //0;
    int     b_cumulativeScaleIndices  = BEAGLE_OP_NONE;
    int     b_count                   = 1;
    double  b_outSumLogLikelihood; //     = NULL; //0;
    double* b_outSumFirstDerivative   = NULL;
    double* b_outSumSecondDerivative  = NULL;

    //-- Return codes for BEAGLE operations.
    int b_code_update_transitions;
    int b_code_update_partials;
    int b_code_calc_edges;

    this->updateBeagleEigensystem();
    this->updateBeagleSiteRates();

    for ( size_t i = 0; i < this->num_site_mixtures; ++i )
    {
        b_model_idx = this->active_eigen_system[i];

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
	if ( b_code_update_partials != 0 )
	{
        throw RbException( "Could not update partials for models '"
	  		             + this->parseBeagleReturnCode(b_code_update_partials));
	}

    //-- Reset the beagle operations queues (used in original)
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
    
    //-- Reset the beagle operations queues (TESTING)
    //this->b_ops.clear();
    //this->b_branch_lengths.clear();
    //this->b_node_indices.clear();

    #if defined ( RB_BEAGLE_DEBUG )
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
    double branch_length = calculateBranchLength(node, node_index);

    size_t node_idx  = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t left_idx  = left       + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx = right      + this->num_nodes * this->activeLikelihood[right];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials  = (node.getChild(0).isTip()) ? left  : left_idx;
    size_t right_partials = (node.getChild(1).isTip()) ? right : right_idx;

    BeagleOperation b_operation =
        { .destinationPartials    = (int) node_idx
        , .destinationScaleWrite  = BEAGLE_OP_NONE
        , .destinationScaleRead   = BEAGLE_OP_NONE
        , .child1Partials         = (int) left_partials
        , .child1TransitionMatrix = (int) left_idx
        , .child2Partials         = (int) right_partials
        , .child2TransitionMatrix = (int) right_idx
        };

    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(node_idx);
    this->b_branch_lengths.push_back(branch_length);
}


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
    double branch_length = calculateBranchLength(node, node_index);

    size_t node_idx   = node_index + this->num_nodes * this->activeLikelihood[node_index];
    size_t left_idx   = left       + this->num_nodes * this->activeLikelihood[left];
    size_t right_idx  = right      + this->num_nodes * this->activeLikelihood[right];
    size_t middle_idx = middle     + this->num_nodes * this->activeLikelihood[middle];

    //-- Tips are actually just stored once, so we dont need offests.
    size_t left_partials   = (node.getChild(0).isTip()) ? left  : left_idx;
    size_t right_partials  = (node.getChild(1).isTip()) ? right : right_idx;
    size_t middle_partials = (node.getChild(2).isTip()) ? middle : middle_idx;

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

    this->b_ops.push_back(b_operation);
    this->b_node_indices.push_back(node_idx);
    this->b_branch_lengths.push_back(branch_length);
}



template<class charType>
void
RevBayesCore::PhyloCTMCSiteHomogeneousBEAGLE<charType>::computeTipLikelihood ( const TopologyNode &node, size_t node_index )
{
    double branch_length = calculateBranchLength(node, node_index);
    this->b_branch_lengths.push_back(branch_length);

    size_t node_idx = node_index + this->num_nodes * this->activeLikelihood[node_index];
    this->b_node_indices.push_back(node_idx);
}


#endif

