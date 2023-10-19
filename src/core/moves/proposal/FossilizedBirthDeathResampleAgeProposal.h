#ifndef FossilizedBirthDeathResampleAgeProposal_H
#define FossilizedBirthDeathResampleAgeProposal_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "Proposal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "StochasticNode.h"

namespace RevBayesCore {
    
    /**
     * The fossilized birth death resample age proposal.
     *
     * This proposal resamples the oldest occurrence age for a single taxon in a fossilized birth death range process.
     * The age is resampled uniformly between the minimum and maximum age.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
	template<class valType>
    class FossilizedBirthDeathResampleAgeProposal : public Proposal {
        
    public:
        FossilizedBirthDeathResampleAgeProposal( StochasticNode<valType> *n);                                               //!<  constructor
        
        // Basic utility functions
        void                                     		  cleanProposal(void);                                        //!< Clean up proposal
        FossilizedBirthDeathResampleAgeProposal<valType>* clone(void) const;                                          //!< Clone object
        double                                   		  doProposal(void);                                           //!< Perform proposal
        const std::string&                       		  getProposalName(void) const;                                //!< Get the name of the proposal for summary printing
        double                                   		  getProposalTuningParameter(void) const;
        void                                     		  prepareProposal(void);                                      //!< Prepare the proposal
        void                                     		  printParameterSummary(std::ostream &o, bool name_only) const;               //!< Print the parameter summary
        void                                     		  setProposalTuningParameter(double tp);
        void                                     		  tune(double r);                                             //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                     		  undoProposal(void);                                         //!< Reject the proposal
        
    protected:
        
        void                                     		  swapNodeInternal(DagNode *oldN, DagNode *newN);             //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        
        // parameters
        StochasticNode<valType>*                 		  variable;                                                   //!< The variable the Proposal is working on
        
        std::vector<double>								  stored_ages;
    };
    
}


template<class valType>
RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::FossilizedBirthDeathResampleAgeProposal( RevBayesCore::StochasticNode<valType> *n ) : RevBayesCore::Proposal(),
    variable( n )
{
    // tell the base class to add the node
    addNode( variable );

    AbstractFossilizedBirthDeathRangeProcess* dist = dynamic_cast<AbstractFossilizedBirthDeathRangeProcess* >(&variable->getDistribution());

    if ( dist == NULL )
    {
    	throw RbException("FossilizedBirthDeathResampleAgeProposal can only be used with Fossilized Birth Death Range Processes");
    }
}


/**
 * The cleanProposal function may be called to clean up memory allocations after AbstractMove
 * decides whether to accept, reject, etc. the proposed value.
 *
 */
template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::cleanProposal( void )
{
    ; // do nothing
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the proposal.
 */
template<class valType>
RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>* RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::clone( void ) const
{

    return new FossilizedBirthDeathResampleAgeProposal<valType>( *this );
}


/**
 * Get Proposals' name of object
 *
 * \return The Proposals' name.
 */
template<class valType>
const std::string& RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::getProposalName( void ) const
{
    static std::string name = "FossilizedBirthDeathResampleAge";

    return name;
}

template<class valType>
double RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::getProposalTuningParameter( void ) const
{
    // this proposal has no tuning parameter
    return RbConstants::Double::nan;
}


/**
 * Perform the proposal.
 *
 * A Uniform-simplex proposal randomly changes some values of a simplex, although the other values
 * change too because of the renormalization.
 * First, some random indices are drawn. Then, the proposal draws a new somplex
 *   u ~ Uniform(val[index] * alpha)
 * where alpha is the tuning parameter.The new value is set to u.
 * The simplex is then renormalized.
 *
 * \return The hastings ratio.
 */
template<class valType>
double RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::doProposal( void )
{

    // Get random number generator
    RandomNumberGenerator* rng     = GLOBAL_RNG;

    AbstractFossilizedBirthDeathRangeProcess* dist = dynamic_cast<AbstractFossilizedBirthDeathRangeProcess* >(&variable->getDistribution());

    // touching handled by the distribution
    //stored_ages = dist->getAges();

    size_t i = rng->uniform01() * dist->getAges().size();

    dist->resampleAge(i);

    variable->addTouchedElementIndex(i);

    return 0.0;

}


/**
 *
 */
template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::prepareProposal( void )
{

}


/**
 * Print the summary of the Proposal.
 *
 * The summary just contains the current value of the tuning parameter.
 * It is printed to the stream that it passed in.
 *
 * \param[in]     o     The stream to which we print the summary.
 */
template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::printParameterSummary(std::ostream &o, bool name_only) const
{

}


/**
 * Reject the Proposal.
 *
 * Since the Proposal stores the previous value and it is the only place
 * where complex undo operations are known/implement, we need to revert
 * the value of the variable/DAG-node to its original value.
 */
template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::undoProposal( void )
{
	AbstractFossilizedBirthDeathRangeProcess* dist = dynamic_cast<AbstractFossilizedBirthDeathRangeProcess* >(&variable->getDistribution());

	// restoration handled by the distribution
	//dist->getAges() = stored_ages;

	variable->clearTouchedElementIndices();
}


/**
 * Swap the current variable for a new one.
 *
 * \param[in]     oldN     The old variable that needs to be replaced.
 * \param[in]     newN     The new RevVariable.
 */
template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::swapNodeInternal(DagNode *oldN, DagNode *newN)
{

    variable = static_cast<StochasticNode<valType>* >(newN) ;

}

template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::setProposalTuningParameter(double tp)
{
    // this proposal has no tuning parameter: nothing to do
}


/**
 * Tune the Proposal to accept the desired acceptance ratio.
 *
 * The acceptance ratio for this Proposal should be around 0.44.
 * If it is too large, then we increase the proposal size,
 * and if it is too small, then we decrease the proposal size.
 */
template<class valType>
void RevBayesCore::FossilizedBirthDeathResampleAgeProposal<valType>::tune( double rate )
{

}

#endif

