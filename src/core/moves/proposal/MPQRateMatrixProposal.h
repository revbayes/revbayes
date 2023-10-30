#ifndef MPQRateMatrixProposal_H
#define MPQRateMatrixProposal_H

#include <cstddef>
#include <iosfwd>
#include <gmpxx.h>

#include "Polyhedron.h"
#include "Proposal.h"
#include "RateGenerator.h"
#include "RateMatrix_MPQ.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class valueType> class TypedDagNode;
template <class variableType> class StochasticNode;
    
    /**
     * The time-reversible and non-reversible rate matrix proposal.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (John & Sebastian)
     * @since 2009-09-08, version 1.0
     *
     */
    class MPQRateMatrixProposal : public Proposal {
        
    public:

        MPQRateMatrixProposal( StochasticNode<RateGenerator> *n);                           //!<  constructor
                                                            MPQRateMatrixProposal( const MPQRateMatrixProposal& p);                             //!<  copy constructor

        // Basic utility functions
        void                                                cleanProposal(void);                                                                //!< Clean up proposal
        MPQRateMatrixProposal*                              clone(void) const;                                                                  //!< Clone object
        double                                              doProposal(void);                                                                   //!< Perform proposal
        const std::string&                                  getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                              getProposalTuningParameter(void) const;
        void                                                printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                                prepareProposal(void);                                                              //!< Prepare the proposal
        void                                                setProposalTuningParameter(double tp);
        void                                                tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                                swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        double                                              updateExchangabilityRates(void);
        double                                              updateNonreversibleRates(void);
        double                                              updateStationaryFrequencies(void);
        double                                              updateToNonReversible(void);
        double                                              updateToReversible(void);
        Polyhedron                                          poly;
                
        // parameters
        
        StochasticNode<RateGenerator>*                      variable;
        
        // tuning parameters
        double                                              rev_alpha_pi;                                                                       //!< The Sliding parameter of the move (larger lambda -> larger proposals).
        double                                              rev_alpha_er;                                                                       //!< The Sliding parameter of the move (larger lambda -> larger proposals).
        double                                              non_rev_alpha;                                                                      //!< The Sliding parameter of the move (larger lambda -> larger proposals).
        double                                              rj_alpha;                                                                           //!< The Sliding parameter of the move (larger lambda -> larger proposals).

        //!< The two indices of the last modified element.
        RateMatrix_MPQ                                      stored_Q;

        std::vector<mpq_class>                              W;

    };
    
}

#endif

