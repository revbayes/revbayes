#ifndef VectorSlideProposal_H
#define VectorSlideProposal_H

#include <cstddef>
#include <ostream>
#include <vector>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Sliding proposal of a all elements of a vector.
     *
     *
     * This proposal randomly slides all elements of a vector using the same sliding factor.
     * A sliding proposal draws a random uniform number u ~ unif (-0.5,0.5)
     * and slides the current vale by a sliding offset
     * delta  = ( lambda * u )
     * where lambda is the tuning parameter of the proposal to influence the size of the proposals.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     *
     */
    class VectorSlideProposal : public Proposal {
        
    public:
        VectorSlideProposal(StochasticNode<RbVector<double> >* n, const std::vector<std::int64_t> &i,double l);                                 //!< Constructor
        
        void                                        cleanProposal(void);                                                                //!< Clean up proposal
        VectorSlideProposal*                        clone(void) const;                                                                  //!< Clone object
        double                                      doProposal(void);                                                                   //!< Perform proposal
        const std::string&                          getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                      getProposalTuningParameter(void) const;
        void                                        printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                        prepareProposal(void);                                                              //!< Prepare the proposal
        void                                        setProposalTuningParameter(double tp);
        void                                        tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                        undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                        swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        StochasticNode<RbVector<double> >*          variable;
        std::vector<std::int64_t>                           indices;
        double                                      lambda;                                                                             //!< The Slide parameter of the Proposal (larger lambda -> larger proposals).
        size_t                                      length;
        double                                      storedSlidingFactor;                                                                        //!< The stored value of the last modified element.
        
        
        
    };
    
}

#endif

