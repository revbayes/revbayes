#ifndef VectorSingleElementSlideMove_H
#define VectorSingleElementSlideMove_H

#include <ostream>
#include <cstdint>
#include <vector>

#include "Proposal.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class RbVector;
template <class variableType> class StochasticNode;
    
    /**
     * @brief Sliding move of a single element randomly picked from a vector.
     *
     *
     * This move randomly picks an element of a vector of positive real numbers,
     * proposes a sliding offset and then slides the value.
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
    class SynchronizedVectorFixedSingleElementSlideProposal : public Proposal {
        
    public:
        SynchronizedVectorFixedSingleElementSlideProposal(const std::vector<StochasticNode<RbVector<double> >* > &n, double l, const std::vector<std::int64_t> &i);  //!< Constructor
        
        void                                                    cleanProposal(void);                                                                //!< Clean up proposal
        SynchronizedVectorFixedSingleElementSlideProposal*      clone(void) const;                                                                  //!< Clone object
        double                                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                                  getProposalTuningParameter(void) const;
        void                                                    printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                                    setProposalTuningParameter(double tp);
        void                                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        
        void                                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
        
    private:
        // parameters
        
        std::vector< StochasticNode<RbVector<double> >* >       variables;
        
        double                                                  lambda;                                                                             //!< The Slide parameter of the move (larger lambda -> larger proposals).
        std::vector<std::int64_t>                                       indices;                                                                              //!< The index of the last modified element.
        double                                                  stored_delta;                                                                        //!< The stored value of the last modified element.
        
        
        
    };
    
}

#endif

