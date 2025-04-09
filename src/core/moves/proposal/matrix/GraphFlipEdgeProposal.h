#ifndef GraphFlipEdgeProposal_H
#define GraphFlipEdgeProposal_H

#include <set>
#include <cstddef>
#include <iosfwd>

#include "Proposal.h"
#include "RbVector.h"
#include "RbVectorImpl.h"

namespace RevBayesCore {
class DagNode;
class MatrixReal;
template <class variableType> class StochasticNode;
    
    /**
     * The sliding operator.
     *
     * This move randomly picks an element of a matrix of positive real numbers.
     * That means, that we randomly pick the i-th row and j-th column with equal probability.
     * Then, we propose a sliding distance and slide the value.
     * The actual sliding distance is computed by delta = lambda * ( u - 0.5 )
     * where u ~ Uniform(0,1).
     * The proposal is thus m[i][j] += lambda * ( u - 0.5 )
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Nicolas Lartillot)
     * @since 2009-09-08, version 1.0
     *
     */
    class GraphFlipEdgeProposal : public Proposal {
        
    public:
        GraphFlipEdgeProposal( StochasticNode<MatrixReal> *n, const RbVector<RbVector<std::int64_t> >& ind, double l, bool s = false);                                                                      //!<  constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                                                //!< Clean up proposal
        GraphFlipEdgeProposal* clone(void) const;                                                                  //!< Clone object
        double                                  doProposal(void);                                                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;                                       //!< Print the parameter summary
        void                                    prepareProposal(void);                                                              //!< Prepare the proposal
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                                     //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                                                 //!< Reject the proposal
        
    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the DAG nodes on which the Proposal is working on
        
    private:
        // parameters
        
        StochasticNode<RbVector<RbVector<double> > >* array;
        StochasticNode<MatrixReal>*                   matrix;
        
        double                                  sampling_probability;
        //!< The two indices of the last modified element.
        RbVector<RbVector<std::int64_t> >               edges;
        size_t                                  edge_list_length;
        double                                  storedValue;                                                                       //!< The value we propose.
        bool                                    symmetric;
        
        std::set<size_t>                        touched_edge_elements;
    };
    
}

#endif

