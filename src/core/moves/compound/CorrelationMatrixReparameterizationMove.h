#ifndef CorrelationMatrixReparameterizationMove_H
#define CorrelationMatrixReparameterizationMove_H

#include <map>
#include <ostream>
#include <set>
#include <string>

#include "AbstractMove.h"
#include "StochasticNode.h"
#include "PhyloMultivariateBrownianProcessREML.h"
#include "ContinuousCharacterData.h"

namespace RevBayesCore {
    
    class CorrelationMatrixReparameterizationMove : public AbstractMove {
        
    public:
        CorrelationMatrixReparameterizationMove( StochasticNode<MatrixReal> *cm, TypedDagNode<RbVector<double> > *v, StochasticNode<ContinuousCharacterData>* mv, bool t, double w);                                                         //!<  constructor
        virtual                                ~CorrelationMatrixReparameterizationMove(void);                                                             //!< Destructor
        
        // Basic utility functions
        CorrelationMatrixReparameterizationMove*   clone(void) const;                                                                  //!< Clone object
        const std::string&                         getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        double                                     getMoveTuningParameter(void) const;
        size_t                                     getNumberAcceptedCurrentPeriod(void) const;                         //!< Get update weight of InferenceMove
        size_t                                     getNumberAcceptedTotal(void) const;
        void                                       printSummary(std::ostream &o, bool current_period) const;                                                //!< Print the move summary
        void                                       setMoveTuningParameter(double tp);
        void                                       setNumberAcceptedCurrentPeriod(size_t na);
        void                                       setNumberAcceptedTotal(size_t na);
        
    protected:
        void                                       performMcmcMove(double prHeat, double lHeat, double pHeat);                                        //!< Perform move
        void                                       resetMoveCounters(void);                                                            //!< Reset the counters such as numAccepted.
        void                                       swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the pointers to the variable on which the move works on.
        void                                       tune(void);
        
    private:
        void                                       reject(void);
        
        
        // member variables
        StochasticNode<MatrixReal>*                correlation_matrix;
        StochasticNode<ContinuousCharacterData>*   mvbm;
        TypedDagNode<RbVector<double> >*           variance;
        
        // stored objects to undo proposal
        MatrixReal                                 stored_matrix;
        
        size_t                                     num_accepted_current_period;
        size_t                                     num_accepted_total;
        
    };
    
}

#endif
