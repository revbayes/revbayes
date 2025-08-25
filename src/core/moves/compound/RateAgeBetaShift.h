#ifndef RateAgeBetaShift_H
#define RateAgeBetaShift_H

#include <map>
#include <ostream>
#include <set>
#include <string>

#include "AbstractMove.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {
    
    class RateAgeBetaShift : public AbstractMove {
        
    public:
        RateAgeBetaShift( StochasticNode<Tree> *tr, std::vector<StochasticNode<double> *> n, StochasticNode<RbVector<double> > *sn, double delta, bool t, double weight);                                                         //!<  constructor
        virtual                                ~RateAgeBetaShift(void);                                                             //!< Destructor
        
        // Basic utility functions
        RateAgeBetaShift*                       clone(void) const;                                                                  //!< Clone object
        const std::string&                      getMoveName(void) const;                                                            //!< Get the name of the move for summary printing
        double                                  getMoveTuningParameter(void) const;
        size_t                                  getNumberAcceptedCurrentPeriod(void) const;                         //!< Get update weight of InferenceMove
        size_t                                  getNumberAcceptedTotal(void) const;                                 //!< Get update weight of InferenceMove
        void                                    printSummary(std::ostream &o, bool current_period) const;                                                //!< Print the move summary
        void                                    setMoveTuningParameter(double tp);
        void                                    setNumberAcceptedCurrentPeriod(size_t na);
        void                                    setNumberAcceptedTotal(size_t na);
        
    protected:
        void                                    performMcmcMove(double prHeat, double lHeat, double pHeat);                                        //!< Perform move
        void                                    resetMoveCounters(void);                                                            //!< Reset the counters such as numAccepted.
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                                     //!< Swap the pointers to the variable on which the move works on.
        void                                    tune(void);
        
    private:
        void                                    reject(void);
        
        // member variables
        StochasticNode<Tree>*                   tree;
        std::vector<StochasticNode<double>* >   rates_vec;
        StochasticNode<RbVector<double> >*      rates;
        double                                  delta;
        
        // stored objects to undo proposal
        TopologyNode*                           stored_node = NULL;
        double                                  stored_age;
        std::vector<double>                     stored_rates;

        size_t                                  num_accepted_current_period;
        size_t                                  num_accepted_total;

    };
    
}

#endif
