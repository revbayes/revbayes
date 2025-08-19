#ifndef Move_H
#define Move_H

#include <cstddef>
#include <ostream>
#include <vector>
#include <cstdint>

#include "Cloneable.h"

namespace RevBayesCore {
    
    class DagNode;
template <class valueType> class RbOrderedSet;
    
    /**
     * Base class for all move within an MCMC (and all other variants).
     *
     * The base class of all moves only provides the interface for the call to propose a new move.
     * All specifics are implemented in the derived classes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-26, version 1.0
     *
     */
    class Move : public Cloneable {
        
    public:
        virtual                                                ~Move(void);                                                 //!< Destructor
        
        // pure virtual public methods
        virtual void                                            addNode(DagNode* p) = 0;                                    //!< add a node to the proposal
        virtual void                                            autoTune(void) = 0;                                         //!< Automatic tuning of the move.
        virtual Move*                                           clone(void) const = 0;                                      //!< Create a deep copy.
        virtual void                                            decrementTriedCounter(void) = 0;                            //!< Get update weight of InferenceMove
        virtual const RbOrderedSet<DagNode*>&                   getAffectedNodes(void) const = 0;                           //!< Get the nodes vector
        virtual const std::vector<DagNode*>&                    getDagNodes(void) const = 0;                                //!< Get the nodes vector
        virtual const std::string&                              getMoveName(void) const = 0;                                //!< Get the name of the move for summary printing
        virtual double                                          getMoveTuningParameter(void) const = 0;
        virtual size_t                                          getNumberAcceptedCurrentPeriod(void) const = 0;             //!< Get update weight of InferenceMove
        virtual size_t                                          getNumberAcceptedTotal(void) const = 0;                     //!< Get update weight of InferenceMove
        virtual size_t                                          getNumberTriedCurrentPeriod(void) const = 0;                //!< Get the number of tries for this move since the last reset
        virtual size_t                                          getNumberTriedTotal(void) const = 0;                        //!< Get the number of tries for this move since the last reset
        virtual double                                          getUpdateWeight(void) const = 0;                            //!< Get update weight of InferenceMove
        virtual bool                                            isActive(std::uint64_t gen) const = 0;                      //!< Is this move active?
//        virtual void                                            perform(double lHeat, double pHeat) = 0;                  //!< Perform the move.
        virtual void                                            performMcmcStep(double prHeat, double lHeat, double pHeat) = 0;            //!< Perform the move.
        virtual void                                            performHillClimbingStep(double lHeat, double pHeat) = 0;    //!< Perform the move.
        virtual void                                            printSummary(std::ostream &o, bool current_period) const = 0;                    //!< Print the move summary
        virtual void                                            removeNode(DagNode* p) = 0;                                 //!< remove a node from the proposal
        virtual void                                            resetCounters(void) = 0;                                    //!< Reset the counters such as numTried and numAccepted.
        virtual void                                            setMoveTuningParameter(double tp) = 0;
        virtual void                                            setNumberAcceptedCurrentPeriod(size_t na) = 0;
        virtual void                                            setNumberAcceptedTotal(size_t na) = 0;
        virtual void                                            setNumberTriedCurrentPeriod(size_t nt) = 0;
        virtual void                                            setNumberTriedTotal(size_t nt) = 0;
        virtual void                                            swapNode(DagNode *oldN, DagNode *newN) = 0;                 //!< Swap the pointers to the variable on which the move works on.
        
        
    protected:
        Move();                                                                                                             //!< Constructor
        
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Move& x);                                         //!< Overloaded output operator

}


#endif
