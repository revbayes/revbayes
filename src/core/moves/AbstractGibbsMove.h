#ifndef AbstractGibbsMove_H
#define AbstractGibbsMove_H

#include <ostream>

#include "AbstractMove.h"

namespace RevBayesCore {

    /**
     * Base class for all Gibbs moves within an MCMC (and all other variants).
     *
     * The base class of all moves only provides the interface for the call to propose a new move.
     * All specifics are implemented in the derived classes.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2015-07-27, version 1.0
     *
     */
    class AbstractGibbsMove : public AbstractMove {

    public:
        virtual                                                ~AbstractGibbsMove(void);                            //!< Destructor

        void                                                    printSummary(std::ostream &o, bool current_period) const;                //!< Print the move summary

        virtual bool                                            heatsAreAllowable(double prHeat, double lHeat, double pHeat); //<! A check we can call for whether it is okay to use this move with given prior/likelihood/posterior heats
        
        // pure virtual public methods
        virtual AbstractGibbsMove*                              clone(void) const = 0;
        virtual const std::string&                              getMoveName(void) const = 0;                        //!< Get the name of the move for summary printing

        virtual double                                          getMoveTuningParameter(void) const;
        virtual void                                            setMoveTuningParameter(double tp);

    protected:
        AbstractGibbsMove(double w, size_t d=0);                                                         //!< Constructor

        void                                                    performMcmcMove(double prHeat, double lHeat, double pHeat);            //!< Perform the move.
        void                                                    tune(void);                                     //!< Specific tuning of the move

        // pure virtual protected methods
        virtual void                                            performGibbsMove(void) = 0;                         //!< Perform the move.

    };

}


#endif
