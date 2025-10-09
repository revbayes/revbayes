#ifndef StochasticVariableMonitor_H
#define StochasticVariableMonitor_H

#include <iosfwd>

#include "VariableMonitor.h"


namespace RevBayesCore {
class Model;
    
    /**
     * @brief A monitor class that monitors all variables of a StochasticVariable and prints their value into a file.
     *
     * @file
     * The StochasticVariable monitor is a convenience monitor that simply monitors all variables of a StochasticVariable
     * instead of a pre-selected set. Thus, one only needs to specify the StochasticVariable and this monitor
     * extracts all variables that can be monitored.
     * The values will be printed into a file.
     *
     * Note that the copy constructor is necessary because streams need to be handled in a particular way.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-06-21, version 1.0
     *
     */
    class StochasticVariableMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
        StochasticVariableMonitor(std::uint64_t g, const path &fname, const std::string &del);           //!< Constructor
        virtual ~StochasticVariableMonitor(void);
        
        
        
        // basic methods
        StochasticVariableMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // getters and setters
        void                                setModel(Model* m);
        
    private:
        // helper methods
        void                                resetDagNodes(void);                                                //!< Extract the variable to be monitored again.
        
    };
    
}

#endif

