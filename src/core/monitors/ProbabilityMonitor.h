#ifndef ProbabilityMonitor_H
#define ProbabilityMonitor_H

#include <iosfwd>

#include "VariableMonitor.h"

namespace RevBayesCore {
class Model;
    
    /**
     * @brief A monitor class that monitors all stochastic variables of a model and prints their probability into a file.
     *
     * @file
     * The probability monitor is a convenience monitor that simply monitors all stochastic variables of a model.
     * The current probability of the values will be printed into a file.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-06-21, version 1.0
     *
     */
    class ProbabilityMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
        ProbabilityMonitor(std::uint64_t g, const std::string &fname, const std::string &del);                                  //!< Constructor
        virtual ~ProbabilityMonitor(void);
        
        
        
        // basic methods
        ProbabilityMonitor*                 clone(void) const;                                                  //!< Clone the object
        
        // public (overloaded) methods
        void                                monitorVariables(std::uint64_t gen);                                //!< Monitor at generation gen

        
        // getters and setters
        void                                setModel(Model* m);
        
    private:
        // helper methods
        void                                resetDagNodes(void);                                                //!< Extract the variable to be monitored again.
        
    };
    
}

#endif

