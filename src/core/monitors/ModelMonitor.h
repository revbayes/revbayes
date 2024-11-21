#ifndef ModelMonitor_H
#define ModelMonitor_H

#include <iosfwd>
#include <cstdint>
#include <set>

#include "VariableMonitor.h"

namespace RevBayesCore {
class Model;
    
    /**
     * @brief A monitor class that monitors all variables of a model and prints their value into a file.
     *
     * @details
     * The model monitor is a convenience monitor that simply monitors all variables of a model
     * instead of a pre-selected set. Thus, one only needs to specify the model and this monitor
     * extracts all variables that can be monitored.
     * The values will be printed into a file.
     *
     */
    class ModelMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
        ModelMonitor(std::uint64_t g, const path &fname, const SampleFormat &f, std::set<std::string> exclude_list);                        //!< Constructor
        virtual ~ModelMonitor(void);
        
        
        
        // basic methods
        ModelMonitor*                       clone(void) const;                                                  //!< Clone the object
        
        // getters and setters
        void                                setModel(Model* m);
        void                                setStochasticNodesOnly(bool tf);                                    //!< Set if only stochastic nodes should be monitored
    
    private:
        // helper methods
        void                                resetDagNodes(void);                                                //!< Extract the variable to be monitored again.
        
        // members
        bool                                stochastic_nodes_only;                                              //!< Flag if only stochastic nodes should be printed
        std::set<std::string> exclude;  //!< List of variables to exclude from monitoring
    };
    
}

#endif

