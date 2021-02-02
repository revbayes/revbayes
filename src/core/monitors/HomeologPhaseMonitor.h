#ifndef HomeologPhaseMonitor_H
#define HomeologPhaseMonitor_H
#include "AbstractHomologousDiscreteCharacterData.h"
#include "VariableMonitor.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "StochasticNode.h"

#include <fstream>
#include <string>
#include <vector>
#include <typeinfo>

namespace RevBayesCore {
    
    class HomeologPhaseMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
		HomeologPhaseMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* ch, unsigned long g, const std::string &fname, const std::string &del);                                  //!< Constructor
        virtual ~HomeologPhaseMonitor(void);
        
        HomeologPhaseMonitor*          clone(void) const;                                                  //!< Clone the object
        
        // functions you may want to overwrite
        virtual void                                    monitorVariables(unsigned long gen);                                //!< Monitor at generation gen
        virtual void                                    printFileHeader(void);                                              //!< Print header
        virtual void                                    swapNode(DagNode *oldN, DagNode *newN);

    private:
        
        // members
		StochasticNode<AbstractHomologousDiscreteCharacterData>*            ctmc;
    };
    
}


#endif 
