#ifndef TajimasDFunction_H
#define TajimasDFunction_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Function to calculate Patterson's D statistic.
     *
     * This function calculates Patterson's D = (ABBA - BABA) / (ABBA + BABA).
     *
     * Will Freyman 8/25/15
     *
     */
    class PattersonsDFunction : public TypedFunction<double> {
        
    public:
        PattersonsDFunction(const std::string &p1, const std::string &p2, const std::string &p3, const std::string &og, const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a);
        virtual                                                         ~PattersonsDFunction(void);                                                     //!< Virtual destructor
        
        // public member functions
        PattersonsDFunction*                                            clone(void) const;                                                              //!< Create an independent clone
        void                                                            update(void);
        
    protected:
        void                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Implementation of swap parameters
        
    private:
        
        // members
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >*  alignment;
        const std::string                                               P1;
        const std::string                                               P2;
        const std::string                                               P3;
        const std::string                                               OG;
        
    };
    
}

#endif