#ifndef TajimasDFunction_H
#define TajimasDFunction_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Computing Tajima's D function.
     *
     * This function computes Tajima's D = (pi - theta) / C.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-04-30
     *
     */
    class TajimasDFunction : public TypedFunction<double> {
        
    public:
        TajimasDFunction(const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a, bool e);
        virtual                                                        ~TajimasDFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        TajimasDFunction*                                               clone(void) const;                                                              //!< Create an independent clone
        void                                                            update(void);
        
    protected:
        void                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaDng parameters
        
    private:
        
        // members
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >*  alignment;
        bool                                                            exclude_ambiguous_sites;
    };
    
}

#endif
