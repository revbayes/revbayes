#ifndef TajimasDFunction_H
#define TajimasDFunction_H

#include "TypedFunction.h"

namespace RevBayesCore {
class AbstractHomologousDiscreteCharacterData;
class DagNode;
template <class valueType> class TypedDagNode;
    
    /**
     * @brief Computing Tajima's D function.
     *
     * This function computes Tajima's D = (pi - theta) / C.
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
