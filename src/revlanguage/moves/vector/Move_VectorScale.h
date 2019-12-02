#ifndef Move_VectorScale_H
#define Move_VectorScale_H

#include "RlMove.h"

namespace RevLanguage {
    
    /**
     * @copybrief RevBayesCore::VectorScaleProposal
     * @see Move_VectorSingleElementScale to scale a single element at a time.
     * @see RevBayesCore::VectorScaleProposal for the internal object.
     */
    class Move_VectorScale : public Move {
        
    public:
        
        Move_VectorScale(void);
        
        // Basic utility functions
        virtual Move_VectorScale*                   clone(void) const;
        void                                        constructInternalObject(void); //!< Build a new internal MH move with a VectorScaleProposal.
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        std::string                                 getMoveName(void) const;
        const MemberRules&                          getParameterRules(void) const;
        virtual const TypeSpec&                     getTypeSpec(void) const;
        virtual void                                printValue(std::ostream& o) const;
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);
        
        RevPtr<const RevVariable>                   x;  //!< The real valued vector on which to operate.
        RevPtr<const RevVariable>                   lambda;  //!< The scaling (and tuning) parameter.
        RevPtr<const RevVariable>                   tune;  //!< Whether autotuning should be used.
        
    };
    
}

#endif
