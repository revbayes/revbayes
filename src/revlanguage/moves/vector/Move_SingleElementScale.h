#ifndef Move_SingleElementScale_H
#define Move_SingleElementScale_H

#include "RlMove.h"

namespace RevLanguage {

    /**
     * @copybrief RevBayesCore::SingleElementScaleProposal
     * @see Move_VectorSingleElementScale for the same applied to a Stochastic vector.
     * @see RevBayesCore::SingleElementScaleProposal for the internal object.
     **/
    class Move_SingleElementScale : public Move {

    public:

        Move_SingleElementScale(void);

        // Basic utility functions
        virtual Move_SingleElementScale*            clone(void) const;
        void                                        constructInternalObject(void);  //!< Build a new internal MH move with a SingleElementScaleProposal.
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        std::string                                 getMoveName(void) const;
        const MemberRules&                          getParameterRules(void) const;
        virtual const TypeSpec&                     getTypeSpec(void) const;
        virtual void                                printValue(std::ostream& o) const;

    protected:

        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                   x;  //!< The real valued vector to operate on.
        RevPtr<const RevVariable>                   lambda;  //!< The scaling (and tuning) parameter.
        RevPtr<const RevVariable>                   tune;  //!< Whether autotuning should be used.

    };

}

#endif
