#ifndef Move_VectorSimplexSwap_H
#define Move_VectorSimplexSwap_H

#include "RlMove.h"

namespace RevLanguage {

    /**
     * @copybrief RevBayesCore::VectorSimplexSwapProposal
     * @see RevBayesCore::VectorSimplexSwapProposal for the internal object.
     */
    class Move_VectorSimplexSwap : public Move {

    public:

        Move_VectorSimplexSwap(void);

        // Basic utility functions
        virtual Move_VectorSimplexSwap*             clone(void) const;
        void                                        constructInternalObject(void);  //!< Build a new internal MH move with a VectorSimplexSwapProposal.
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        std::string                                 getMoveName(void) const;
        const MemberRules&                          getParameterRules(void) const;
        virtual const TypeSpec&                     getTypeSpec(void) const;
        virtual void                                printValue(std::ostream& o) const;

    protected:

        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                   x;  //!< The vector on which to operate.
        //RevPtr<const RevVariable>                   lambda;  //!< The scaling (and tuning) parameter
        RevPtr<const RevVariable>                   tune;  //!< Whether autotuning should be used.

    };

}

#endif
