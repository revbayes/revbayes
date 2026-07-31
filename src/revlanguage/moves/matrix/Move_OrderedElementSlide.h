#ifndef Move_OrderedElementSlide_H
#define Move_OrderedElementSlide_H

#include <ostream>

#include "RlMove.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
class TypeSpec;

    class Move_OrderedElementSlide : public Move {

    public:

        Move_OrderedElementSlide(void);

        virtual Move_OrderedElementSlide*           clone(void) const;
        void                                        constructInternalObject(void);
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        std::string                                 getMoveName(void) const;
        const MemberRules&                          getParameterRules(void) const;
        virtual const TypeSpec&                     getTypeSpec(void) const;
        virtual void                                printValue(std::ostream& o) const;

    protected:

        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

    private:

        RevPtr<const RevVariable>                   x;
    };

}

#endif
