#ifndef Func_factorial_h
#define Func_factorial_h

#include <string>

#include "RealPos.h"
#include "RlTypedFunction.h"
#include "TypedFunction.h"

namespace RevLanguage {

    class ArgumentRules;
    class TypeSpec;

    class Func_factorial : public TypedFunction<RealPos> {

    public:
        Func_factorial(void);

        Func_factorial* clone(void) const;

        static const std::string& getClassType(void);
        static const TypeSpec& getClassTypeSpec(void);

        std::string getFunctionName(void) const;
        const TypeSpec& getTypeSpec(void) const;

        RevBayesCore::TypedFunction<double>* createFunction(void) const;
        const ArgumentRules& getArgumentRules(void) const;
    };
}

#endif
