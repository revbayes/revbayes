#include "Func_factorial.h"

#include <boost/math/special_functions/factorials.hpp>
#include <cstdint>

#include "Natural.h"
#include "GenericFunction.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RevVariable.h"
#include "TypeSpec.h"

using namespace RevLanguage;

double* my_factorial(std::int64_t n)
{
    return new double(
        boost::math::factorial<double>(static_cast<unsigned int>(n))
    );
}

Func_factorial::Func_factorial(void) : TypedFunction<RealPos>()
{
}

Func_factorial* Func_factorial::clone(void) const
{
    return new Func_factorial(*this);
}

RevBayesCore::TypedFunction<double>* Func_factorial::createFunction(void) const
{
    RevBayesCore::TypedDagNode<std::int64_t>* n =
        static_cast<const Natural&>(
            this->args[0].getVariable()->getRevObject()
        ).getDagNode();

    return RevBayesCore::generic_function_ptr<double>(my_factorial, n);
}

const ArgumentRules& Func_factorial::getArgumentRules(void) const
{
    static ArgumentRules argumentRules;
    static bool rules_set = false;

    if (!rules_set)
    {
        argumentRules.push_back(
            new ArgumentRule(
                "n",
                Natural::getClassTypeSpec(),
                "The natural number whose factorial is calculated.",
                ArgumentRule::BY_CONSTANT_REFERENCE,
                ArgumentRule::ANY
            )
        );

        rules_set = true;
    }

    return argumentRules;
}

const std::string& Func_factorial::getClassType(void)
{
    static std::string rev_type = "Func_factorial";
    return rev_type;
}

const TypeSpec& Func_factorial::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec(
        getClassType(),
        new TypeSpec(Function::getClassTypeSpec())
    );

    return rev_type_spec;
}

std::string Func_factorial::getFunctionName(void) const
{
    return "factorial";
}

const TypeSpec& Func_factorial::getTypeSpec(void) const
{
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}
