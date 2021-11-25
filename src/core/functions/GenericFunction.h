#ifndef GenericFunction_H
#define GenericFunction_H

#include "TypedDagNode.h"
#include <boost/mp11.hpp>

/* Overview:
 *
 * We want to make a TypedFunction< > that links arg1, arg2, ... to f(arg1->getValue(), arg2->getValue() ...).
 * The arguments are all TypeDagNode<T> objects for different types T.
 *
 * Problem #1: Different functions take  different numbers of arguments.
 * Solution #1: Use template parameters packs like "...Args".
 *
 * Problem #2:  We store the arguments in a tuple.  But we can't use a for loop over a tuple.
 *              So, how do we do something to each argument, such as get its value?
 *              It would be easy to do this by hand -- how do we do it programmatically?
 * Solution #2: The Boost.mp11 library allows us to do something to each element in a tuple.
 *
 * Problem #3: How do we describe the "thing" that we want to do to each element of the tuple?
 * Solution #3: We use lambda functions with an argument of type "auto".
 *              This makes the argument type into a template parameter.
 *
 */

/*
 * C++11/14 technologies that we use:
 *
 * - std::tuple    See: https://en.cppreference.com/w/cpp/utility/tuple
 *   + extends std::pair to more than two members.
 *   + Example: std::tuple<int, double, int, char> t;
 *
 * - template parameter packs:    See: https://en.cppreference.com/w/cpp/language/parameter_pack
 *   + allows templates to have a variable number of arguments, like variadic functions.
 *   + Example: template <typename ReturnType, typename ...ArgumentTypes> ReturnType function(ArgumentTypes ...args) {  }
 *
 * - lambda functions:     See: https://en.cppreference.com/w/cpp/language/lambda
 *   + allows passing one-off functions to other functions.
 *   + similar to a function pointer, but can refer to ("capture") local variables.
 *   + Example: std::sort (v.begin(), v.end(), [](auto& x, auto& y) {return x.second < y.second} )
 *
 * - initializer lists for aggregate initialization    See: https://en.cppreference.com/w/cpp/language/aggregate_initialization
 *   + Example: std::tuple<double,int,int> t = {1.0, 2, 3}
 *
 */

namespace RevBayesCore
{
    class DagNode;
    template <class valueType> class TypedDagNode;

    template <class R, class ...Args>
    class GenericFunction : public TypedFunction<R>
    {
        // The arguments
        std::tuple<const TypedDagNode<Args>*...>  arguments;

        // The function we are wrapping.
        R (*func)(Args...);

    public:

        GenericFunction(R (*f)(Args...), const TypedDagNode<Args>*... args):
            TypedFunction<R>(new R),
            arguments({args...}),
            func(f)
        {
            // for each i: addParameter(get<i>(arguments))
            boost::mp11::tuple_for_each(arguments, [this](const auto& arg) {
                this->addParameter(arg);
            });
        }

        GenericFunction<R, Args...>*  clone(void) const
        {
            return new GenericFunction<R, Args...>(*this);
        }

        void update(void)
        {
            // 1. Get the argument values
            // for each i: get<i>(values) = get<i>(arguments)->getValue()
            auto values = boost::mp11::tuple_transform([](auto& node) {return node->getValue();}, arguments);

            // 2. Then call the function and store the result.
            // *value = func( get<0>(values), get<1>(values), get<2>(values), ...)
            *this->value = boost::mp11::tuple_apply(*func, values);
        }

    protected:
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP)
        {
            auto F2 = [=](auto& arg) {
                using TDN = typename std::remove_reference<decltype(arg)>::type;
                if (arg == oldP)
                    // We have to cast from `const DagNode*` to `const TypedDagNode<T>*`
                    arg = static_cast<TDN>(newP);
            };
            boost::mp11::tuple_for_each(arguments, F2);
        }

    };

    // These functions are supposed to allow us to write generic_function(sqrt,arg) instead of
    // generic_function<double,double>(sqrt, arg).

    // However, they don't work on OS X.  Why not?
    template <class R, class ...Args>
    GenericFunction<R,Args...> generic_function(R (*f)(Args...), const TypedDagNode<Args>*... args)
    {
        return GenericFunction<R,Args...>(f, args...);
    }

    template <class R, class ...Args>
    GenericFunction<R,Args...>* generic_function_ptr(R (*f)(Args...), const TypedDagNode<Args>*... args)
    {
        return new GenericFunction<R,Args...>(f, args...);
    }
}


#endif
