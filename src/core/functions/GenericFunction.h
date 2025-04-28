#ifndef GenericFunction_H
#define GenericFunction_H

#include "TypedDagNode.h"
#include "TypedFunction.h"
#include <boost/mp11.hpp>
#include <boost/callable_traits/return_type.hpp>

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
 * Problem #4: We don't want to specify that the function takes each parameter by value or by reference.
 *             We just want the call f(arg1->getValue(), arg2->getValue() ...) to be well-formed.
 * Solution #4: Take the function object as a template parameter.
 *
 *
 * Problem #5:  Template error messages are verbose.
 * Solution #5: Not sure.
 *              C++20 "concepts" allow constraints on template type parameters,
 *                such as that f(arg1->getValue(), arg2->getValue() ...) makes sense.
 *
 * Problem #6: In theory, we want to allow both function POINTERS (like sqrt) and 
 *               callable OBJECTS (like lambda functions).
 *             In practice, we only want to handle functions pointers right now.
 * Solution #6: Provide two different overloaded defs for InvokeFunction.
 *
 * Problem #7: Sometimes we want a TypedFunction<T1>, even though the callable
 *             produces a T2 that is a derived class of T1.
 * Questions: Is seems like we need 3 versions:
 *            - one does *this->value = compute_value();
 *            - one does delete this->value ; this->value = compute_value();
 *            - one does delete this->value ; this->value = compute_value().clone();
 *            - can we complain if, for example, the return type is not a pointer?
 * Solution #7: 
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


// TODO: Remove this ugly block after we can require boost version 1.74
// Define tuple_transform( ) if we have BOOOST version < 1.74
#if BOOST_VERSION < 107400
#include <boost/mp11/integer_sequence.hpp>

namespace boost
{

    namespace mp11
    {

        namespace detail
        {
            template<class... T>
            constexpr auto tp_forward_r( T&&... t )
            {
                return std::tuple<T&&...>( std::forward<T>( t )... );
            }

            template<class... T>
            constexpr auto tp_forward_v( T&&... t )
            {
                return std::tuple<T...>( std::forward<T>( t )... );
            }

            template<std::size_t J, class... Tp>
            constexpr auto tp_extract( Tp&&... tp )
            {
                return tp_forward_r( std::get<J>( std::forward<Tp>( tp ) )... );
            }

            template<class F, class... Tp, std::size_t... J>
            constexpr auto tuple_transform_impl( integer_sequence<std::size_t, J...>, F const& f, Tp&&... tp )
            {
                return tp_forward_v( tuple_apply( f, tp_extract<J>( std::forward<Tp>(tp)... ) )... );
            }

        } // detail

        template<class F, class... Tp,
                 class Z = mp_list<mp_size_t<std::tuple_size<typename std::remove_reference<Tp>::type>::value>...>,
                 class E = mp_if<mp_apply<mp_same, Z>, mp_front<Z>>,
                 class Seq = make_index_sequence<E::value>>
        constexpr auto tuple_transform( F const& f, Tp&&... tp )
        {
            return detail::tuple_transform_impl( Seq(), f, std::forward<Tp>(tp)... );
        }

    } // mp11

} // boost

#endif // BOOST_VERSION < 107400

namespace RevBayesCore
{
    class DagNode;
    template <class valueType> class TypedDagNode;


    // This version handles function POINTERS ( C++ functions are passed as pointers )
    template <class F, class ArgTuple>
    auto InvokeFunction(F* func, const ArgTuple& args_tuple)
    {
        // 1. Get the argument values -- BY REFERENCE.
        // for each i: get<i>(values) = get<i>(arguments)->getValue()
        auto values = boost::mp11::tuple_transform([](auto& node) ->decltype(auto) {return node->getValue();}, args_tuple);

        // 2. Then call the function and store the result.
        // *value = func( get<0>(values), get<1>(values), get<2>(values), ...)
        return boost::mp11::tuple_apply(*func, values);
    }

    // This version handles callable objects ( C++ objects with operator()(....) )
    template <class F, class ArgTuple>
    auto InvokeFunction(const F& func, const ArgTuple& args_tuple)
    {
        // 1. Get the argument values -- BY REFERENCE.
        // for each i: get<i>(values) = get<i>(arguments)->getValue()
        auto values = boost::mp11::tuple_transform([](auto& node) -> decltype(auto) {return node->getValue();}, args_tuple);

        // 2. Then call the function and store the result.
        // *value = func( get<0>(values), get<1>(values), get<2>(values), ...)
        return boost::mp11::tuple_apply(func, values);
    }

    // The callable must return a pointer!
    template <class ValueType, class F, class ...Args>
    class GenericFunction : public TypedFunction< ValueType >
    {
        // The arguments
        std::tuple<const TypedDagNode<Args>*...>  arguments;

        // The callable that we are wrapping.
        F func;

        typedef boost::callable_traits::return_type_t<F> ReturnType;

    public:

        static_assert( std::is_pointer<ReturnType>::value , "Function must return a pointer");

        GenericFunction(F f, const TypedDagNode<Args>*... args):
            TypedFunction<ValueType>(nullptr),
            arguments({args...}),
            func(f)
        {
            // for each i: addParameter(get<i>(arguments))
            boost::mp11::tuple_for_each(arguments, [this](const auto& arg) {
                this->addParameter(arg);
            });
        }

        GenericFunction<ValueType, F, Args...>*  clone(void) const
        {
            return new GenericFunction<ValueType, F, Args...>(*this);
        }

        ValueType* compute_value()
        {
            return InvokeFunction(func, arguments);
        }

        void update(void)
        {
            auto new_value = compute_value();

            // If compute_value() throws an exception, then we won't get here.
            delete this->value;
            this->value = new_value;
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

    // Passing overloaded functions like sqrt means that the compiler needs to know
    // whether to pass
    //    double sqrt(double)
    //    float sqrt(float)
    // generic_function<double,double>(sqrt, arg).

    // So, we will need to either make a NON-overloaded function, or add a type-cast.

    template <typename ValueType, class F, class ...Args>
    GenericFunction<ValueType, F, Args...> generic_function(F&& f, const TypedDagNode<Args>*... args)
    {
        return GenericFunction<ValueType, F, Args...>(f, args...);
    }

    template <typename ValueType, class F, class ...Args>
    GenericFunction<ValueType, F, Args...>* generic_function_ptr(F&& f, const TypedDagNode<Args>*... args)
    {
        return new GenericFunction<ValueType, F,Args...>(f, args...);
    }
}


#endif
