#ifndef GenericFunction_H
#define GenericFunction_H

#include "TypedDagNode.h"
#include <boost/mp11.hpp>

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

        // WORKAROUND: lambdas with template arguments are only available in C++14
        struct addParameterLambda
        {
            GenericFunction<R,Args...>& F;

            template <typename T>
            void operator()(const T& arg)
            {
                F.addParameter(arg);
            };

            addParameterLambda(GenericFunction<R,Args...>& FF):F(FF) {}
        };

        // WORKAROUND: lambdas with template arguments are only available in C++14
        struct getValueLambda
        {
            template<typename T>
            const T& operator()(const TypedDagNode<T>* node) const
            {
                return node->getValue();
            }
        };

        // WORKAROUND: lambdas with template arguments are only available in C++14
        struct swapParameterLambda
        {
            const DagNode* oldP;
            const DagNode* newP;

            template <typename TDN>
            void operator()(TDN& arg)
            {
                if (arg == oldP)
                    // We have to cast from `const DagNode*` to `const TypedDagNode<T>*`
                    arg = static_cast<TDN>(newP);
            };
        };

    public:

        GenericFunction(R (*f)(Args...), const TypedDagNode<Args>*... args):
            TypedFunction<R>(new R),
            arguments({args...}),
            func(f)
        {
            boost::mp11::tuple_for_each(arguments, addParameterLambda(*this));
        }

        GenericFunction<R, Args...>*  clone(void) const
        {
            return new GenericFunction<R, Args...>(*this);
        }

        void update(void)
        {
            using namespace boost::mp11;

            auto values = tuple_transform(getValueLambda(), arguments);
            *TypedFunction<R>::value = tuple_apply(*func, values);
        }

    protected:
        void swapParameterInternal(const DagNode *oldP, const DagNode *newP)
        {
            boost::mp11::tuple_for_each(arguments, swapParameterLambda{oldP,newP});
        }

    };

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
