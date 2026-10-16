#ifndef RevVariableUtilities_H
#define RevVariableUtilities_H

#include "RevVariable.h"
#include "TypeUtilities.h"

namespace RevLanguage {

    // These helpers require a non-null variable containing the assumed Rev wrapper type T.
    // Returned references borrow from its object and do not extend the object's lifetime.

    // Access the object as const T, checking the assumed type in debug builds.
    template <class T>
    const T& objectFor(const RevVariable* variable)
    {
        return RevBayesCore::assume_cast<const T&>(variable->getRevObject());
    }

    // Read through the const object accessor, preserving its value or reference return type.
    template <class T>
    decltype(auto) valueFor(const RevVariable* variable)
    {
        return objectFor<T>(variable).getValue();
    }

    // Retrieve the object's DAG node, retaining the accessor's existing mutable pointer type.
    template <class T>
    auto dagNodeFor(const RevVariable* variable)
    {
        return objectFor<T>(variable).getDagNode();
    }

}

#endif
