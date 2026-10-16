#ifndef RevTypeUtilities_H
#define RevTypeUtilities_H

#include "RevVariable.h"
#include "TypeUtilities.h"

#include <concepts>

namespace RevLanguage {

// Supply the Rev type as the template argument. For example, RealPos is the Rev type
// for an object whose underlying C++ value is a double.
template <class T>
concept RevType = std::derived_from<T, RevObject>;

// Access a non-null variable's object as the assumed Rev type T, without converting it.
// The const reference borrows the object; it does not extend its lifetime.
template <RevType T>
const T& objectForRevType(const RevVariable* variable)
{
    return cast_to_derived<const T&>(variable->getRevObject());
}

// Read a non-null variable through its assumed Rev type, preserving the const accessor's
// value or reference return type. A returned reference does not extend the value's lifetime.
template <RevType T>
decltype(auto) valueForRevType(const RevVariable* variable)
{
    return objectForRevType<T>(variable).getValue();
}

// Retrieve a non-null variable's DAG node through its assumed Rev type. Preserve the
// accessor's pointer type and mutability; ownership remains with the existing object.
template <RevType T>
decltype(auto) dagNodeForRevType(const RevVariable* variable)
{
    return objectForRevType<T>(variable).getDagNode();
}

}

#endif
