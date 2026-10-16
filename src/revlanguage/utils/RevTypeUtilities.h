#ifndef RevTypeUtilities_H
#define RevTypeUtilities_H

#include "RevVariable.h"
#include "TypeUtilities.h"

namespace RevLanguage {

// Access a non-null variable's object as the assumed Rev type T, without converting it.
// The const reference borrows the object; it does not extend its lifetime.
template <class T>
    requires std::is_base_of_v<RevObject, T> && std::is_same_v<T, std::remove_cvref_t<T>>
const T& objectForRevType(const RevVariable* variable)
{
    return assume_cast<const T&>(variable->getRevObject());
}

// Read a non-null variable through its assumed Rev type, preserving the const accessor's
// value or reference return type. A returned reference does not extend the value's lifetime.
template <class T>
    requires requires (const RevVariable* variable) { objectForRevType<T>(variable).getValue(); }
decltype(auto) valueForRevType(const RevVariable* variable)
{
    return objectForRevType<T>(variable).getValue();
}

// Retrieve a non-null variable's DAG node through its assumed Rev type. Preserve the
// accessor's pointer type and mutability; ownership remains with the existing object.
template <class T>
    requires requires (const RevVariable* variable) { objectForRevType<T>(variable).getDagNode(); }
decltype(auto) dagNodeForRevType(const RevVariable* variable)
{
    return objectForRevType<T>(variable).getDagNode();
}

}

#endif
