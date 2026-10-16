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

}

#endif
