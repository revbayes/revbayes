#ifndef UTIL_VARIANT_H
#define UTIL_VARIANT_H

#include <variant>

template <typename T, typename U>
const T* to(const U& u)
{
    if (std::holds_alternative<T>(u))
        return &std::get<T>(u);
    else
        return nullptr;
}

template <typename T, typename U>
T* to(U& u)
{
    if (std::holds_alternative<T>(u))
        return &std::get<T>(u);
    else
        return nullptr;
}

#endif
