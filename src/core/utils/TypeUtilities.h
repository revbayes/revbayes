#ifndef TypeUtilities_H
#define TypeUtilities_H

#include "RbException.h"

#include <memory>
#include <type_traits>
#include <typeinfo>

// Use cast_to_derived when you have a pointer or reference to a base class and already know
// the object is an instance of the derived class you need. For example, if a RevObject&
// refers to a RealPos, cast_to_derived<RealPos&>(object) refers to that same object as a RealPos.
//
// Like boost::polymorphic_downcast, cast_to_derived checks the cast in debug builds.
// Boost uses BOOST_ASSERT on failure. We throw RbException and report the object's
// dynamic type along with the requested type.
// Release builds use static_cast without a runtime check.

namespace RevBayesCore {

// Debug builds use dynamic_cast and release builds use static_cast, so both casts must compile.
template <class Base, class DerivedPtr>
concept CanCastBasePointerTo = std::is_pointer_v<DerivedPtr> && std::is_class_v<std::remove_pointer_t<DerivedPtr>> &&
    requires (Base* value) {
        static_cast<DerivedPtr>(value);
        dynamic_cast<DerivedPtr>(value);
    };

}

// A null pointer stays null.
template <class DerivedPtr, class Base>
    requires RevBayesCore::CanCastBasePointerTo<Base, DerivedPtr>
DerivedPtr cast_to_derived(Base* value)
{
#ifndef NDEBUG
    // With multiple inheritance, dynamic_cast can reach the target from a different base
    // subobject. Comparing the pointers checks that the target leads back to the original base.
    if (dynamic_cast<DerivedPtr>(value) != value)
        throw RbException(RbException::BUG) << "cast_to_derived: cannot cast from " << typeid(Base).name()
                                            << " to " << typeid(std::remove_pointer_t<DerivedPtr>).name()
                                            << " for object of dynamic type " << typeid(*value).name();
#endif
    return static_cast<DerivedPtr>(value);
}

// Require the corresponding pointer cast because a reference cast alone can create a temporary
// and leave the caller with a dangling reference.
template <class DerivedRef, class Base>
    requires std::is_lvalue_reference_v<DerivedRef> &&
             RevBayesCore::CanCastBasePointerTo<Base, std::add_pointer_t<DerivedRef>>
DerivedRef cast_to_derived(Base& value)
{
    return *cast_to_derived<std::add_pointer_t<DerivedRef>>(std::addressof(value));
}

#endif
