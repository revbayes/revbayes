//
//  Func_vectorFlatten.h
//  revbayes-proj
//
//  Created by Michael Landis on 4/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__Func_vectorFlatten__
#define __revbayes_proj__Func_vectorFlatten__

#include "ModelVector.h"
#include "RlTypedFunction.h"

#include <string>

#include "ArgumentRule.h"
#include "GenericFunction.h"
#include "ModelVector.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlTypedFunction.h"


using namespace RevLanguage;

namespace RevLanguage {

namespace Core = RevBayesCore;

template <typename T>
class Func_vectorFlatten :  public TypedFunction<ModelVector<T> > {
        
public:
    Func_vectorFlatten();
        
    // Basic utility functions
    Func_vectorFlatten<T>*                                                     clone(void) const;                                          //!< Clone the object
    static const std::string&                                               getClassType(void);                                         //!< Get Rev type (static)
    static const TypeSpec&                                                  getClassTypeSpec(void);                                     //!< Get Rev type spec (static)
    std::string                                                             getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
    const TypeSpec&                                                         getTypeSpec(void) const;                                    //!< Get Rev type spec (from instance)
        
    // Regular functions
    Core::TypedFunction< Core::RbVector< typename T::valueType> >*         createFunction(void) const;                                                     //!< Create internal function object
    const ArgumentRules&                                                    getArgumentRules(void) const;                               //!< Get argument rules
        
};
    

template <typename T>
Core::RbVector<T>* flatten(const Core::RbVector< Core::RbVector<T>>& vs)
{
    auto f = new Core::RbVector<T>;
    for(auto& v: vs)
	for(auto& x: v)
	    f->push_back(x);
    return f;
}


/** Default constructor */
template <typename T>
Func_vectorFlatten<T>::Func_vectorFlatten( void ) :
TypedFunction<ModelVector<T> >()
{
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
template <typename T>
Func_vectorFlatten<T>* Func_vectorFlatten<T>::clone( void ) const
{
    return new Func_vectorFlatten( *this );
}


/** Execute function: Compute simplex from vector of RealPos values */
template <typename T>
Core::TypedFunction< Core::RbVector<typename T::valueType> >* Func_vectorFlatten<T>::createFunction( void ) const
{
    typedef typename T::valueType CT;
    const Core::TypedDagNode< Core::RbVector<Core::RbVector<CT> > >* vec;
    vec = static_cast< const ModelVector<ModelVector<T> >& >( this->args[0].getVariable()->getRevObject() ).getDagNode();

    return Core::generic_function_ptr<Core::RbVector<CT>>(flatten<CT>, vec);
}


/** Get argument rules */
template <typename T>
const ArgumentRules& Func_vectorFlatten<T>::getArgumentRules( void ) const
{
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "x", ModelVector<ModelVector<T> >::getClassTypeSpec(), "A vector of a vector.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }

    return argumentRules;
}


/** Get Rev type of object (static version) */
template <typename T>
const std::string& Func_vectorFlatten<T>::getClassType( void )
{
    static std::string rev_type = "Func_vectorFlatten";

    return rev_type;
}


/** Get Rev type spec of object (static version) */
template <typename T>
const TypeSpec& Func_vectorFlatten<T>::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), &Function::getClassTypeSpec() );

    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
template <typename T>
std::string Func_vectorFlatten<T>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "vectorFlatten";

    return f_name;
}


/** Get Rev type spec from an instance of the object */
template <typename T>
const TypeSpec& Func_vectorFlatten<T>::getTypeSpec( void ) const
{
    return getClassTypeSpec();
}

}


#endif /* defined(__revbayes_proj__Func_vectorFlatten__) */
