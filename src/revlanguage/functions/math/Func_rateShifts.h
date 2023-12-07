#ifndef Func_shiftEvents_H
#define Func_shiftEvents_H

#include "RlOrderedEvents.h"
#include "RlTypedFunction.h"
#include "ShiftEventsFunction.h"

#include <map>
#include <string>


namespace RevLanguage {
    
	template <typename valType>
    class Func_shiftEvents : public TypedFunction< RlOrderedEvents<valType> > {
        
    public:
        Func_shiftEvents();
        
        // Basic utility functions
        Func_shiftEvents*                                                    clone(void) const;                                          //!< Clone the object
        static const std::string&                                            getClassName(void);                                         //!< Get class name
        static const TypeSpec&                                               getClassTypeSpec(void);                                     //!< Get class type spec
        std::string                                                          getFunctionName(void) const;                                //!< Get the primary name of the function in Rev
        const TypeSpec&                                                      getTypeSpec(void) const;                                    //!< Get language type of the object
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::OrderedEvents<typename valType::valueType> >* createFunction(void) const;                                 //!< Create internal function object
        const ArgumentRules&                                                                     getArgumentRules(void) const;                               //!< Get argument rules
        
    };
    
}

template<typename valType>
RevLanguage::Func_shiftEvents<valType>::Func_shiftEvents() : TypedFunction< RlOrderedEvents<valType> >() {

}

template<typename valType>
Func_shiftEvents<valType>* RevLanguage::Func_shiftEvents<valType>::clone( void ) const
{

    return new Func_shiftEvents<valType>( *this );
}

template<typename valType>
RevBayesCore::TypedFunction< RevBayesCore::OrderedEvents<typename valType::valueType> >* RevLanguage::Func_shiftEvents<valType>::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<typename valType::valueType>*                              init   = static_cast<const valType &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
	const RevBayesCore::TypedDagNode<RevBayesCore::OrderedEvents<typename valType::valueType>>* shifts = static_cast<const RlOrderedEvents<valType> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();

    RevBayesCore::ShiftEventsFunction<typename valType::valueType> *func = new RevBayesCore::ShiftEventsFunction<typename valType::valueType>( init, shifts );

	return func;
}

/** Get argument rules */
template<typename valType>
const ArgumentRules& RevLanguage::Func_shiftEvents<valType>::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;

    if ( !rules_set )
    {
    	argumentRules.push_back( new ArgumentRule( "initialValue",  valType::getClassTypeSpec(),                  "The intitial values.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    	argumentRules.push_back( new ArgumentRule( "shiftEvents",   RlOrderedEvents<valType>::getClassTypeSpec(), "The  shift events.",   ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }

    return argumentRules;
}


/** Get class name of object */
template<typename valType>
const std::string& RevLanguage::Func_shiftEvents<valType>::getClassName(void)
{

    static std::string rbClassName = "Func_shiftEvents";

    return rbClassName;
}

/** Get class type spec describing type of object */
template<typename valType>
const RevLanguage::TypeSpec& RevLanguage::Func_shiftEvents<valType>::getClassTypeSpec(void)
{

    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Function::getClassTypeSpec() ) );

    return rbClass;
}

/**
 * Get the primary Rev name for this function.
 */
template<typename valType>
std::string RevLanguage::Func_shiftEvents<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnShiftEvents";

    return f_name;
}

/** Get type spec */
template<typename valType>
const TypeSpec& RevLanguage::Func_shiftEvents<valType>::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}










#endif


