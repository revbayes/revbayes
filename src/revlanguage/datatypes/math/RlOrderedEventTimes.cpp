/*
 * RlOrderedEventTimes.cpp
 *
 *  Created on: Apr 9, 2020
 *      Author: mrmay
 */

#include "FunctionTable.h"
#include "Natural.h"
#include "TypeSpec.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "RlOrderedEventTimes.h"
#include "RlMemberFunction.h"

namespace RevLanguage {

/** Default constructor */
RlOrderedEventTimes::RlOrderedEventTimes(void) : ModelObject<RevBayesCore::OrderedEventTimes>()
{
    initMethods();
}

/** Construct from core OrderedEventTimes */
RlOrderedEventTimes::RlOrderedEventTimes(RevBayesCore::OrderedEventTimes *c) : ModelObject<RevBayesCore::OrderedEventTimes>( c )
{
    initMethods();
}

/** Construct from core OrderedEventTimes */
RlOrderedEventTimes::RlOrderedEventTimes(const RevBayesCore::OrderedEventTimes &t) : ModelObject<RevBayesCore::OrderedEventTimes>( new RevBayesCore::OrderedEventTimes( t ) )
{
    initMethods();
}

/** Construct from DAG node */
RlOrderedEventTimes::RlOrderedEventTimes(RevBayesCore::TypedDagNode<RevBayesCore::OrderedEventTimes> *n) : ModelObject<RevBayesCore::OrderedEventTimes>( n )
{
    initMethods();
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
RlOrderedEventTimes* RlOrderedEventTimes::clone(void) const
{
    return new RlOrderedEventTimes(*this);
}

/* Map calls to member methods */
RevLanguage::RevPtr<RevLanguage::RevVariable> RlOrderedEventTimes::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    return ModelObject<RevBayesCore::OrderedEventTimes>::executeMethod( name, args, found );
}

/** Get Rev type of object */
const std::string& RlOrderedEventTimes::getClassType(void)
{
    static std::string rev_type = "OrderedEventTimes";

    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& RlOrderedEventTimes::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );

    return rev_type_spec;
}

/** Get type spec */
const TypeSpec& RlOrderedEventTimes::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}

/**
 * Initialize the member methods.
 */
void RlOrderedEventTimes::initMethods( void )
{

//    ArgumentRules* num_events_arg_rules = new ArgumentRules();
//    this->methods.addFunction( new MemberFunction<RlOrderedEventTimes, Natural >("numEvents", this, num_events_arg_rules ) );

}












} /* namespace RevLanguage */
