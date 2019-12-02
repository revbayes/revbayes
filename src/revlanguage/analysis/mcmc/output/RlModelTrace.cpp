//
//  RlModelTrace.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 3/27/13.
//  Copyright 2013 __MyCompanyName__. All rights reserved.
//

#include "RlModelTrace.h"

#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "Natural.h"
#include "RlUtils.h"

using namespace RevLanguage;

ModelTrace::ModelTrace() : WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace>()
{
    
}


ModelTrace::ModelTrace(const RevBayesCore::ModelTrace &t) : WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace>( new RevBayesCore::ModelTrace( t ) )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
ModelTrace* ModelTrace::clone(void) const
{
    
    return new ModelTrace(*this);
}


void ModelTrace::constructInternalObject( void )
{
    throw RbException("We do not support a constructor function for ModelTrace.");
}


/* Map calls to member methods */
RevPtr<RevVariable> ModelTrace::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& ModelTrace::getClassType(void) {
    
    static std::string rev_type = "ModelTrace";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& ModelTrace::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::ModelTrace>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}



/** Return member rules (no members) */
const MemberRules& ModelTrace::getParameterRules(void) const {
    
    static MemberRules modelMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        rules_set = true;
    }
    
    return modelMemberRules;
}


/** Get type spec */
const TypeSpec& ModelTrace::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void ModelTrace::printValue(std::ostream &o) const {
    
    o << "Filename:                  " << value->getFileName() << std::endl;
    o << "Parameter:                 " << value->getParameterName() << std::endl;
    o << "#Samples:                  " << value->size() << std::endl;
}


/** Set a member variable */
void ModelTrace::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "xxx") {
        
    }
    else {
        RevObject::setConstParameter(name, var);
    }
}


std::ostream& RevLanguage::operator<<(std::ostream& o, const RevLanguage::ModelTrace& x) {
    //    o << x.getParameterName();
    o << "ModelTrace (";
    //    const std::vector<double>& values = x.getValues();
    //    for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it) {
    //        if ( it != values.begin() ) {
    //            o << ", ";
    //        }
    //        o << *it;
    //    }
    o << ")";
    
    return o;
}
