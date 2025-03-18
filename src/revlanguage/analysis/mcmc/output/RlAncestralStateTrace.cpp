/**
 * @file
 * This file contains the implementation of the RevLanguage wrapper of the AncestralStateTrace class.
 *
 * @brief Implementation of RlAncestralStateTrace
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: RlAncestralStateTrace.h $
 */

#include "RlAncestralStateTrace.h"

#include <string>

#include "ArgumentRules.h"
#include "Natural.h"
#include "ArgumentRule.h"
#include "RbException.h"
#include "RevVariable.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

AncestralStateTrace::AncestralStateTrace() : WorkspaceToCoreWrapperObject<RevBayesCore::AncestralStateTrace>()
{

	ArgumentRules* summarizeArgRules = new ArgumentRules();
    summarizeArgRules->push_back( new ArgumentRule("burnin", Natural::getClassTypeSpec(), "The number of iterations to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0LL)) );
    
}



AncestralStateTrace::AncestralStateTrace(const RevBayesCore::AncestralStateTrace& m) : WorkspaceToCoreWrapperObject<RevBayesCore::AncestralStateTrace>( new RevBayesCore::AncestralStateTrace( m ) ) {
    
	ArgumentRules* summarizeArgRules = new ArgumentRules();
    summarizeArgRules->push_back( new ArgumentRule("burnin", Natural::getClassTypeSpec(), "The number of iterations to discard as burnin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0LL)) );

}

AncestralStateTrace::AncestralStateTrace(const AncestralStateTrace& m) : WorkspaceToCoreWrapperObject<RevBayesCore::AncestralStateTrace>( m ) {
    
}



/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */

AncestralStateTrace* AncestralStateTrace::clone(void) const
{
    
	return new AncestralStateTrace(*this);
}



void AncestralStateTrace::constructInternalObject( void )
{
    throw RbException("We do not support a constructor function for AncestralStateTrace.");
}


/* Map calls to member methods */

RevPtr<RevVariable> AncestralStateTrace::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found) {
    
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */

const std::string& AncestralStateTrace::getClassType(void) { 
    
    static std::string rev_type = "AncestralStateTrace";
    
	return rev_type; 
}

/** Get class type spec describing type of object */

const TypeSpec& AncestralStateTrace::getClassTypeSpec(void) { 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::AncestralStateTrace>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}



/** Return member rules (no members) */

const MemberRules& AncestralStateTrace::getParameterRules(void) const {
    
    static MemberRules modelMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set ) {        
        rules_set = true;
    }
    
    return modelMemberRules;
}


/** Get type spec */

const TypeSpec& AncestralStateTrace::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */

void AncestralStateTrace::printValue(std::ostream &o) const {
    
    o << "AncestralStateTrace";
}


/** Set a member variable */

void AncestralStateTrace::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "xxx") {
        
    } 
    else {
        RevObject::setConstParameter(name, var);
    }
}

std::ostream& RevLanguage::operator<<(std::ostream& o, const RevLanguage::AncestralStateTrace& x) {
	//    o << x.getParameterName();
    o << "AncestralStateTrace (";
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
