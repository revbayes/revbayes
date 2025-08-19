#include <cstddef>
#include <sstream>
#include <string>

#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "RbException.h"
#include "TypeSpec.h"
#include "DagNode.h"
#include "Procedure.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Constructor */
MemberProcedure::MemberProcedure(const std::string &n, const TypeSpec retType, ArgumentRules* argRules) : Procedure(),
    argument_rules( argRules ),
    proc_name( n ),
    object( NULL ),
    return_type( retType )
{
    
}



MemberProcedure::MemberProcedure(const MemberProcedure &mp) : Procedure(mp),
    argument_rules( new ArgumentRules(*mp.argument_rules) ),
    proc_name( mp.proc_name ),
    object( mp.object ),
    return_type( mp.return_type )
{
    
}


MemberProcedure::~MemberProcedure( void )
{
    
    delete argument_rules;
}


MemberProcedure& MemberProcedure::operator=(const MemberProcedure &mp)
{
    
    // check for self-assignment
    if ( this != &mp )
    {
        Procedure::operator=( mp );
        
        delete argument_rules;
        
        argument_rules  = new ArgumentRules(*mp.argument_rules);
        proc_name       = mp.proc_name;
        object          = mp.object;
        return_type      = mp.return_type;

    }
    
    
    return *this;
}


/** Clone the object */
MemberProcedure* MemberProcedure::clone(void) const
{
    
    return new MemberProcedure(*this);
}


/** Execute function: call the object's internal implementation through executeOperation */
RevPtr<RevVariable> MemberProcedure::execute( void )
{
    
    bool found = false;
    RevPtr<RevVariable> retValue = object->getRevObject().executeMethod( getFunctionName(), args, found );
    
    if ( found == false )
    {
        throw RbException() << "Couldn't find member procedure called '" << getFunctionName() << "'";
    }
    
    if (object->getRevObject().hasDagNode())
    {
        RevBayesCore::DagNode* the_node = object->getRevObject().getDagNode();
	assert(the_node);

	the_node->touch();
    }
    
    return retValue;
    
}


/** Get class name of object */
const std::string& MemberProcedure::getClassType(void)
{
    
    static std::string revClassType = "MemberProcedure";
    
	return revClassType;
}

/** Get class type spec describing type of object */
const TypeSpec& MemberProcedure::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/** Get type spec */
const TypeSpec& MemberProcedure::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get argument rules */
const ArgumentRules& MemberProcedure::getArgumentRules(void) const
{
    
    return *argument_rules;
}


/**
 * Get the name for this procedure.
 */
std::string MemberProcedure::getFunctionName( void ) const
{
    
    return proc_name;
}


/** Get return type */
const TypeSpec& MemberProcedure::getReturnType(void) const
{
    
    return return_type;
}


/** This is a procedure, so return true for this function call */
bool MemberProcedure::isProcedure( void ) const
{
    return true;
}


void MemberProcedure::setMemberObject( const RevPtr<RevVariable> &obj)
{
    
    // we do not own the object itself because one object can have multiple member functions
    object = obj;
}

