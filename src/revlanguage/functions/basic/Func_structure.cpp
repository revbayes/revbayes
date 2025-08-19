#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "DagNode.h"
#include "Func_structure.h"
#include "RlBoolean.h"
#include "RlUserInterface.h"
#include "RlUtils.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "Procedure.h"
#include "RbBoolean.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"

using namespace RevLanguage;

/** Default constructor */
Func_structure::Func_structure( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_structure* Func_structure::clone( void ) const
{
    
    return new Func_structure( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_structure::execute( void )
{
    
    std::ostringstream o;
    
    bool verbose = static_cast< const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();

    // get a reference to the variable from which we want to know the struct
    RevPtr<RevVariable> &the_var = args[0].getVariable();
    
    o << std::endl;
    if ( verbose == true )
    {
        
        // print the name and the memory address
        if ( the_var->getName() == "" )
        {
            o << "_variable     = <unnamed> <" << the_var << ">" << std::endl;
        }
        else
        {
            o << "_variable     = " << the_var->getName() << " <" << the_var << ">" << std::endl;
        }
        
    }
    else
    {
        
        // print only the name
        if ( the_var->getName() == "" )
        {
            o << "_variable     = <unnamed>" << std::endl;
        }
        else
        {
            o << "_variable     = " << the_var->getName() << std::endl;
        }
        
    }
    
    o << "_RevType      = " << the_var->getRevObject().getType() << std::endl;
    o << "_value        = ";
    the_var->printValue(o, true);
    o << std::endl;
    
    if ( the_var->isWorkspaceVariable() && verbose == true )
    {
        o << "_varType      = workspace (control)" << std::endl;
    }
    else if ( the_var->isReferenceVariable() && verbose == true  )
    {
        o << "_varType      = reference" << std::endl;
        o << "_refVar       = " << the_var->getName() << " <" << the_var << ">" << std::endl;
    }
    
    if ( the_var->getRevObject().isModelObject() == true )
    {
        
        const RevObject& obj = the_var->getRevObject();
        RevBayesCore::DagNode* dag_node = obj.getDagNode();
        dag_node->printStructureInfo(o, verbose);
    }
    
    o << std::endl;

    RBOUT( o.str() );

    // no return value
    return NULL;
}


/** Get argument rules */
const ArgumentRules& Func_structure::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "x"      , RevObject::getClassTypeSpec(), "The variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "verbose", RlBoolean::getClassTypeSpec(), "Do you want all the information?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_structure::getClassType(void)
{
    
    static std::string rev_type = "Func_type";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Func_structure::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the alternative Rev names (aliases) for the constructor function.
 *
 * \return Rev aliases of constructor function.
 */
std::vector<std::string> Func_structure::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "str" );
    
    return a_names;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_structure::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "structure";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_structure::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_structure::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = RlUtils::Void;
    
    return return_typeSpec;
}

