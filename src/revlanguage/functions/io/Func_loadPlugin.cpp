#include <sstream>
#include <vector>

#include "Func_loadPlugin.h"

#include "Loader.h"
#include "RlString.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "RlFunction.h"
#include "TypeSpec.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_loadPlugin* Func_loadPlugin::clone( void ) const
{

	return new Func_loadPlugin( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_loadPlugin::execute( void )
{

	// get the information from the arguments for reading the file
	const RlString&    pn  = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
	const RlString&    fn  = static_cast<const RlString&>( args[1].getVariable()->getRevObject() );

	if(pn.getValue() == "TensorPhylo") {

		bool successfulLoad = false;
		if(fn.getValue() == "") {
			successfulLoad = Plugin::loader().loadTensorPhylo();
		} else {
			successfulLoad = Plugin::loader().loadTensorPhylo(fn.getValue());
		}
		if(!successfulLoad) {
			throw RbException("Failed to locate TensorPhylo.");
		}
	} else {
		throw RbException("This plugin is not registered in RevBayes (did you mean TensorPhylo?).");
	}

    return new RevVariable( RlUtils::Void );

}


/** Get argument rules */
const ArgumentRules& Func_loadPlugin::getArgumentRules( void ) const
{

	static ArgumentRules argumentRules = ArgumentRules();
	static bool rules_set = false;

	if (!rules_set)
	{
		argumentRules.push_back( new ArgumentRule( "name", RlString::getClassTypeSpec(), "Name of the plugin (e.g., TensorPhylo).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
		argumentRules.push_back( new ArgumentRule( "path", RlString::getClassTypeSpec(), "Relative or absolute path of the plugin.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
		rules_set = true;

	}

	return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_loadPlugin::getClassType(void)
{

	static std::string rev_type = "Func_loadPlugin";

	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_loadPlugin::getClassTypeSpec(void)
{

	static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_loadPlugin::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "loadPlugin";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_loadPlugin::getTypeSpec( void ) const
{

	static TypeSpec type_spec = getClassTypeSpec();

	return type_spec;
}


/** Get return type */
const TypeSpec& Func_loadPlugin::getReturnType( void ) const
{

    static TypeSpec return_typeSpec = RlUtils::Void;
	return return_typeSpec;
}
