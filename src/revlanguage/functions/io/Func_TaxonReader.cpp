#include "Func_TaxonReader.h"

#include <sstream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "Delimiter.h"
#include "ModelVector.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "TaxonReader.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_TaxonReader* Func_TaxonReader::clone( void ) const
{
    
    return new Func_TaxonReader( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_TaxonReader::execute( void )
{
    
    // get the information from the arguments for reading the file
    const RlString& fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() );
    std::string del = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::TaxonReader tr = RevBayesCore::TaxonReader(fn.getValue(), del);
    const std::vector<RevBayesCore::Taxon>& taxa = tr.getTaxa();
    
    return new RevVariable( new ModelVector<Taxon>( taxa ) );
}


/** Get argument rules */
const ArgumentRules& Func_TaxonReader::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;
    
    if (!rules_set)
    {
        argumentRules.push_back( new ArgumentRule( "filename", RlString::getClassTypeSpec(), "Relative or absolute file name.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new Delimiter() );
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_TaxonReader::getClassType(void)
{
    
    static std::string rev_type = "Func_TaxonReader";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_TaxonReader::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_TaxonReader::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readTaxonData";
    
    return f_name;
}


/** Get type spec */
const TypeSpec& Func_TaxonReader::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_TaxonReader::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = ModelVector<Taxon>::getClassTypeSpec();
    return return_typeSpec;
}




