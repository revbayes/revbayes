/**
 * @file
 * This file contains the implementation of Func_range, which is
 * the function used to compute range (colon) vectors.
 *
 * @brief Implementation of Func_range
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-05-04 18:03:37 +0200 (Fri, 04 May 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @interface RbFunction
 * @package functions
 * @since Version 1.0, 2009-09-03
 *
 * $Id: Func_source.cpp 1485 2012-05-04 16:03:37Z hoehna $
 */

#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_range.h"
#include "Integer.h"
#include "ModelVector.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Natural.h"
#include "Procedure.h"
#include "RbHelpReference.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
Func_range::Func_range( void ) : Procedure()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_range* Func_range::clone( void ) const
{
    
    return new Func_range( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_range::execute( void )
{
    
    std::int64_t f = static_cast<const Integer &>( args[0].getVariable()->getRevObject() ).getValue();
    std::int64_t l = static_cast<const Integer &>( args[1].getVariable()->getRevObject() ).getValue();
    
    if ( f >= 0 && l >= 0 )
    {

        ModelVector<Natural> *range = new ModelVector<Natural>( RevBayesCore::RbVector<std::int64_t>() );
        if (f < l)
        {
            for ( std::int64_t i = f; i <= l; i++ )
            {
                range->push_back( Natural(i) );
            }
        }
        else
        {
            for ( std::int64_t i = f; i >= l; i-- )
            {
                range->push_back( Natural(i) );
            }
        }
        
        return new RevVariable( range );

    }
    else
    {

        ModelVector<Integer> *range = new ModelVector<Integer>( RevBayesCore::RbVector<std::int64_t>() );
        if (f < l)
        {
            for ( std::int64_t i = f; i <= l; i++ )
            {
                range->push_back( Integer(i) );
            }
        }
        else
        {
            for ( std::int64_t i = f; i >= l; i-- )
            {
                range->push_back( Integer(i) );
            }
        }
    
        return new RevVariable( range );
        
    }
    
}


/** Get argument rules */
const ArgumentRules& Func_range::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if (!rules_set)
    {
        
        argumentRules.push_back( new ArgumentRule( "first", Integer::getClassTypeSpec(), "Lower value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "last" , Integer::getClassTypeSpec(), "Upper value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_range::getClassType(void)
{
    
    static std::string rev_type = "Func_range";
    
	return rev_type; 
}


/** Get class type spec describing type of object */
const TypeSpec& Func_range::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_range::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "range";
    
    return f_name;
}



/** Get type spec */
const TypeSpec& Func_range::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get return type */
const TypeSpec& Func_range::getReturnType( void ) const
{
    
    static TypeSpec return_typeSpec = ModelVector<Integer>::getClassTypeSpec();
    
    return return_typeSpec;
}

