/**
 * @file
 * This file contains the implementation of RlBoolean, which is
 * a RevBayes wrapper around a regular bool.
 *
 * @brief Implementation of RlBoolean
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-09-04 20:14:58 +0200 (Tue, 04 Sep 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 * @extends RbObject
 *
 * $Id: RlBoolean.cpp 1793 2012-09-04 18:14:58Z hoehna $
 */


#include <sstream>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "RlStandardState.h"
#include "TypeSpec.h"
#include "CharacterState.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "StandardState.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
StandardState::StandardState(void) : ModelObject<RevBayesCore::StandardState>()
{
    
}

/** Construct from bool */
StandardState::StandardState(const RevBayesCore::StandardState &d) : ModelObject<RevBayesCore::StandardState>( new RevBayesCore::StandardState(d) ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
StandardState* StandardState::clone(void) const
{
    
	return new StandardState(*this);
}


/** Get Rev type of object */
const std::string& StandardState::getClassType(void)
{
    
    static std::string rev_type = "Standard";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& StandardState::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& StandardState::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

