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
#include "RlDnaState.h"
#include "TypeSpec.h"
#include "CharacterState.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DnaState.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "StringUtilities.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** Default constructor */
DnaState::DnaState(void) : ModelObject<RevBayesCore::DnaState>() {
    
}

/** Construct from bool */
DnaState::DnaState(const RevBayesCore::DnaState &d) : ModelObject<RevBayesCore::DnaState>( new RevBayesCore::DnaState(d) ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
DnaState* DnaState::clone(void) const {
    
	return new DnaState(*this);
}


/** Get Rev type of object */
const std::string& DnaState::getClassType(void) { 
    
    static std::string rev_type = "DNA";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& DnaState::getClassTypeSpec(void) { 
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/** Get type spec */
const TypeSpec& DnaState::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

