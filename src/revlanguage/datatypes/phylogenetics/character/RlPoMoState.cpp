/**
 * @file
 * This file contains the implementation of RlPoMoState, which is
 * a RevBayes wrapper around a regular PoMoState.
 *
 * @brief Implementation of RlPoMoState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-09-04 20:14:58 +0200 (Tue, 04 Sep 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2009-11-20, version 1.0
 * @extends RbObject
 *
 * $Id: RlPoMoState.cpp 1793 2012-09-04 18:14:58Z Boussau $
 */


#include "ConstantNode.h"
#include "RlPoMoState.h"
#include "RbUtil.h"
#include "TypeSpec.h"

#include <sstream>

using namespace RevLanguage;

/** Default constructor */
PoMoState::PoMoState(void) : ModelObject<RevBayesCore::PoMoState>() {
    
}

/** Construct from bool */
PoMoState::PoMoState(const RevBayesCore::PoMoState &d) : ModelObject<RevBayesCore::PoMoState>( new RevBayesCore::PoMoState(d) ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
PoMoState* PoMoState::clone(void) const {
    
	return new PoMoState(*this);
}

/** Get Rev type of object */
const std::string& PoMoState::getClassType(void) {
    
    static std::string rev_type = "PoMo";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& PoMoState::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/** Get type spec */
const TypeSpec& PoMoState::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}

