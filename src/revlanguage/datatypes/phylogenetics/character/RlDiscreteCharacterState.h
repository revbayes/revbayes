/**
 * @file
 * This file contains the declaration of RlDiscreteCharacterState, which is
 * a RevBayes wrapper around a regular DiscreteCharacter character.
 *
 * @brief Declaration of RlBoolean
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-23 09:27:32 +0200 (Thu, 23 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-09-08, version 1.0
 * @extends RbObject
 *
 * $Id: RlBoolean.h 1766 2012-08-23 07:27:32Z hoehna $
 */

#ifndef RlDiscreteCharacterState_H
#define RlDiscreteCharacterState_H

#include "DiscreteCharacterState.h"
#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class DiscreteCharacterState : public ModelObject<RevBayesCore::DiscreteCharacterState> {
        
    public:
        DiscreteCharacterState(void);                                                                               //!< Default constructor
        DiscreteCharacterState(const RevBayesCore::DiscreteCharacterState &d);                                      //!< Construct from DiscreteCharacter
        
        // Operators
        
        // Basic utility functions
        DiscreteCharacterState*             clone(void) const;                                                      //!< Clone object
        RevObject*                          convertTo(const TypeSpec& type) const;                                  //!< Convert to type
        static const std::string&           getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&              getClassTypeSpec(void);                                                 //!< Get class type spec
        const TypeSpec&                     getTypeSpec(void) const;                                                //!< Get language type of the object
        double                              isConvertibleTo(const TypeSpec& type, bool convert_by_value) const;                 //!< Is convertible to type?
        
        std::string                         getGuiName(void) { return ""; }
        std::string                         getGuiUnicodeSymbol(void) { return ""; }
        std::string                         getGuiInfo(void) { return ""; }
    };
    
}

#endif

