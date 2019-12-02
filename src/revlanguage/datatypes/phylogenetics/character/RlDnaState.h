/**
 * @file
 * This file contains the declaration of RlDnaState, which is
 * a RevBayes wrapper around a regular DNA character.
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

#ifndef RlDnaState_H
#define RlDnaState_H

#include "DnaState.h"
#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class DnaState : public ModelObject<RevBayesCore::DnaState> {
        
    public:
        DnaState(void);                                                                                         //!< Default constructor
        DnaState(const RevBayesCore::DnaState &d);                                                              //!< Construct from DNA
        
        // Operators
        
        // Basic utility functions
        DnaState*                       clone(void) const;                                                      //!< Clone object
        static const std::string&       getClassType(void);                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                 //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                                                //!< Get language type of the object
        
        std::string                     getGuiName(void) { return "DNA State"; }
        std::string                     getGuiUnicodeSymbol(void) { return "S(DNA)"; }
        std::string                     getGuiInfo(void) { return ""; }
    };
    
}

#endif

