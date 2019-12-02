/**
 * @file
 * This file contains the declaration of RlRnaState, which is
 * a RevBayes wrapper around a regular RNA character.
 *
 * @brief Declaration of RlRnaState
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

#ifndef RlRnaState_H
#define RlRnaState_H

#include "ModelObject.h"
#include "RnaState.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class RnaState : public ModelObject<RevBayesCore::RnaState> {
        
    public:
                                        RnaState(void);                                 //!< Default constructor
                                        RnaState(const RevBayesCore::RnaState &d);      //!< Construct from Rna
        
        // Operators
        
        // Basic utility functions
        RnaState*                       clone(void) const;                              //!< Clone object
        static const std::string&       getClassType(void);                             //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                         //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                        //!< Get language type of the object
        
        std::string                     getGuiName(void) { return "RNA State"; }
        std::string                     getGuiUnicodeSymbol(void) { return "S(RNA)"; }
        std::string                     getGuiInfo(void) { return ""; }
    };
}

#endif

