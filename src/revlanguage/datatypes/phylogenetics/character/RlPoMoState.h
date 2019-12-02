/**
 * @file
 * This file contains the declaration of RlPoMoState, which is
 * a RevBayes wrapper around a regular PoMo character.
 *
 * @brief Declaration of RlPoMoState
 *
 * (c) Copyright 2009-
 * @date Last modified: $Date: 2012-08-23 09:27:32 +0200 (Thu, 23 Aug 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-09-08, version 1.0
 * @extends RbObject
 *
 * $Id: RlPoMoState.h 1766 2012-08-23 07:27:32Z Boussau $
 */

#ifndef RlPoMoState_H
#define RlPoMoState_H

#include "PoMoState.h"
#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class PoMoState : public ModelObject<RevBayesCore::PoMoState> {
        
    public:
                                        PoMoState(void);                                //!< Default constructor
                                        PoMoState(const RevBayesCore::PoMoState &d);    //!< Construct from PoMoState
        
        // Operators
        
        // Basic utility functions
        PoMoState*                      clone(void) const;                              //!< Clone object
        static const std::string&       getClassType(void);                             //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                         //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                        //!< Get language type of the object
        
        std::string                     getGuiName(void) { return ""; }
        std::string                     getGuiUnicodeSymbol(void) { return ""; }
        std::string                     getGuiInfo(void) { return ""; }
    };
}

#endif

