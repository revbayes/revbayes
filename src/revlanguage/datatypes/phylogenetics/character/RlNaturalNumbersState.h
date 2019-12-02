/**
 * @file
 * This file contains the declaration of RlNaturalNumbersState, which is
 * a RevBayes wrapper around a regular NaturalNumbers character.
 *
 * @brief Declaration of RlNaturalNumbersState
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: RlNaturalNumbersState.h $
 */

#ifndef RlNaturalNumbersState_H
#define RlNaturalNumbersState_H

#include "NaturalNumbersState.h"
#include "ModelObject.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class NaturalNumbersState : public ModelObject<RevBayesCore::NaturalNumbersState> {
        
    public:
                                        NaturalNumbersState(void);                                                              //!< Default constructor
                                        NaturalNumbersState(RevBayesCore::TypedDagNode<RevBayesCore::NaturalNumbersState> *v);  //!< Constructor from DAG node
                                        NaturalNumbersState(const RevBayesCore::NaturalNumbersState &d);                        //!< Construct from NaturalNumbersState
        
        // Operators
        
        // Basic utility functions
        NaturalNumbersState*            clone(void) const;                                                                      //!< Clone object
        static const std::string&       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                                                                //!< Get language type of the object
       
        std::string                     getGuiName(void) { return ""; }
        std::string                     getGuiUnicodeSymbol(void) { return ""; }
        std::string                     getGuiInfo(void) { return ""; }
    };
}

#endif
