/**
 * @file
 * This file contains the declaration of RlDiscretizedContinuousState, which is
 * a RevBayes wrapper around a regular DiscretizedContinuous character.
 *
 * @brief Declaration of RlDiscretizedContinuousState
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: RlDiscretizedContinuousState.h $
 */

#ifndef RlDiscretizedContinuousState_H
#define RlDiscretizedContinuousState_H

#include <ostream>
#include <string>
#include <vector>

#include "DiscretizedContinuousState.h"
#include "ModelObject.h"
#include "TypedDagNode.h"
#include "CharacterState.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "RevPtr.h"
#include "RlConstantNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage {
class TypeSpec;
    
    class DiscretizedContinuousState : public ModelObject<RevBayesCore::DiscretizedContinuousState> {
        
    public:
                                        DiscretizedContinuousState(void);                                                              //!< Default constructor
                                        DiscretizedContinuousState(RevBayesCore::TypedDagNode<RevBayesCore::DiscretizedContinuousState> *v);  //!< Constructor from DAG node
                                        DiscretizedContinuousState(const RevBayesCore::DiscretizedContinuousState &d);                        //!< Construct from DiscretizedContinuousState
        
        // Operators
        
        // Basic utility functions
        DiscretizedContinuousState*     clone(void) const;                                                                      //!< Clone object
        static const std::string&       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&          getClassTypeSpec(void);                                                                 //!< Get class type spec
        const TypeSpec&                 getTypeSpec(void) const;                                                                //!< Get language type of the object
       
        std::string                     getGuiName(void) { return ""; }
        std::string                     getGuiUnicodeSymbol(void) { return ""; }
        std::string                     getGuiInfo(void) { return ""; }
    };
}

#endif
