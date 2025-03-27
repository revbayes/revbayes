/**
 * @file
 * This file contains the declaration of RlExpandedStandard, which is
 * a RevBayes wrapper around a regular ExpandedStandard character.
 *
 * @brief Declaration of RlExpandedStandardState
 *
 * (c) Copyright 2014-
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 *
 * $Id: RlExplandedStandard.h $
 */

#ifndef RlExpandedStandardState_H
#define RlExpandedStandardState_H

#include <ostream>
#include <string>
#include <vector>

#include "ExpandedStandardState.h"
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

namespace RevLanguage
{
  class TypeSpec;

  class ExpandedStandardState : public ModelObject<RevBayesCore::ExpandedStandardState>
  {

  public:
    ExpandedStandardState(void);                                                               //!< Default constructor
    ExpandedStandardState(RevBayesCore::TypedDagNode<RevBayesCore::ExpandedStandardState> *v); //!< Constructor from DAG node
    ExpandedStandardState(const RevBayesCore::ExpandedStandardState &d);                       //!< Construct from NaturalNumbersState

    // Operators

    // Basic utility functions
    ExpandedStandardState *clone(void) const;      //!< Clone object
    static const std::string &getClassType(void);  //!< Get Rev type
    static const TypeSpec &getClassTypeSpec(void); //!< Get class type spec
    const TypeSpec &getTypeSpec(void) const;       //!< Get language type of the object

    std::string getGuiName(void) { return ""; }
    std::string getGuiUnicodeSymbol(void) { return ""; }
    std::string getGuiInfo(void) { return ""; }
  };
}

#endif
