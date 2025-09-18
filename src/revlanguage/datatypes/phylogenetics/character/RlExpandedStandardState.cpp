#include "RlExpandedStandardState.h"

#include "ConstantNode.h"
#include "TypeSpec.h"
#include "RevObject.h"

using namespace RevLanguage;

/** Default constructor */
ExpandedStandardState::ExpandedStandardState(void) : ModelObject<RevBayesCore::ExpandedStandardState>()
{
}

/** Construct from bool */
ExpandedStandardState::ExpandedStandardState(const RevBayesCore::ExpandedStandardState &d) : ModelObject<RevBayesCore::ExpandedStandardState>(new RevBayesCore::ExpandedStandardState(d))
{
}

ExpandedStandardState::ExpandedStandardState(RevBayesCore::TypedDagNode<RevBayesCore::ExpandedStandardState> *v) : ModelObject<RevBayesCore::ExpandedStandardState>(v)
{
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
ExpandedStandardState *ExpandedStandardState::clone(void) const
{

  return new ExpandedStandardState(*this);
}

/** Get Rev type of object */
const std::string &ExpandedStandardState::getClassType(void)
{

  static std::string rev_type = "ExpandedStandard";

  return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec &ExpandedStandardState::getClassTypeSpec(void)
{

  static TypeSpec rev_type_spec = TypeSpec(getClassType(), new TypeSpec(RevObject::getClassTypeSpec()));

  return rev_type_spec;
}

/** Get type spec */
const TypeSpec &ExpandedStandardState::getTypeSpec(void) const
{

  static TypeSpec type_spec = getClassTypeSpec();

  return type_spec;
}
