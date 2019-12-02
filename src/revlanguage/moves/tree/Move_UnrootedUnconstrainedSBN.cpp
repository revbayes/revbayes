#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "UnrootedUnconstrainedSBNProposal.h"
#include "RevObject.h"
#include "RbException.h"
#include "RealPos.h"
#include "MetropolisHastingsMove.h"
#include "Move_UnrootedUnconstrainedSBN.h"
#include "RlBranchLengthTree.h"
#include "RlSBNParameters.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Move_UnrootedUnconstrainedSBN::Move_UnrootedUnconstrainedSBN() : Move()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_UnrootedUnconstrainedSBN* Move_UnrootedUnconstrainedSBN::clone(void) const
{

	return new Move_UnrootedUnconstrainedSBN(*this);
}


void Move_UnrootedUnconstrainedSBN::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const BranchLengthTree &>( tree->getRevObject() ).getDagNode();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    const RevBayesCore::SBNParameters &sbn = static_cast<const SBNParameters &>(parameters->getRevObject()).getValue();

    RevBayesCore::Proposal *s = new RevBayesCore::UnrootedUnconstrainedSBNProposal(t,sbn);
    value = new RevBayesCore::MetropolisHastingsMove(s,w);

}


/** Get class name of object */
const std::string& Move_UnrootedUnconstrainedSBN::getClassName(void)
{

    static std::string rbClassName = "Move_UnrootedUnconstrainedSBN";

	return rbClassName;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_UnrootedUnconstrainedSBN::getClassTypeSpec(void)
{

    static TypeSpec rbClass = TypeSpec( getClassName(), new TypeSpec( Move::getClassTypeSpec() ) );

	return rbClass;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_UnrootedUnconstrainedSBN::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "UnrootedUnconstrainedSBN";

    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_UnrootedUnconstrainedSBN::getParameterRules(void) const
{

    static MemberRules USBNMemberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {

        USBNMemberRules.push_back( new ArgumentRule( "tree",          BranchLengthTree::getClassTypeSpec(),     "The tree variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        USBNMemberRules.push_back( new ArgumentRule( "parameters",    SBNParameters::getClassTypeSpec() ,       "The parameters of the SBN",                ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        USBNMemberRules.insert( USBNMemberRules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return USBNMemberRules;
}

/** Get type spec */
const TypeSpec& Move_UnrootedUnconstrainedSBN::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}



/** Get type spec */
void Move_UnrootedUnconstrainedSBN::printValue(std::ostream &o) const
{

    o << "UnrootedUnconstrainedSBN(";
    if (tree != NULL)
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void Move_UnrootedUnconstrainedSBN::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "parameters" )
    {
        parameters = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
