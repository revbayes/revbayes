#include "TypeUtilities.h"
#include "RevTypeUtilities.h"
#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "ContinuousStochasticNode.h"
#include "MetropolisHastingsMove.h"
#include "Move_LevyJumpSum.h"
#include "LevyJumpSumProposal.h"
#include "Probability.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }



using namespace RevLanguage;

Move_LevyJumpSum::Move_LevyJumpSum() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_LevyJumpSum* Move_LevyJumpSum::clone(void) const
{
    
	return new Move_LevyJumpSum(*this);
}


void Move_LevyJumpSum::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double sf = valueForRevType<RealPos>(slideFactor);
    double w = valueForRevType<RealPos>(weight);
    
    RevBayesCore::TypedDagNode<double>* tmpsv1 = dagNodeForRevType<Real>(slideUp);
    RevBayesCore::ContinuousStochasticNode *sv1 = assume_cast<RevBayesCore::ContinuousStochasticNode *>( tmpsv1 );
    
    RevBayesCore::TypedDagNode<double>* tmpsv2 = dagNodeForRevType<Real>(slideDown);
    RevBayesCore::ContinuousStochasticNode *sv2 = assume_cast<RevBayesCore::ContinuousStochasticNode *>( tmpsv2 );
    
    bool tv = valueForRevType<RlBoolean>(tune);
    double tt = valueForRevType<Probability>(tuneTarget);
    
    RevBayesCore::Proposal *p = new RevBayesCore::LevyJumpSumProposal(sv1, sv2, sf);
    p->setTargetAcceptanceRate(tt);
    
    value = new RevBayesCore::MetropolisHastingsMove(p, w, tv);
    
}


/** Get Rev type of object */
const std::string& Move_LevyJumpSum::getClassType(void)
{
    
    static std::string rev_type = "Move_LevyJumpSum";
    
	return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_LevyJumpSum::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_LevyJumpSum::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "LevyJumpSum";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_LevyJumpSum::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "value_1", Real::getClassTypeSpec()     , "", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "value_2", Real::getClassTypeSpec()     , "", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "slide"  , RealPos::getClassTypeSpec()  , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
        memberRules.push_back( new ArgumentRule( "tune"   , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Move_LevyJumpSum::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_LevyJumpSum::printValue(std::ostream &o) const
{
    
    o << "Move_LevyJumpSum(";
    if (slideUp != NULL) {
        o << slideUp->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_LevyJumpSum::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "value_1" ) {
        slideUp = var;
    }
    else if ( name == "value_2" ) {
        slideDown = var;
    }
    else if ( name == "slide" ) {
        slideFactor = var;
    }
    else if ( name == "weight" ) {
        weight = var;
    }
    else if ( name == "tune" ) {
        tune = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}
