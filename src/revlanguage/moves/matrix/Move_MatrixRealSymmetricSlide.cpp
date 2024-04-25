#include "Move_MatrixRealSymmetricSlide.h"

#include <cstddef>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "MatrixRealSymmetricSingleElementSlidingProposal.h"
#include "RlMatrixRealSymmetric.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "StochasticNode.h"
#include "StringUtilities.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_MatrixRealSymmetricSlide::Move_MatrixRealSymmetricSlide() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_MatrixRealSymmetricSlide* Move_MatrixRealSymmetricSlide::clone(void) const
{
    
	return new Move_MatrixRealSymmetricSlide(*this);
}


void Move_MatrixRealSymmetricSlide::constructInternalObject( void )
{
    // we free the memory first
    delete value;
  
    // now allocate a new wishart simple move
    double l = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* tmp = static_cast<const MatrixRealSymmetric &>( mat->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *matrix = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::MatrixRealSymmetricSingleElementSlidingProposal(matrix,l);
    value = new RevBayesCore::MetropolisHastingsMove(p,w,del,t);
        
}


/** Get class name of object */
const std::string& Move_MatrixRealSymmetricSlide::getClassType(void)
{
    
    static std::string revClassType = "Move_MatrixRealSymmetricSlideMove";
    
	return revClassType;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_MatrixRealSymmetricSlide::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_MatrixRealSymmetricSlide::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SymmetricMatrixElementSlide";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_MatrixRealSymmetricSlide::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , MatrixRealSymmetric::getClassTypeSpec(), "The matrix variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "delta", RealPos::getClassTypeSpec()            , "The sliding window size.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()          , "Should we tune the move during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_MatrixRealSymmetricSlide::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_MatrixRealSymmetricSlide::printValue(std::ostream &o) const
{
    
    o << "Move_MatrixRealSymmetricSlide(";
    if (mat != NULL)
    {
        o << mat->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_MatrixRealSymmetricSlide::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" ) {
        mat = var;
    }
    else if ( name == "delta" ) {
        delta = var;
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
