#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ConjugateInverseWishartMove.h"
#include "Move_ConjugateInverseWishart.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixRealSymmetric.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_ConjugateInverseWishart::Move_ConjugateInverseWishart() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_ConjugateInverseWishart* Move_ConjugateInverseWishart::clone(void) const
{
    
    return new Move_ConjugateInverseWishart(*this);
}


void Move_ConjugateInverseWishart::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixRealSymmetric &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::MatrixReal > *z = static_cast< RevBayesCore::StochasticNode< RevBayesCore::MatrixReal >* >( tmp );
    
    value = new RevBayesCore::ConjugateInverseWishartMove(z, w);
}


/** Get Rev type of object */
const std::string& Move_ConjugateInverseWishart::getClassType(void)
{
    
    static std::string rev_type = "Move_ConjugateInverseWishart";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_ConjugateInverseWishart::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_ConjugateInverseWishart::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ConjugateInverseWishart";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_ConjugateInverseWishart::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , MatrixRealSymmetric::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_ConjugateInverseWishart::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_ConjugateInverseWishart::printValue(std::ostream &o) const
{
    
    o << "ConjugateInverseWishartBrownian(";
    if (x != NULL)
    {
        o << x->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_ConjugateInverseWishart::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}
