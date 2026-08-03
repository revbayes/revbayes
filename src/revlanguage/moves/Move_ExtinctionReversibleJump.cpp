#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ExtinctionReversibleJumpProposal.h"
#include "MetropolisHastingsMove.h"
#include "Move_ExtinctionReversibleJump.h"
#include "RealPos.h"
#include "RlMatrixReal.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class MatrixReal; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_ExtinctionReversibleJump::Move_ExtinctionReversibleJump() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_ExtinctionReversibleJump* Move_ExtinctionReversibleJump::clone(void) const
{
    
    return new Move_ExtinctionReversibleJump(*this);
}


/**
 * Create a new internal move object.
 *
 * This function simply dynamically allocates a new internal move object that is
 * associated with the variable (DAG-node). The internal move object is created by calling its
 * constructor and passing the move-parameters (the variable and other parameters) as arguments of the
 * constructor. The move constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
void Move_ExtinctionReversibleJump::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    RevBayesCore::Proposal *p = NULL;

    // the range matrix of a dnFBDRP, or the tree of an extended dnFBDSP
    if ( x->getRevObject().isType( TimeTree::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( x->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::Tree> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

        p = new RevBayesCore::ExtinctionReversibleJumpProposal( n );
    }
    else
    {
        RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* tmp = static_cast<const MatrixReal &>( x->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );

        p = new RevBayesCore::ExtinctionReversibleJumpProposal( n );
    }

    // now allocate a new sliding move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_ExtinctionReversibleJump::getClassType(void)
{
    
    static std::string rev_type = "Move_ExtinctionReversibleJump";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_ExtinctionReversibleJump::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_ExtinctionReversibleJump::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ExtinctionRJSwitch";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules are:
 * (1) the variable, which must be a range matrix.
 *
 * \return The member rules.
 */
const MemberRules& Move_ExtinctionReversibleJump::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        std::vector<TypeSpec> types;
        types.push_back( MatrixReal::getClassTypeSpec() );
        types.push_back( TimeTree::getClassTypeSpec() );
        memberRules.push_back( new ArgumentRule( "x", types, "The fossilized birth death range process whose extinction times to jump.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        
        /* Inherit weight (but not tuneTarget!) from Move and put it after the arguments created above */
        const MemberRules& inheritedRules = Move::getParameterRules();
        for (size_t i = 0; i < inheritedRules.size(); ++i)
        {
            if ( inheritedRules[i].getArgumentLabel() == "weight" )
            {
                memberRules.push_back( inheritedRules[i].clone() );
            }
        }
        
        rules_set = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_ExtinctionReversibleJump::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Print the value for the user.
 */
void Move_ExtinctionReversibleJump::printValue(std::ostream &o) const
{
    
    o << "ExtinctionReversibleJump(";
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


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
void Move_ExtinctionReversibleJump::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}




