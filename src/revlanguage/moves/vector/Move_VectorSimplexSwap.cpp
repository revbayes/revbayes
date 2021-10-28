#include <stddef.h>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_VectorSimplexSwap.h"
#include "VectorSimplexSwapProposal.h"
#include "RealPos.h"
#include "ModelObject.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMove.h"
#include "DeterministicNode.h"
//#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "RlSimplex.h"


//-- TODO : Fix rules to remove lambda


namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevLanguage;

Move_VectorSimplexSwap::Move_VectorSimplexSwap() : Move() {}


Move_VectorSimplexSwap* Move_VectorSimplexSwap::clone(void) const
{    
    return new Move_VectorSimplexSwap(*this);
}


void Move_VectorSimplexSwap::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate the new move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< RevBayesCore::RbVector< RevBayesCore::Simplex > >* tmp = static_cast< const ModelVector< Simplex > &>( x->getRevObject() ).getDagNode();
    //RevBayesCore::StochasticNode< RevBayesCore::RbVector< RevBayesCore::Simplex > >* n = static_cast< RevBayesCore::StochasticNode< RevBayesCore::RbVector< RevBayesCore::Simplex > > * >( tmp );
    RevBayesCore::DeterministicNode< RevBayesCore::RbVector< RevBayesCore::Simplex > >* n = static_cast< RevBayesCore::DeterministicNode< RevBayesCore::RbVector< RevBayesCore::Simplex > > * >( tmp );

    RevBayesCore::Proposal *p = new RevBayesCore::VectorSimplexSwapProposal(n);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, false);    
}


const std::string& Move_VectorSimplexSwap::getClassType(void)
{   
    static std::string rev_type = "Move_VectorSimplexSwap";
    return rev_type;
}


const TypeSpec& Move_VectorSimplexSwap::getClassTypeSpec(void)
{    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );    
    return rev_type_spec;
}


std::string Move_VectorSimplexSwap::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "VectorSimplexSwap";    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the VectorSimplexSwap move are:
 * (1) the variable which must be a positive real vector.
 * (2) the tuning parameter lambda that defines the size of the proposal (positive real)
 * (3) a flag whether auto-tuning should be used.
 *
 * \return The member rules.
 */
const MemberRules& Move_VectorSimplexSwap::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"     , ModelVector<Simplex>::getClassTypeSpec()  , "The vector on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        //move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()               , "The scaling factor (strength) of the proposal.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec()             , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}


const TypeSpec& Move_VectorSimplexSwap::getTypeSpec( void ) const
{    
    static TypeSpec type_spec = getClassTypeSpec();    
    return type_spec;
}


void Move_VectorSimplexSwap::printValue(std::ostream &o) const
{    
    o << "VectorSimplexSwap(";
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


void Move_VectorSimplexSwap::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{   
    if ( name == "x" )
    {
        x = var;
    }
    //else if ( name == "lambda" )
    //{
    //    lambda = var;
    //}
    else if ( name == "tune" )
    {
        tune = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }    
}
