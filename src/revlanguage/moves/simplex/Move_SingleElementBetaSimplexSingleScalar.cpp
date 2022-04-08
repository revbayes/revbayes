#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "Natural.h"
#include "MetropolisHastingsMove.h"
#include "Move_SingleElementBetaSimplexSingleScalar.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlSimplex.h"
#include "SingleElementBetaSimplexSingleScalarProposal.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Simplex; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

Move_SingleElementBetaSimplexSingleScalar::Move_SingleElementBetaSimplexSingleScalar() : Move()
{
    // first, the argument rules
    ArgumentRules* addIndexRules = new ArgumentRules();
    ArgumentRules* addScalarRules = new ArgumentRules();
    
    // next, set the specific arguments
    addIndexRules->push_back( new ArgumentRule( "elementIndex", Natural::getClassTypeSpec(), "The index of the element to scale.", ArgumentRule::BY_VALUE , ArgumentRule::ANY, new Natural( 1 ) ) );
    addScalarRules->push_back( new ArgumentRule( "scalar" , RealPos::getClassTypeSpec(), "The scalar.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addIndex", RlUtils::Void, addIndexRules) );
    methods.addFunction( new MemberProcedure( "addScalar", RlUtils::Void, addScalarRules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_SingleElementBetaSimplexSingleScalar* Move_SingleElementBetaSimplexSingleScalar::clone(void) const
{
    
	return new Move_SingleElementBetaSimplexSingleScalar(*this);
}


void Move_SingleElementBetaSimplexSingleScalar::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double a = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const RealPos &>( tuneTarget->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< RevBayesCore::Simplex>* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::Simplex> *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *prop = new RevBayesCore::SingleElementBetaSimplexSingleScalarProposal(n, a, r);
    value = new RevBayesCore::MetropolisHastingsMove(prop, w, t);

}


RevPtr<RevVariable> Move_SingleElementBetaSimplexSingleScalar::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addIndex" )
    {
        found = true;
        
        size_t idx = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        
        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::SingleElementBetaSimplexSingleScalarProposal &prop = static_cast<RevBayesCore::SingleElementBetaSimplexSingleScalarProposal&>( m->getProposal() );
        
        RevBayesCore::TypedDagNode< RevBayesCore::Simplex>* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode< RevBayesCore::Simplex> *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
        
        if ( idx >= 1 && idx <= n->getValue().size())
        {
            prop.addIndex( idx - 1 );
        }
        else
        {
            std::stringstream ss_err;
            ss_err << "Index out of bounds: The vector of size " << n->getValue().size() << " does not have an element for index " << idx << ".";
            throw RbException(ss_err.str());
        }

        return NULL;
    }
    else if ( name == "addScalar" )
    {
        found = true;
        
        RealPos* tmp_s = dynamic_cast<RealPos *>( &args[0].getVariable()->getRevObject() );
        
        if ( tmp_s != NULL )
        {
            RevBayesCore::StochasticNode<double> *s = dynamic_cast< RevBayesCore::StochasticNode<double> * >( tmp_s->getDagNode() );
            
            RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
            RevBayesCore::SingleElementBetaSimplexSingleScalarProposal &prop = static_cast<RevBayesCore::SingleElementBetaSimplexSingleScalarProposal&>( m->getProposal() );
            
            if ( s != NULL )
            {
                prop.addScalar( s );
            }
            else
            {
                throw RbException("Could not add the node because it isn't a stochastic nodes.");
            }
        }
        
        else
        {
            throw RbException("A problem occured when trying to add " + args[0].getVariable()->getName() + " to the move.");
        }
        
        return NULL;
    }
    
    return Move::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& Move_SingleElementBetaSimplexSingleScalar::getClassType(void)
{
    
    static std::string rev_type = "Move_SingleElementBetaSimplexSingleScalar";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_SingleElementBetaSimplexSingleScalar::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SingleElementBetaSimplexSingleScalar::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SingleElementBetaSimplexSingleScalar";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_SingleElementBetaSimplexSingleScalar::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "x"    , Simplex::getClassTypeSpec()  , "The variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "alpha", RealPos::getClassTypeSpec()  , "The concentration parameter on the current value.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune" , RlBoolean::getClassTypeSpec(), "Should we tune the concentration parameter during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_SingleElementBetaSimplexSingleScalar::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_SingleElementBetaSimplexSingleScalar::printValue(std::ostream &o) const
{
    
    o << "Move_SingleElementBetaSimplexSingleScalar(";
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
void Move_SingleElementBetaSimplexSingleScalar::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "x" )
    {
        x = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
    }
    else if ( name == "tune" )
    {
        tune = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
