#include <stddef.h>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "Natural.h"
#include "MetropolisHastingsMove.h"
#include "ModelVector.h"
#include "Move_SingleElementScaleSimplexAllScalarVector.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlSimplex.h"
#include "SingleElementScaleSimplexAllScalarVectorProposal.h"
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

Move_SingleElementScaleSimplexAllScalarVector::Move_SingleElementScaleSimplexAllScalarVector() : Move()
{
    // first, the argument rules
    ArgumentRules* addIndexRules = new ArgumentRules();
    
    // next, set the specific arguments
    addIndexRules->push_back( new ArgumentRule( "elementIndex", Natural::getClassTypeSpec(), "The index of the element to scale.", ArgumentRule::BY_VALUE , ArgumentRule::ANY, new Natural( 1 ) ) );
    
    // finally, create the methods
    methods.addFunction( new MemberProcedure( "addIndex", RlUtils::Void, addIndexRules) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_SingleElementScaleSimplexAllScalarVector* Move_SingleElementScaleSimplexAllScalarVector::clone(void) const
{
    
	return new Move_SingleElementScaleSimplexAllScalarVector(*this);
}


void Move_SingleElementScaleSimplexAllScalarVector::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double l = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const RealPos &>( tuneTarget->getRevObject() ).getValue();
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();

    RevBayesCore::TypedDagNode< RevBayesCore::Simplex>* tmp = static_cast<const Simplex &>( x->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode< RevBayesCore::Simplex> *n = static_cast<RevBayesCore::StochasticNode< RevBayesCore::Simplex > *>( tmp );
    
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* tmp_sv = static_cast<const ModelVector<RealPos> &>( scalar_vector->getRevObject() ).getDagNode();
    std::vector<const RevBayesCore::DagNode*> p = tmp_sv->getParents();
    std::vector< RevBayesCore::StochasticNode<double> *> sv;
    for (std::vector<const RevBayesCore::DagNode*>::const_iterator it = p.begin(); it != p.end(); ++it)
    {
        const RevBayesCore::StochasticNode<double> *the_node = dynamic_cast< const RevBayesCore::StochasticNode<double>* >( *it );
        if ( the_node != NULL )
        {
            sv.push_back( const_cast< RevBayesCore::StochasticNode<double>* >( the_node ) );
        }
        else
        {
            throw RbException("Could not create a mvElementScale because the node isn't a vector of stochastic nodes.");
        }
    }
    
    RevBayesCore::Proposal *prop = new RevBayesCore::SingleElementScaleSimplexAllScalarVectorProposal(n, sv, l, r);
    value = new RevBayesCore::MetropolisHastingsMove(prop, w, t);

}


RevPtr<RevVariable> Move_SingleElementScaleSimplexAllScalarVector::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{
    
    if ( name == "addIndex" )
    {
        found = true;
        
        size_t idx = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        
        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::SingleElementScaleSimplexAllScalarVectorProposal &prop = static_cast<RevBayesCore::SingleElementScaleSimplexAllScalarVectorProposal&>( m->getProposal() );
        
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
    
    return Move::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& Move_SingleElementScaleSimplexAllScalarVector::getClassType(void)
{
    
    static std::string rev_type = "Move_SingleElementScaleSimplexAllScalarVector";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_SingleElementScaleSimplexAllScalarVector::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SingleElementScaleSimplexAllScalarVector::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SingleElementScaleSimplexAllScalarVector";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_SingleElementScaleSimplexAllScalarVector::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        move_member_rules.push_back( new ArgumentRule( "simplex"    , Simplex::getClassTypeSpec()  , "The variable this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "scalarVector", ModelVector<RealPos>::getClassTypeSpec(), "The vector contains all the scalar(s) the simplex will be multiplied with.", ArgumentRule::BY_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec()  , "The strength of the proposal.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune" , RlBoolean::getClassTypeSpec(), "Should we tune this move during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_SingleElementScaleSimplexAllScalarVector::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_SingleElementScaleSimplexAllScalarVector::printValue(std::ostream &o) const
{
    
    o << "Move_SingleElementScaleSimplexAllScalarVector(";
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
void Move_SingleElementScaleSimplexAllScalarVector::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "simplex" )
    {
        x = var;
    }
    else if ( name == "scalarVector" )
    {
        scalar_vector = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
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
