#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "Move_RateAgeBetaShift.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "VectorFunction.h"
#include "DeterministicNode.h"
#include "ModelObject.h"
#include "Move.h"
#include "RateAgeBetaShift.h"
#include "RbBoolean.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "TypedFunction.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class RbVector; }


using namespace RevLanguage;

Move_RateAgeBetaShift::Move_RateAgeBetaShift() : Move() {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_RateAgeBetaShift* Move_RateAgeBetaShift::clone(void) const
{
    
	return new Move_RateAgeBetaShift(*this);
}


void Move_RateAgeBetaShift::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    double d = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    bool at = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* tmpRates = static_cast<const ModelVector<RealPos> &>( rates->getRevObject() ).getDagNode();
    std::vector< RevBayesCore::StochasticNode<double> *> rates;
    RevBayesCore::StochasticNode< RevBayesCore::RbVector<double> >* snode_rates = dynamic_cast<RevBayesCore::StochasticNode< RevBayesCore::RbVector<double> > *>( tmpRates );
    if ( tmpRates->isStochastic() == false )
    {
        RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> >*dnode = static_cast< RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> > *>( tmpRates );

        RevBayesCore::VectorFunction<double>* func_vec = dynamic_cast<RevBayesCore::VectorFunction<double>*>( &dnode->getFunction() );
        if ( func_vec == NULL )
        {
            throw RbException("Problem in RateAgeBetaShift move. Wrong argument type for the rates vector. We expect a vector of iid elements.");
        }
        const std::vector<const RevBayesCore::TypedDagNode<double>* >& pars = func_vec->getVectorParameters();

        for (std::vector<const RevBayesCore::TypedDagNode<double>* >::const_iterator it = pars.begin(); it != pars.end(); ++it)
        {
            rates.push_back( const_cast<RevBayesCore::StochasticNode<double>* >(static_cast<const RevBayesCore::StochasticNode<double>* >( *it ) ) );
        }
    }
    else
    {
        
    }
    
    value = new RevBayesCore::RateAgeBetaShift(t, rates, snode_rates, d, at, w);
}


/** Get Rev type of object */
const std::string& Move_RateAgeBetaShift::getClassType(void)
{
    
    static std::string rev_type = "Move_RateAgeBetaShift";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_RateAgeBetaShift::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev aliases for the constructor function.
 *
 * \return Rev aliasses of constructor function.
 */
std::vector<std::string> Move_RateAgeBetaShift::getMoveAliases(void) const
{
    std::vector<std::string> aliases;

    aliases.push_back("NodeRateTimeSlideBeta");

    return aliases;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_RateAgeBetaShift::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "RateAgeBetaShift";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_RateAgeBetaShift::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "tree" , Tree::getClassTypeSpec()                , "The tree on which this move operates on.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "rates", ModelVector<RealPos>::getClassTypeSpec(), "The vector of per-branch rates (from a relaxed clock).", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY)  );
        move_member_rules.push_back( new ArgumentRule( "delta", RealPos::getClassTypeSpec()             , "The concentration of the move on the previous age.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Real(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune" , RlBoolean::getClassTypeSpec()           , "Should we tune this move during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_RateAgeBetaShift::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_RateAgeBetaShift::printValue(std::ostream &o) const {
    
    o << "Move_RateAgeBetaShift(";
    if (tree != NULL) {
        o << tree->getName();
    }
    else {
        o << "?";
    }
    o << ")";
}


/** Set a NearestNeighborInterchange variable */
void Move_RateAgeBetaShift::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "tree" ) {
        tree = var;
    }
    else if ( name == "rates" ) {
        rates = var;
    }
    else if ( name == "delta" ) {
        delta = var;
    }
    else if ( name == "tune" ) {
        tune = var;
    }
    else {
        Move::setConstParameter(name, var);
    }
}
