/* 
 * File:   Move_CorrelationMatrixSingleElementBeta.cpp
 * Author: Michael R. May
 *
 * Created on 5 August 2017
 */

#include "Move_CorrelationMatrixSingleElementBeta.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "RlBoolean.h"
#include "MatrixReal.h"
#include "CorrelationMatrixElementBetaProposal.h"
#include "MetropolisHastingsMove.h"
#include "Natural.h"
#include "Probability.h"
#include "RbException.h"
#include "Real.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlMatrixReal.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"


using namespace RevLanguage;

Move_CorrelationMatrixSingleElementBeta::Move_CorrelationMatrixSingleElementBeta() : Move()
{
    
}

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_CorrelationMatrixSingleElementBeta* Move_CorrelationMatrixSingleElementBeta::clone(void) const
{
    
	return new Move_CorrelationMatrixSingleElementBeta(*this);
}


void Move_CorrelationMatrixSingleElementBeta::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double a = static_cast<const RealPos &>( alpha->getRevObject() ).getValue();
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double r = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal >* tmp = static_cast<const MatrixReal &>( v->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal > *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal> *>( tmp );
    bool t = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    
    RevBayesCore::Proposal *p = new RevBayesCore::CorrelationMatrixElementBetaProposal(n,a,r);
    value = new RevBayesCore::MetropolisHastingsMove(p,w,t);

}


/** Get class name of object */
const std::string& Move_CorrelationMatrixSingleElementBeta::getClassType(void) { 
    
    static std::string revClassType = "Move_CorrelationMatrixSingleElementBeta";
    
	return revClassType; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_CorrelationMatrixSingleElementBeta::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return revClassTypeSpec; 
}


/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Move_CorrelationMatrixSingleElementBeta::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Michael R. May" );
    
    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Move_CorrelationMatrixSingleElementBeta::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string descriptions = "";
    descriptions += "Beta proposal on a random element of a correlation matrix.";
    
    return descriptions;
}


/**
 * Get the more detailed description of the function
 */
std::string Move_CorrelationMatrixSingleElementBeta::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "";
    details += "This move chooses a single element of the correlation matrix at random, and draws a proposed value from a Beta distribution centered on the current value (and stretched to range from -1 to 1).";
    
    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Move_CorrelationMatrixSingleElementBeta::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "\n";
    
    example += "# we draw a correlation matrix from an LKJ distribution\n";
    example += "R ~ dnLKJ(eta=1, dim=5)\n\n";
    example += "# we specify a beta move on the correlation matrix\n";
    example += "moves[1] = mvCorrelationMatrixSingleElementBeta(R, alpha=10.0)\n";
    
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Move_CorrelationMatrixSingleElementBeta::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;
    
    
    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Move_CorrelationMatrixSingleElementBeta::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;
    see_also.push_back( "mvCorrelationMatrixSpecificElementBeta" );
    see_also.push_back( "mvCorrelationMatrixRandomWalk" );
    
    
    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Move_CorrelationMatrixSingleElementBeta::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Correlation Matrix Beta proposal.";
    
    return title;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_CorrelationMatrixSingleElementBeta::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CorrelationMatrixSingleElementBeta";
    
    return c_name;
}



/** Return member rules (no members) */
const MemberRules& Move_CorrelationMatrixSingleElementBeta::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        move_member_rules.push_back( new ArgumentRule( "x"     , MatrixReal::getClassTypeSpec(), "The variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "alpha",  RealPos::getClassTypeSpec()   , "The shape of the beta distribution to draw from.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Real(10.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec() , "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_CorrelationMatrixSingleElementBeta::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Move_CorrelationMatrixSingleElementBeta::printValue(std::ostream &o) const
{
    
    o << "Move_CorrelationMatrixSingleElementBeta(";
    if (v != NULL)
    {
        o << v->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_CorrelationMatrixSingleElementBeta::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "x" )
    {
        v = var;
    }
    else if ( name == "alpha" )
    {
        alpha = var;
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

