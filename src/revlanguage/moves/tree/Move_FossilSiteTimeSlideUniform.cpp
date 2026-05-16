#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_FossilSiteTimeSlideUniform.h"
#include "FossilSiteTimeSlideUniformProposal.h"
#include "Probability.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlTaxon.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevNullObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "RlClade.h"
#include "RbVector.h"
#include "ModelVector.h"


using namespace RevLanguage;

Move_FossilSiteTimeSlideUniform::Move_FossilSiteTimeSlideUniform() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_FossilSiteTimeSlideUniform* Move_FossilSiteTimeSlideUniform::clone(void) const
{
    
    return new Move_FossilSiteTimeSlideUniform(*this);
}


void Move_FossilSiteTimeSlideUniform::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    double de = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    double we = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    bool   tu = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    double tt = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    
    RevBayesCore::TypedDagNode<double> *org = NULL;
    if ( origin != NULL && origin->getRevObject() != RevNullObject::getInstance() )
    {
        org = static_cast<const RealPos &>( origin->getRevObject() ).getDagNode();
    }
    
    // use NULL as defaults
    RevBayesCore::TypedDagNode<double> *ma = NULL;
    RevBayesCore::TypedDagNode<double> *mi = NULL;
    
    if ( max->getRevObject() != RevLanguage::RevNullObject::getInstance() )
    {
        ma = static_cast<const RealPos &>( max->getRevObject() ).getDagNode();
    }
    if ( min->getRevObject() != RevLanguage::RevNullObject::getInstance() )
    {
        mi = static_cast<const RealPos &>( min->getRevObject() ).getDagNode();
    }
    
    const RevBayesCore::Clade &c = static_cast<const Clade &>( clade->getRevObject() ).getValue();

    RevBayesCore::FossilSiteTimeSlideUniformProposal *p = new RevBayesCore::FossilSiteTimeSlideUniformProposal( t, org, ma, mi, c, de, tt );

    value = new RevBayesCore::MetropolisHastingsMove(p, we, tu);
}


/** Get Rev type of object */
const std::string& Move_FossilSiteTimeSlideUniform::getClassType(void)
{
    
    static std::string rev_type = "Move_FossilSiteTimeSlideUniform";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_FossilSiteTimeSlideUniform::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_FossilSiteTimeSlideUniform::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "FossilSiteTimeSlideUniform";
    
    return c_name;
}

std::vector<std::string> Move_FossilSiteTimeSlideUniform::getMoveAliases(void) const
{
    std::vector<std::string> aliases;
    
    aliases.push_back("FossilSiteSlideUniform");
    aliases.push_back("SiteTimeSlideUniform");

    return aliases;
}


/** Return member rules (no members) */
const MemberRules& Move_FossilSiteTimeSlideUniform::getParameterRules(void) const
{
    
    static MemberRules move_member_rules;
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        move_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "origin", RealPos::getClassTypeSpec() , "The variable for the origin of the process giving a maximum age.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL) );
        move_member_rules.push_back( new ArgumentRule( "max", RealPos::getClassTypeSpec() , "The variable for the maximum age of this fossil tip.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL) );
        move_member_rules.push_back( new ArgumentRule( "min", RealPos::getClassTypeSpec() , "The variable for the minimun age of this fossil tip.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL) );
        
        std::vector<TypeSpec> tip_index_arg_types;
        tip_index_arg_types.push_back( RlString::getClassTypeSpec() );
        tip_index_arg_types.push_back( Taxon::getClassTypeSpec() );
        move_member_rules.push_back( new ArgumentRule( "clade", Clade::getClassTypeSpec(), "A clade object with fossils from the same site", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        move_member_rules.push_back( new ArgumentRule( "delta" , RealPos::getClassTypeSpec()  , "The window size parameter.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune"  , RlBoolean::getClassTypeSpec(), "Should we tune the window size during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY       , new RlBoolean( true ) ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_FossilSiteTimeSlideUniform::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_FossilSiteTimeSlideUniform::printValue(std::ostream &o) const
{
    
    o << "Move_FossilSiteTimeSlideUniform(";
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
void Move_FossilSiteTimeSlideUniform::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if (name == "origin")
    {
        origin = var;
    }
    else if (name == "clade")
    {
        clade = var;
    }
    else if (name == "max")
    {
        max = var;
    }
    else if (name == "min")
    {
        min = var;
    }
    else if ( name == "delta" )
    {
        delta = var;
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
