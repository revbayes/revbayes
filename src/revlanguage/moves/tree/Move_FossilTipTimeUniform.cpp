#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_FossilTipTimeUniform.h"
#include "Natural.h"
#include "FossilTipTimeUniformProposal.h"
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


using namespace RevLanguage;

Move_FossilTipTimeUniform::Move_FossilTipTimeUniform() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_FossilTipTimeUniform* Move_FossilTipTimeUniform::clone(void) const
{
    
    return new Move_FossilTipTimeUniform(*this);
}


void Move_FossilTipTimeUniform::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    double we = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    
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
    
    std::string tip_name = "";
    if ( tip->getRevObject().getType() == RlString::getClassType() )
    {
        tip_name = static_cast<const RlString&>( tip->getRevObject() ).getValue();
    }
    else if ( tip->getRevObject().getType() == Taxon::getClassType() )
    {
        tip_name = static_cast<const Taxon&>( tip->getRevObject() ).getValue().getSpeciesName();
    }

    RevBayesCore::FossilTipTimeUniformProposal *p = new RevBayesCore::FossilTipTimeUniformProposal( t, org, ma, mi, tip_name );
    
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    value = new RevBayesCore::MetropolisHastingsMove(p, we, del, false);
}


/** Get Rev type of object */
const std::string& Move_FossilTipTimeUniform::getClassType(void)
{
    
    static std::string rev_type = "Move_FossilTipTimeUniform";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Move_FossilTipTimeUniform::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_FossilTipTimeUniform::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "FossilTipTimeUniform";
    
    return c_name;
}

std::vector<std::string> Move_FossilTipTimeUniform::getMoveAliases(void) const
{
    std::vector<std::string> aliases;

    return aliases;
}


/** Return member rules (no members) */
const MemberRules& Move_FossilTipTimeUniform::getParameterRules(void) const
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
        move_member_rules.push_back( new ArgumentRule( "tip", tip_index_arg_types, "The name of a specific tip/taxon.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RlString("") ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rules_set = true;
    }
    
    return move_member_rules;
}

/** Get type spec */
const TypeSpec& Move_FossilTipTimeUniform::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_FossilTipTimeUniform::printValue(std::ostream &o) const
{
    
    o << "Move_FossilTipTimeUniform(";
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
void Move_FossilTipTimeUniform::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if (name == "origin")
    {
        origin = var;
    }
    else if (name == "tip")
    {
        tip = var;
    }
    else if (name == "max")
    {
        max = var;
    }
    else if (name == "min")
    {
        min = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
