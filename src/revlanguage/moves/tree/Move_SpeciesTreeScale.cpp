#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_SpeciesTreeScale.h"
#include "Natural.h"
#include "SpeciesTreeScaleProposal.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "Move.h"
#include "RbBoolean.h"
#include "RbHelpReference.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RlUtils.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_SpeciesTreeScale::Move_SpeciesTreeScale() : Move()
{

    // add method for call "addGeneTreeVariable" as a function
    ArgumentRules* addGeneTreeArgRules = new ArgumentRules();
    addGeneTreeArgRules->push_back( new ArgumentRule( "geneTree" , TimeTree::getClassTypeSpec(), "A gene tree variable.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
    methods.addFunction( new MemberProcedure( "addGeneTreeVariable", RlUtils::Void, addGeneTreeArgRules) );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_SpeciesTreeScale* Move_SpeciesTreeScale::clone(void) const
{

    return new Move_SpeciesTreeScale(*this);
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
void Move_SpeciesTreeScale::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( speciesTree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *st = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    RevBayesCore::TypedDagNode<double>* tmp_a = static_cast<const RealPos &>( rootAge->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<double> *ra = static_cast<RevBayesCore::StochasticNode<double> *>( tmp_a );
    double d = static_cast<const RealPos &>( delta->getRevObject() ).getValue();
    bool tune = static_cast<const RlBoolean &>( tuning->getRevObject() ).getValue();

    RevBayesCore::Proposal *p = new RevBayesCore::SpeciesTreeScaleProposal(st,ra,d);
    value = new RevBayesCore::MetropolisHastingsMove(p, w, del, tune);

}


RevPtr<RevVariable> Move_SpeciesTreeScale::executeMethod(const std::string& name, const std::vector<Argument>& args, bool &found)
{

    if ( name == "addGeneTreeVariable" )
    {
        found = true;

        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::Tree> *gt = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::SpeciesTreeScaleProposal &p = static_cast<RevBayesCore::SpeciesTreeScaleProposal&>( m->getProposal() );
        p.addGeneTree( gt );

        return NULL;
    }
    else if ( name == "removeGeneTreeVariable" )
    {
        found = true;

        RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<RevBayesCore::Tree> *gt = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

        RevBayesCore::MetropolisHastingsMove *m = static_cast<RevBayesCore::MetropolisHastingsMove*>(this->value);
        RevBayesCore::SpeciesTreeScaleProposal &p = static_cast<RevBayesCore::SpeciesTreeScaleProposal&>( m->getProposal() );
        p.removeGeneTree( gt );

        return NULL;
    }

    return Move::executeMethod( name, args, found );
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_SpeciesTreeScale::getClassType(void)
{

    static std::string rev_type = "Move_SpeciesTreeScale";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_SpeciesTreeScale::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_SpeciesTreeScale::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SpeciesTreeScale";

    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be a time-tree.
 *
 * \return The member rules.
 */
const MemberRules& Move_SpeciesTreeScale::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "speciesTree", TimeTree::getClassTypeSpec()             , "The species tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC    ) );
        memberRules.push_back( new ArgumentRule( "rootAge"    , RealPos::getClassTypeSpec()              , "The root age variable.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC    ) );
        memberRules.push_back( new ArgumentRule( "delta"      , RealPos::getClassTypeSpec()              , "The strength of the proposal", ArgumentRule::BY_VALUE    , ArgumentRule::ANY          , new RealPos( 1.0 ) ) );
        memberRules.push_back( new ArgumentRule( "tune"       , RlBoolean::getClassTypeSpec()            , "Should we tune the strength during burnin?", ArgumentRule::BY_VALUE    , ArgumentRule::ANY          , new RlBoolean( true ) ) );

        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Move_SpeciesTreeScale::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/**
 * Print the value for the user.
 */
void Move_SpeciesTreeScale::printValue(std::ostream &o) const
{

    o << "SpeciesTreeScale(";
    if (speciesTree != NULL)
    {
        o << speciesTree->getName();
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
void Move_SpeciesTreeScale::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "speciesTree" )
    {
        speciesTree = var;
    }
    else if ( name == "geneTrees" )
    {
        geneTrees = var;
    }
    else if ( name == "rootAge" )
    {
        rootAge = var;
    }
    else if ( name == "delta" )
    {
        delta = var;
    }
    else if ( name == "tune" )
    {
        tuning = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }

}
