#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_FNPRCharacterHistory.h"
#include "FixedNodeheightPruneAndRegraftCharacterHistoryProposal.h"
#include "RealPos.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlRateGenerator.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Move.h"
#include "PathRejectionSampleProposal.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "RnaState.h" // IWYU pragma: keep
#include "AminoAcidState.h" // IWYU pragma: keep
#include "StandardState.h" // IWYU pragma: keep

namespace RevBayesCore { class DnaState; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class RateGenerator; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Move_FNPRCharacterHistory::Move_FNPRCharacterHistory() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_FNPRCharacterHistory* Move_FNPRCharacterHistory::clone(void) const
{
    
    return new Move_FNPRCharacterHistory(*this);
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
void Move_FNPRCharacterHistory::constructInternalObject( void )
{
    
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    double w = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *n = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );
    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn   = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* qmap_tdn                             = static_cast<const RateGenerator&>( rate_generator->getRevObject() ).getDagNode();
    
    std::string mt = ctmc_sn->getValue().getDataType();
    RevBayesCore::Proposal *p = NULL;
    
    if (mt == "DNA")
    {
        RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::DnaState> *tmp_p = new RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::DnaState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    else if (mt == "RNA")
    {
        RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::RnaState> *tmp_p = new RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::RnaState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    else if (mt == "AA" || mt == "Protein")
    {
        RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::AminoAcidState> *tmp_p = new RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::AminoAcidState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    else if (mt == "Standard")
    {
        RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::FixedNodeheightPruneAndRegraftCharacterHistoryProposal<RevBayesCore::StandardState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
    
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Move_FNPRCharacterHistory::getClassType(void)
{
    
    static std::string rev_type = "Move_FNPRCharacterHistory";
    
    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Move_FNPRCharacterHistory::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_FNPRCharacterHistory::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "FNPRCharacterHistory";
    
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
const MemberRules& Move_FNPRCharacterHistory::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree variable on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "ctmc", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The ctmc.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "qmap", RateGenerator::getClassTypeSpec(), "The rate matrix.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
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
const TypeSpec& Move_FNPRCharacterHistory::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/**
 * Print the value for the user.
 */
void Move_FNPRCharacterHistory::printValue(std::ostream &o) const
{
    
    o << "FNPRCharacterHistory(";
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
void Move_FNPRCharacterHistory::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "ctmc" )
    {
        ctmc = var;
    }
    else if ( name == "qmap")
    {
        rate_generator = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}




