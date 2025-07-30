#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "HomeologPhaseProposal.h"
#include "MetropolisHastingsMove.h"
#include "Move_HomeologPhase.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
RevLanguage::Move_HomeologPhase::Move_HomeologPhase() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_HomeologPhase* RevLanguage::Move_HomeologPhase::clone(void) const
{
    
    return new Move_HomeologPhase(*this);
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
void RevLanguage::Move_HomeologPhase::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // move/proposal arguments
    double w        = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    size_t del = static_cast<const Natural &>( delay->getRevObject() ).getValue();
    
    // move/proposal parameters
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn   = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
    std::string t1  = static_cast<const RlString &>( tip1->getRevObject() ).getValue();
    std::string t2  = static_cast<const RlString &>( tip2->getRevObject() ).getValue();
    
    // create the proposal
    RevBayesCore::HomeologPhaseProposal *p = new RevBayesCore::HomeologPhaseProposal(ctmc_sn, t1, t2);
    
    value = new RevBayesCore::MetropolisHastingsMove(p, w, del, false);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& RevLanguage::Move_HomeologPhase::getClassType(void)
{
    
    static std::string revType = "Move_HomeologPhase";
    
    return revType;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& RevLanguage::Move_HomeologPhase::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_HomeologPhase::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "HomeologPhase";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * \return The member rules.
 */
const MemberRules& RevLanguage::Move_HomeologPhase::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
        
        memberRules.push_back( new ArgumentRule( "ctmc", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The PhyloCTMC variable.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        memberRules.push_back( new ArgumentRule( "tip1", RlString::getClassTypeSpec(), "Tip 1 to switch.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule( "tip2", RlString::getClassTypeSpec(), "Tip 2 to switch.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        memberRules.insert( memberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rulesSet = true;
    }
    
    return memberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& RevLanguage::Move_HomeologPhase::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}

void RevLanguage::Move_HomeologPhase::printValue(std::ostream &o) const {
    
    o << "HomeologPhase()";
    
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
void RevLanguage::Move_HomeologPhase::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "ctmc" )
    {
        ctmc = var;
    }
    else if ( name == "tip1" )
    {
        tip1 = var;
    }
    else if ( name == "tip2" )
    {
        tip2 = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}
