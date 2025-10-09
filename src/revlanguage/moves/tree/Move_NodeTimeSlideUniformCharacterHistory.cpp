#include <cstddef>
#include <ostream>
#include <string>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_NodeTimeSlideUniformCharacterHistory.h"
#include "NodeTimeSlideUniformCharacterHistoryProposal.h"
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

Move_NodeTimeSlideUniformCharacterHistory::Move_NodeTimeSlideUniformCharacterHistory() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Move_NodeTimeSlideUniformCharacterHistory* Move_NodeTimeSlideUniformCharacterHistory::clone(void) const
{
    
	return new Move_NodeTimeSlideUniformCharacterHistory(*this);
}


void Move_NodeTimeSlideUniformCharacterHistory::constructInternalObject( void )
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
        RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::DnaState> *tmp_p = new RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::DnaState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    else if (mt == "RNA")
    {
        RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::RnaState> *tmp_p = new RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::RnaState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    else if (mt == "AA" || mt == "Protein")
    {
        RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::AminoAcidState> *tmp_p = new RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::AminoAcidState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    else if (mt == "Standard")
    {
        RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::NodeTimeSlideUniformCharacterHistoryProposal<RevBayesCore::StandardState>(n,ctmc_sn);
        tmp_p->setRateGenerator( qmap_tdn );
        p = tmp_p;
    }
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w);
}


/** Get Rev type of object */
const std::string& Move_NodeTimeSlideUniformCharacterHistory::getClassType(void)
{
    
    static std::string rev_type = "Move_NodeTimeSlideUniformCharacterHistory";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Move_NodeTimeSlideUniformCharacterHistory::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_NodeTimeSlideUniformCharacterHistory::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "NodeTimeSlideUniformCharacterHistory";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_NodeTimeSlideUniformCharacterHistory::getParameterRules(void) const
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

/** Get type spec */
const TypeSpec& Move_NodeTimeSlideUniformCharacterHistory::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



/** Get type spec */
void Move_NodeTimeSlideUniformCharacterHistory::printValue(std::ostream &o) const
{
    
    o << "Move_NodeTimeSlideUniformCharacterHistory(";
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
void Move_NodeTimeSlideUniformCharacterHistory::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
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
