#include <cstddef>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
//#include "BiogeographicTreeHistoryCtmc.h"
#include "BiogeographyNodeRejectionSampleProposal.h"
#include "BiogeographyCladogeneticRejectionShiftProposal.h"
#include "BiogeographyCladogeneticRejectionSampleProposal.h"
//#include "BiogeographyPathRejectionSampleProposal.h"
#include "NodeRejectionSampleProposal.h"
#include "MetropolisHastingsMove.h"
#include "Move_CharacterHistory.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlRateGenerator.h"
#include "RlRateGeneratorSequence.h"
#include "RlString.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "Move.h"
#include "PathRejectionSampleProposal.h"
#include "RevNullObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "RnaState.h" // IWYU pragma: keep
#include "AminoAcidState.h" // IWYU pragma: keep
#include "StandardState.h" // IWYU pragma: keep

namespace RevBayesCore { class DnaState; }
namespace RevBayesCore { class Proposal; }
namespace RevBayesCore { class RateGenerator; }
namespace RevBayesCore { class RateGeneratorSequence; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
RevLanguage::Move_CharacterHistory::Move_CharacterHistory() : Move()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the move.
 */
Move_CharacterHistory* RevLanguage::Move_CharacterHistory::clone(void) const
{
    
    return new Move_CharacterHistory(*this);
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
void RevLanguage::Move_CharacterHistory::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // move/proposal arguments
    double w        = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double l        = static_cast<const Probability &>( lambda->getRevObject() ).getValue();
    std::string gt  = static_cast<const RlString &>( graph->getRevObject() ).getValue();
    std::string pt  = static_cast<const RlString &>( proposal->getRevObject() ).getValue();
    double r        = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    
    // move/proposal parameters
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn   = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
    
    bool use_site = false;
    bool use_seq = false;
    
    RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator>* qmap_site_tdn = NULL;
    RevBayesCore::TypedDagNode<RevBayesCore::RateGeneratorSequence>* qmap_seq_tdn = NULL;
    //    if (qmap_site != NULL ) {
    if (qmap_site->getRevObject() != RevNullObject::getInstance()) {
        use_site = true;
        qmap_site_tdn = static_cast<const RateGenerator&>( qmap_site->getRevObject() ).getDagNode();
    }
    else if (qmap_seq->getRevObject() != RevNullObject::getInstance()) {
        use_seq = true;
        qmap_seq_tdn = static_cast<const RateGeneratorSequence&>( qmap_seq->getRevObject() ).getDagNode();
    }
    else {
        throw RbException("qmap_site or qmap_seq must be provided!");
    }
    
    // get data type
    std::string mt  = ctmc_tdn->getValue().getDataType();
    
    // create the proposal
    RevBayesCore::Proposal *p = NULL;
    
    if (mt == "DNA")
    {
        if (gt == "node" && pt == "rejection")
        {
            RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::DnaState> *tmp_p = new RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::DnaState>(ctmc_sn, l, r);
            
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "branch" && pt == "rejection")
        {
            RevBayesCore::PathRejectionSampleProposal<RevBayesCore::DnaState> *tmp_p = new RevBayesCore::PathRejectionSampleProposal<RevBayesCore::DnaState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
    }
    else if (mt == "RNA")
    {
        if (gt == "node" && pt == "rejection")
        {
            RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::RnaState> *tmp_p = new RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::RnaState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "branch" && pt == "rejection")
        {
            RevBayesCore::PathRejectionSampleProposal<RevBayesCore::RnaState> *tmp_p = new RevBayesCore::PathRejectionSampleProposal<RevBayesCore::RnaState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
    }
    else if (mt == "AA" || mt == "Protein")
    {
        if (gt == "node" && pt == "rejection")
        {
            RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::AminoAcidState> *tmp_p = new RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::AminoAcidState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "branch" && pt == "rejection")
        {
            RevBayesCore::PathRejectionSampleProposal<RevBayesCore::AminoAcidState> *tmp_p = new RevBayesCore::PathRejectionSampleProposal<RevBayesCore::AminoAcidState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn);
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
    }
    else if (mt == "Standard")
    {
        if (gt == "node" && pt == "rejection")
        {
            RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::NodeRejectionSampleProposal<RevBayesCore::StandardState>(ctmc_sn, l, r);
            
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "cladogenetic" && pt == "rejection")
        {
            RevBayesCore::BiogeographicCladogeneticRejectionSampleProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::BiogeographicCladogeneticRejectionSampleProposal<RevBayesCore::StandardState>(ctmc_sn, l, r);
            
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "cladogenetic2" && pt == "rejection")
        {
            RevBayesCore::BiogeographicNodeRejectionSampleProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::BiogeographicNodeRejectionSampleProposal<RevBayesCore::StandardState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "cladogenetic" && pt == "rejection_shift")
        {
            RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::BiogeographicCladogeneticRejectionShiftProposal<RevBayesCore::StandardState>(ctmc_sn, l, r);
            //            tmp_p->setRateGenerator( qmap_tdn );
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
        else if (gt == "branch" && pt == "rejection")
        {
            RevBayesCore::PathRejectionSampleProposal<RevBayesCore::StandardState> *tmp_p = new RevBayesCore::PathRejectionSampleProposal<RevBayesCore::StandardState>(ctmc_sn, l, r);
            
            if (use_site) {
                tmp_p->setRateGenerator( qmap_site_tdn );
            } else if (use_seq) {
                tmp_p->setRateGenerator( qmap_seq_tdn );
            }
            p = tmp_p;
        }
    }
    
    
    
    value = new RevBayesCore::MetropolisHastingsMove(p,w,false);
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& RevLanguage::Move_CharacterHistory::getClassType(void)
{
    
    static std::string revType = "Move_CharacterHistory"; // <" + treeType::getClassType() + ">";
    
    return revType;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& RevLanguage::Move_CharacterHistory::getClassTypeSpec(void)
{
    
    static TypeSpec revTypeSpec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );
    
    return revTypeSpec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Move_CharacterHistory::getMoveName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CharacterHistory";
    
    return c_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the scale move are:
 * (1) the variable which must be a positive real.
 * (2) the tuning parameter lambda that defines the size of the proposal (positive real)
 * (3) a flag whether auto-tuning should be used.
 *
 * \return The member rules.
 */
const MemberRules& RevLanguage::Move_CharacterHistory::getParameterRules(void) const
{
    
    static MemberRules nodeChrsMoveMemberRules;
    static bool rulesSet = false;
    
    if ( !rulesSet )
    {
        
        nodeChrsMoveMemberRules.push_back( new ArgumentRule( "ctmc", AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "The PhyloCTMC variable.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        nodeChrsMoveMemberRules.push_back( new ArgumentRule( "qmap_site", RateGenerator::getClassTypeSpec(),         "Per-site rate generator.",     ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        nodeChrsMoveMemberRules.push_back( new ArgumentRule( "qmap_seq",  RateGeneratorSequence::getClassTypeSpec(), "Per-sequence rate generator.", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        nodeChrsMoveMemberRules.push_back( new ArgumentRule( "lambda", Probability::getClassTypeSpec(), "Tuning probability to propose new site history.", ArgumentRule::BY_VALUE    , ArgumentRule::ANY, new Probability(1.0) ) );
        
        //        std::vector<std::string> optionsType;
        //        optionsType.push_back( "Biogeo" );
        //        optionsType.push_back( "DNA" );
        //        optionsType.push_back( "RNA" );
        //        optionsType.push_back( "AA" );
        //        optionsType.push_back( "Protein" );
        //        optionsType.push_back( "Standard" );
        //        nodeChrsMoveMemberRules.push_back( new OptionRule( "type", new RlString("Standard"), optionsType, "The data type." ) );
        
        std::vector<std::string> optionsGraph;
        optionsGraph.push_back( "node" );
        optionsGraph.push_back( "branch" );
        optionsGraph.push_back( "cladogenetic" );
        optionsGraph.push_back( "cladogenetic2" );
        nodeChrsMoveMemberRules.push_back( new OptionRule( "graph", new RlString("node"), optionsGraph, "" ) );
        
        std::vector<std::string> optionsProposal;
        optionsProposal.push_back( "rejection" );
        optionsProposal.push_back( "rejection_shift" );
        optionsProposal.push_back( "uniformization" );
        nodeChrsMoveMemberRules.push_back( new OptionRule( "proposal", new RlString("rejection"), optionsProposal, "" ) );
        
        /* Inherit weight from Move, put it after variable */
        
        /* Inherit weight from Move, put it after variable */
        const MemberRules& inheritedRules = Move::getParameterRules();
        nodeChrsMoveMemberRules.insert( nodeChrsMoveMemberRules.end(), inheritedRules.begin(), inheritedRules.end() );
        
        rulesSet = true;
    }
    
    return nodeChrsMoveMemberRules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& RevLanguage::Move_CharacterHistory::getTypeSpec( void ) const
{
    
    static TypeSpec typeSpec = getClassTypeSpec();
    
    return typeSpec;
}

void RevLanguage::Move_CharacterHistory::printValue(std::ostream &o) const {
    
    o << "CharacterHistory(";
    if (qmap_site != NULL)
    {
        o << qmap_site->getName();
    }
    else if (qmap_seq != NULL)
    {
        o << qmap_seq->getName();
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
void RevLanguage::Move_CharacterHistory::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "ctmc" )
    {
        ctmc = var;
    }
    else if ( name == "qmap_site" )
    {
        qmap_site = var;
    }
    else if ( name == "qmap_seq" )
    {
        qmap_seq = var;
    }
    else if ( name == "graph" )
    {
        graph = var;
    }
    else if ( name == "proposal" )
    {
        proposal = var;
    }
    else if ( name == "lambda" )
    {
        lambda = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
    
}
