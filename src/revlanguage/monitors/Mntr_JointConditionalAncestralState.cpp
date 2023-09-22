#include "Mntr_JointConditionalAncestralState.h"

#include <cstddef>
#include <string>

#include "BinaryState.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "IntegerPos.h"
#include "RbException.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlTimeTree.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "NaturalNumbersState.h"
#include "DnaState.h"
#include "StandardState.h"
#include "AminoAcidState.h"
#include "PoMoState.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "DiscreteTaxonData.h"
#include "JointConditionalAncestralStateMonitor.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Mntr_JointConditionalAncestralState::Mntr_JointConditionalAncestralState(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_JointConditionalAncestralState* Mntr_JointConditionalAncestralState::clone(void) const
{
    
    return new Mntr_JointConditionalAncestralState(*this);
}


void Mntr_JointConditionalAncestralState::constructInternalObject( void )
{
    const std::string&                  fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned int                        g       = (int)static_cast<const IntegerPos  &>( printgen->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* t = static_cast<const Tree &>( tree->getRevObject() ).getDagNode();
    
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                wt      = static_cast<const RlBoolean &>( withTips->getRevObject() ).getValue();
    bool                                wss     = static_cast<const RlBoolean &>( withStartStates->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    std::string                         character = static_cast<const RlString &>( monitorType->getRevObject() ).getValue();
    

    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn = NULL;
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* cdbdp_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::Tree>* cdbdp_sn = NULL;
    
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* glhbdsp_tdn = NULL;
    RevBayesCore::StochasticNode<RevBayesCore::Tree>* glhbdsp_sn = NULL;

    if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() )
    {
        ctmc_tdn = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
        ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
        
        if ( ctmc_sn->getValue().getDataType() != character )
        {
            throw RbException("mnJointConditionalAncestralStateMonitor requires the CTMC to be of same type as the specified character.");
        }
    }
    else if ( static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).isModelObject() )
    {
        cdbdp_tdn = static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).getDagNode();
        cdbdp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( cdbdp_tdn );
    }
    else if ( static_cast<const RevLanguage::Tree&>( glhbdsp->getRevObject() ).isModelObject() )
    {
    	glhbdsp_tdn = static_cast<const RevLanguage::Tree&>( glhbdsp->getRevObject() ).getDagNode();
    	glhbdsp_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree>* >( glhbdsp_tdn );
    }
    else
    {
        throw RbException("mnJointConditionalAncestralStateMonitor requires either a CTMC or a character-dependent birth death process (CDBDP).");
    }
    
    
    delete value;
    if (character == "AA" || character == "Protein")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::AminoAcidState> *m;
        m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::AminoAcidState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "DNA")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::DnaState> *m;
        m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::DnaState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "NaturalNumbers")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::NaturalNumbersState> *m;
        if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() )
        {
            m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::NaturalNumbersState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        }
        else if ( static_cast<const RevLanguage::Tree&>( cdbdp->getRevObject() ).isModelObject() )
        {
            m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::NaturalNumbersState>(cdbdp_sn, (unsigned long)g, fn, sep, wt, wss);
        }
        else
        {
        	m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::NaturalNumbersState>(glhbdsp_sn, (unsigned long)g, fn, sep, wt, wss);
        }
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "PoMo")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::PoMoState> *m;
        m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::PoMoState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "RNA")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::DnaState> *m;
        m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::DnaState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "Standard")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::StandardState> *m;
        if ( static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).isModelObject() )
        {
            m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::StandardState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        }
        else
        {
            m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::StandardState>(cdbdp_sn, (unsigned long)g, fn, sep, wt, wss);
        }
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "Binary" || character == "Restriction")
    {
        RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::BinaryState> *m;
        m = new RevBayesCore::JointConditionalAncestralStateMonitor<RevBayesCore::BinaryState>(t, ctmc_sn, (unsigned long)g, fn, sep, wt, wss);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else
    {
        throw RbException( "Incorrect character type specified. Valid options are: AA, DNA, NaturalNumbers, PoMo, Protein, RNA, Standard, Binary/Restriction" );
    }
    
}


/** Get Rev type of object */
const std::string& Mntr_JointConditionalAncestralState::getClassType(void)
{
    
    static std::string rev_type = "Mntr_JointConditionalAncestralState";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_JointConditionalAncestralState::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_JointConditionalAncestralState::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "JointConditionalAncestralState";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_JointConditionalAncestralState::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        memberRules.push_back( new ArgumentRule("tree"           , Tree::getClassTypeSpec() , "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("ctmc"           , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        memberRules.push_back( new ArgumentRule("cdbdp"          , TimeTree::getClassTypeSpec(), "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL) );
        memberRules.push_back( new ArgumentRule("glhbdsp"        , TimeTree::getClassTypeSpec(), "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL) );
        memberRules.push_back( new ArgumentRule("type"           , RlString::getClassTypeSpec() , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("withTips"       , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("withStartStates", RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Mntr_JointConditionalAncestralState::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_JointConditionalAncestralState::printValue(std::ostream &o) const
{
    
    o << "Mntr_JointConditionalAncestralState";
}


/** Set a member variable */
void Mntr_JointConditionalAncestralState::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "type" )
    {
        monitorType = var;
    }
    else if ( name == "ctmc" )
    {
        ctmc = var;
    }
    else if ( name == "cdbdp" )
    {
        cdbdp = var;
    }
    else if ( name == "glhbdsp" )
    {
    	glhbdsp = var;
    }
    else if ( name == "withTips" )
    {
        withTips = var;
    }
    else if ( name == "withStartStates" )
    {
        withStartStates = var;
    }
    else
    {
        FileMonitor::setConstParameter(name, var);
    }
    
}
