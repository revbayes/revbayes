#include "Mntr_SiteMixtureAllocation.h"

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
#include "CodonState.h"
#include "StandardState.h"
#include "AminoAcidState.h"
#include "PoMoState.h"
#include "AbstractHomologousDiscreteCharacterData.h"
#include "DiscreteTaxonData.h"
#include "SiteMixtureAllocationMonitor.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RlBoolean.h"
#include "RlTree.h"
#include "StochasticNode.h"

namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Mntr_SiteMixtureAllocation::Mntr_SiteMixtureAllocation(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_SiteMixtureAllocation* Mntr_SiteMixtureAllocation::clone(void) const
{
    
    return new Mntr_SiteMixtureAllocation(*this);
}


void Mntr_SiteMixtureAllocation::constructInternalObject( void )
{
    const std::string&                  fn      = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string&                  sep     = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned int                        g       = (int)static_cast<const IntegerPos  &>( printgen->getRevObject() ).getValue();
    
    bool                                ap      = static_cast<const RlBoolean &>( append->getRevObject() ).getValue();
    bool                                wv      = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    std::string                         character = static_cast<const RlString &>( monitorType->getRevObject() ).getValue();
    

    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
    
    delete value;
    if (character == "AA" || character == "Protein")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::AminoAcidState> *m;
        m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::AminoAcidState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "DNA")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::DnaState> *m;
        m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::DnaState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "Codon")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::CodonState> *m;
        m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::CodonState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "NaturalNumbers")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::NaturalNumbersState> *m;
		m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::NaturalNumbersState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "PoMo")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::PoMoState> *m;
        m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::PoMoState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "RNA")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::DnaState> *m;
        m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::DnaState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "Standard")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::StandardState> *m;
		m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::StandardState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else if (character == "Binary" || character == "Restriction")
    {
        RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::BinaryState> *m;
        m = new RevBayesCore::SiteMixtureAllocationMonitor<RevBayesCore::BinaryState>(ctmc_sn, (std::uint64_t)g, fn, sep);
        m->setAppend( ap );
        m->setPrintVersion(wv);
        value = m;
    }
    else
    {
        throw RbException( "Incorrect character type specified. Valid options are: AA, DNA, Codon, NaturalNumbers, PoMo, Protein, RNA, Standard, Binary/Restriction" );
    }
    
}


/** Get Rev type of object */
const std::string& Mntr_SiteMixtureAllocation::getClassType(void)
{
    
    static std::string rev_type = "Mntr_SiteMixtureAllocation";
    
    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Mntr_SiteMixtureAllocation::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_SiteMixtureAllocation::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "SiteMixtureAllocation";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_SiteMixtureAllocation::getParameterRules(void) const
{
    
    static MemberRules asMonitorMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        asMonitorMemberRules.push_back( new ArgumentRule("ctmc"           , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "", ArgumentRule::BY_REFERENCE, ArgumentRule::ANY, NULL ) );
        asMonitorMemberRules.push_back( new ArgumentRule("type"           , RlString::getClassTypeSpec() , "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        asMonitorMemberRules.insert(asMonitorMemberRules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }
    
    return asMonitorMemberRules;
}

/** Get type spec */
const TypeSpec& Mntr_SiteMixtureAllocation::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_SiteMixtureAllocation::printValue(std::ostream &o) const
{
    
    o << "Mntr_SiteMixtureAllocation";
}


/** Set a member variable */
void Mntr_SiteMixtureAllocation::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "type" )
    {
        monitorType = var;
    }
    else if ( name == "ctmc" )
    {
        ctmc = var;
    }
    else
    {
        FileMonitor::setConstParameter(name, var);
    }
    
}
