#include "RevVariableUtilities.h"
#include <cstddef>
#include "TypeUtilities.h"
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "DagNodeTypeUtilities.h"
#include "ArgumentRules.h"
//#include "TreeCharacterHistoryNhxMonitor.h"
#include "TreeCharacterHistoryNodeMonitor.h"
#include "Mntr_CharacterHistoryNewickFile.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "IntegerPos.h"
#include "RevObject.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "StochasticNode.h"
#include "Real.h" // IWYU pragma: keep
#include "StandardState.h" // IWYU pragma: keep

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class RbVector; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Mntr_CharacterHistoryNewickFile::Mntr_CharacterHistoryNewickFile(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_CharacterHistoryNewickFile* Mntr_CharacterHistoryNewickFile::clone(void) const
{
    
	return new Mntr_CharacterHistoryNewickFile(*this);
}


void Mntr_CharacterHistoryNewickFile::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const std::string& fn = valueFor<RlString>(filename);
    const std::string& sep = valueFor<RlString>(separator);
    unsigned int g = (int)valueFor<Natural>(printgen);
   
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *t = dagNodeFor<TimeTree>(tree);

    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn   = dagNodeFor<RevLanguage::AbstractHomologousDiscreteCharacterData>(ctmc);
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn  = RevBayesCore::assumeStochastic<RevBayesCore::AbstractHomologousDiscreteCharacterData>( ctmc_tdn );
    
    bool pp = valueFor<RlBoolean>(posterior);
    bool l = valueFor<RlBoolean>(likelihood);
    bool pr = valueFor<RlBoolean>(prior);
//    bool sc = static_cast<const RlBoolean &>( counts->getRevObject() ).getValue();
//    bool se = static_cast<const RlBoolean &>( events->getRevObject() ).getValue();
    bool sm = true; // show metadata disabled for now
    bool ap = false; //static_cast<const RlBoolean &>( append->getRevObject() ).getValue();

    
    std::string ms = valueFor<RlString>(style);
    bool sc = false;
    if (ms == "counts")
        sc = true;
    bool se = !sc;
    
    std::string mt = valueFor<RlString>(type);
    if (mt == "std")
        ; // value = XXXXXX
    else if (mt == "biogeo")
        value = new RevBayesCore::TreeCharacterHistoryNodeMonitor<RevBayesCore::StandardState>(ctmc_sn, t, size_t(g), fn, sep, pp, l, pr, ap, sm, sc, se);
}


/** Get class name of object */
const std::string& Mntr_CharacterHistoryNewickFile::getClassType(void)
{
    
    static std::string revClassType = "Mntr_CharacterHistoryNewickFile";
    
	return revClassType;
}

/** Get class type spec describing type of object */
const TypeSpec& Mntr_CharacterHistoryNewickFile::getClassTypeSpec(void)
{
    
    static TypeSpec revClassTypeSpec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
	return revClassTypeSpec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_CharacterHistoryNewickFile::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CharHistoryNewick";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_CharacterHistoryNewickFile::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        memberRules.push_back( new ArgumentRule("ctmc"      , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("tree"      , TimeTree::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("posterior" , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("likelihood", RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("prior"     , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );

        std::vector<std::string> options_style;
        options_style.push_back( "events" );
        options_style.push_back( "counts" );
        memberRules.push_back( new OptionRule( "style", new RlString("events"), options_style, "" ) );

        
        std::vector<std::string> options;
        options.push_back( "biogeo" );
        memberRules.push_back( new OptionRule( "type", new RlString("biogeo"), options, "" ) );

        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());
        
        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Mntr_CharacterHistoryNewickFile::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_CharacterHistoryNewickFile::printValue(std::ostream &o) const
{
    
    o << "Mntr_CharacterHistoryNewickFile";
}


/** Set a member variable */
void Mntr_CharacterHistoryNewickFile::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "ctmc" )
    {
        ctmc = var;
    }
    else if ( name == "prior" )
    {
        prior = var;
    }
    else if ( name == "posterior" )
    {
        posterior = var;
    }
    else if ( name == "likelihood" )
    {
        likelihood = var;
    }
    else if ( name == "type" )
    {
        type = var;
    }
    else if ( name == "style" )
    {
        style = var;
    }
    else
    {
        FileMonitor::setConstParameter(name, var);
    }
}
