#include <stddef.h>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
//#include "TreeCharacterHistoryNhxMonitor.h"
#include "TreeCharacterHistoryNodeMonitor.h"
#include "Mntr_CharacterHistoryNewickFile.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "Natural.h"
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
#include "RlMonitor.h"
#include "StochasticNode.h"
#include "Real.h" // IWYU pragma: keep
#include "StandardState.h" // IWYU pragma: keep

namespace RevBayesCore { class AbstractHomologousDiscreteCharacterData; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class RbVector; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevLanguage;

Mntr_CharacterHistoryNewickFile::Mntr_CharacterHistoryNewickFile(void) : Monitor() {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_CharacterHistoryNewickFile* Mntr_CharacterHistoryNewickFile::clone(void) const {
    
	return new Mntr_CharacterHistoryNewickFile(*this);
}


void Mntr_CharacterHistoryNewickFile::constructInternalObject( void ) {
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const std::string& fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string& sep = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned int g = (int)static_cast<const Natural &>( printgen->getRevObject() ).getValue();
   
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *t = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    std::set<RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> > *> n;
    for (std::set<RevPtr<const RevVariable> >::iterator i = vars.begin(); i != vars.end(); ++i) {
        RevBayesCore::TypedDagNode< RevBayesCore::RbVector<double> >* node = static_cast< const ModelVector<Real> & >((*i)->getRevObject()).getDagNode();
        n.insert( node );
    }
    
    RevBayesCore::TypedDagNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_tdn   = static_cast<const RevLanguage::AbstractHomologousDiscreteCharacterData&>( ctmc->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* ctmc_sn  = static_cast<RevBayesCore::StochasticNode<RevBayesCore::AbstractHomologousDiscreteCharacterData>* >(ctmc_tdn);
    
    bool pp = static_cast<const RlBoolean &>( posterior->getRevObject() ).getValue();
    bool l = static_cast<const RlBoolean &>( likelihood->getRevObject() ).getValue();
    bool pr = static_cast<const RlBoolean &>( prior->getRevObject() ).getValue();
//    bool sc = static_cast<const RlBoolean &>( counts->getRevObject() ).getValue();
//    bool se = static_cast<const RlBoolean &>( events->getRevObject() ).getValue();
    bool sm = true; // show metadata disabled for now
    bool ap = false; //static_cast<const RlBoolean &>( append->getRevObject() ).getValue();

    
    std::string ms = static_cast<const RlString&>( style->getRevObject() ).getValue();
    bool sc = false;
    if (ms == "counts")
        sc = true;
    bool se = !sc;
    
    std::string mt = static_cast<const RlString&>( type->getRevObject() ).getValue();
    if (mt == "std")
        ; // value = XXXXXX
    else if (mt == "biogeo")
        value = new RevBayesCore::TreeCharacterHistoryNodeMonitor<RevBayesCore::StandardState>(ctmc_sn, t, size_t(g), fn, sep, pp, l, pr, ap, sm, sc, se);
}


/** Get class name of object */
const std::string& Mntr_CharacterHistoryNewickFile::getClassType(void) {
    
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
    
    static MemberRules Mntr_CharacterHistoryNewickFileMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("filename"  , RlString::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("ctmc"      , AbstractHomologousDiscreteCharacterData::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("tree"      , TimeTree::getClassTypeSpec(), "", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("printgen"  , IntegerPos::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(1) ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("separator" , RlString::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("\t") ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("posterior" , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("likelihood", RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("prior"     , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
//        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("counts", true, RlBoolean::getClassTypeSpec(), new RlBoolean(false) ) );
//        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("events", true, RlBoolean::getClassTypeSpec(), new RlBoolean(true) ) );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new ArgumentRule("append"    , RlBoolean::getClassTypeSpec(), "", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        std::vector<std::string> options_style;
        //        options.push_back( RlString("std") );
        options_style.push_back( "events" );
        options_style.push_back( "counts" );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new OptionRule( "style", new RlString("events"), options_style, "" ) );

        
        std::vector<std::string> options;
//        options.push_back( RlString("std") );
        options.push_back( "biogeo" );
        Mntr_CharacterHistoryNewickFileMemberRules.push_back( new OptionRule( "type", new RlString("biogeo"), options, "" ) );

        
        
        rules_set = true;
    }
    
    return Mntr_CharacterHistoryNewickFileMemberRules;
}

/** Get type spec */
const TypeSpec& Mntr_CharacterHistoryNewickFile::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_CharacterHistoryNewickFile::printValue(std::ostream &o) const {
    
    o << "Mntr_CharacterHistoryNewickFile";
}


/** Set a member variable */
void Mntr_CharacterHistoryNewickFile::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "" ) {
        vars.insert( var );
    }
    else if ( name == "filename" ) {
        filename = var;
    }
    else if ( name == "tree" ) {
        tree = var;
    }
    else if ( name == "ctmc" ) {
        ctmc = var;
    }
    else if ( name == "separator" ) {
        separator = var;
    }
    else if ( name == "printgen" ) {
        printgen = var;
    }
    else if ( name == "prior" ) {
        prior = var;
    }
    else if ( name == "posterior" ) {
        posterior = var;
    }
    else if ( name == "likelihood" ) {
        likelihood = var;
    }
    else if ( name == "type" ) {
        type = var;
    }
    else if ( name == "style" ) {
        style = var;
    }
    else if ( name == "append" ) {
        append = var;
    }
    else {
        RevObject::setConstParameter(name, var);
    }
}
