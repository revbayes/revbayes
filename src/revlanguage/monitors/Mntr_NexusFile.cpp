#include "Mntr_NexusFile.h"

#include <cstddef>
#include <algorithm>
#include <ostream>
#include <string>

#include "Ellipsis.h"
#include "IntegerPos.h"
#include "NexusMonitor.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "TypeSpec.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

namespace RevLanguage {

Mntr_NexusFile::Mntr_NexusFile() : Monitor() {}

Mntr_NexusFile* Mntr_NexusFile::clone(void) const {
    return new Mntr_NexusFile(*this);
}


void Mntr_NexusFile::constructInternalObject( void ) {
    // we free the memory first
    delete value;

    // now allocate a new monitor
    const std::string& fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    unsigned int g = (int)static_cast<const IntegerPos &>( printgen->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *t = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();

    vars.erase( unique( vars.begin(), vars.end() ), vars.end() );
    sort( vars.begin(), vars.end(), compareVarNames );
    std::vector<RevBayesCore::DagNode *> n;
    for (std::vector<RevPtr<const RevVariable> >::iterator i = vars.begin(); i != vars.end(); ++i)
    {
        RevBayesCore::DagNode* node = (*i)->getRevObject().getDagNode();
        n.push_back( node );
    }

    bool np = static_cast<const RlBoolean &>( isNodeParameter->getRevObject() ).getValue();
    bool taxa = static_cast<const RlBoolean &>( writeTaxa->getRevObject() ).getValue();

    RevBayesCore::NexusMonitor* m = new RevBayesCore::NexusMonitor(t, n, np, size_t(g), fn, false, taxa);
    value = m;
}


const std::string& Mntr_NexusFile::getClassType(void) {

    static std::string rev_type = "Mntr_NexusFile";
    return rev_type;
}


const TypeSpec& Mntr_NexusFile::getClassTypeSpec(void) {

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    return rev_type_spec;
}


std::string Mntr_NexusFile::getMonitorName( void ) const {
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "Nexus";
    return c_name;
}


const MemberRules& Mntr_NexusFile::getParameterRules(void) const {

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set ) {

        memberRules.push_back( new ArgumentRule("filename", RlString::getClassTypeSpec(), "The name of the file.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("tree"    , TimeTree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new Ellipsis( "Variables at nodes or branches.", RevObject::getClassTypeSpec() ) );
        memberRules.push_back( new ArgumentRule("isNodeParameter" , RlBoolean::getClassTypeSpec(), "Is this a node or branch parameter?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("writeTaxa" , RlBoolean::getClassTypeSpec(), "Should a taxa block be written?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        memberRules.push_back( new ArgumentRule("printgen"  , IntegerPos::getClassTypeSpec()  , "The number of generations between stored samples.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new IntegerPos(1) ) );

        rules_set = true;
    }

    return memberRules;
}


const TypeSpec& Mntr_NexusFile::getTypeSpec( void ) const {

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}

void Mntr_NexusFile::printValue(std::ostream &o) const {

    o << "Mntr_NexusFile";
}


void Mntr_NexusFile::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {

    if ( name == "" )
    {
        vars.push_back(var);
    }
    else if ( name == "filename" )
    {
        filename = var;
    }
    else if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "isNodeParameter" )
    {
        isNodeParameter = var;
    }
    else if ( name == "printgen" )
    {
        printgen = var;
    }
    else if ( name == "writeTaxa" )
    {
        writeTaxa = var;
    }
    else
    {
        Monitor::setConstParameter(name, var);
    }
}

} /* namespace RevLanguage */
