
#include <cstddef>
#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Ellipsis.h"
#include "ExtendedNewickTreeMonitor.h"
#include "Mntr_ExtendedNewickFile.h"
#include "IntegerPos.h"
#include "RevObject.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "Monitor.h"
#include "RbBoolean.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlMonitor.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Tree; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }


using namespace RevLanguage;



Mntr_ExtendedNewickFile::Mntr_ExtendedNewickFile(void) : FileMonitor()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Mntr_ExtendedNewickFile* Mntr_ExtendedNewickFile::clone(void) const
{
    
	return new Mntr_ExtendedNewickFile(*this);
}


void Mntr_ExtendedNewickFile::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new sliding move
    const std::string& fn = static_cast<const RlString &>( filename->getRevObject() ).getValue();
    const std::string& sep = static_cast<const RlString &>( separator->getRevObject() ).getValue();
    unsigned int g = (int)static_cast<const IntegerPos &>( printgen->getRevObject() ).getValue();
    bool wv = static_cast<const RlBoolean &>( version->getRevObject() ).getValue();
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *t = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    
    vars.erase( unique( vars.begin(), vars.end() ), vars.end() );
    sort( vars.begin(), vars.end(), compareVarNames );
    std::vector<RevBayesCore::DagNode *> n;
    for (std::vector<RevPtr<const RevVariable> >::iterator i = vars.begin(); i != vars.end(); ++i)
    {
        RevBayesCore::DagNode* node = (*i)->getRevObject().getDagNode();
        n.push_back( node );
    }
//    
//    
//    std::set<RevBayesCore::DagNode*> n;
//    for (std::set<RevPtr<const RevVariable> >::iterator i = vars.begin(); i != vars.end(); ++i)
//    {
//        RevBayesCore::DagNode* node = (*i)->getRevObject().getDagNode();
//        n.insert( node );
//    }
    bool np = static_cast<const RlBoolean &>( isNodeParameter->getRevObject() ).getValue();
    bool pp = static_cast<const RlBoolean &>( posterior->getRevObject() ).getValue();
    bool l = static_cast<const RlBoolean &>( likelihood->getRevObject() ).getValue();
    bool pr = static_cast<const RlBoolean &>( prior->getRevObject() ).getValue();

    RevBayesCore::ExtendedNewickTreeMonitor* m = new RevBayesCore::ExtendedNewickTreeMonitor(t, n, np, size_t(g), fn, sep, pp, l, pr);

    m->setPrintVersion( wv );

    value = m;
}


/** Get Rev type of object */
const std::string& Mntr_ExtendedNewickFile::getClassType(void)
{
    
    static std::string rev_type = "Mntr_ExtendedNewickFile";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Mntr_ExtendedNewickFile::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Monitor::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Mntr_ExtendedNewickFile::getMonitorName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "ExtNewick";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Mntr_ExtendedNewickFile::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
    
        memberRules.push_back( new ArgumentRule("tree"    , TimeTree::getClassTypeSpec(), "The tree variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new Ellipsis( "Variables at nodes or branches.", RevObject::getClassTypeSpec() ) );
        memberRules.push_back( new ArgumentRule("isNodeParameter" , RlBoolean::getClassTypeSpec(), "Is this a node or branch parameter?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("posterior" , RlBoolean::getClassTypeSpec(), "Should we print the posterior probability as well.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("likelihood", RlBoolean::getClassTypeSpec(), "Should we print the likelihood as well?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        memberRules.push_back( new ArgumentRule("prior"     , RlBoolean::getClassTypeSpec(), "Should we print the prior probability as well?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
        
        // add the rules from the base class
        const MemberRules &parentRules = FileMonitor::getParameterRules();
        memberRules.insert(memberRules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }
    
    return memberRules;
}

/** Get type spec */
const TypeSpec& Mntr_ExtendedNewickFile::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void Mntr_ExtendedNewickFile::printValue(std::ostream &o) const
{
    
    o << "Mntr_ExtendedNewickFile";
}


/** Set a member variable */
void Mntr_ExtendedNewickFile::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) {
    
    if ( name == "" )
    {
        vars.push_back( var );
    }
    else if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "isNodeParameter" )
    {
        isNodeParameter = var;
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
    else
    {
        FileMonitor::setConstParameter(name, var);
    }
}
