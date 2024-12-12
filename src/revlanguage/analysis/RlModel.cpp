
#include <cstddef>
#include <vector>
#include <algorithm>
#include <string>
#include <set>
#include <sstream>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Ellipsis.h"
#include "Model.h"
#include "RevObject.h"
#include "RlModel.h"
#include "RlString.h"
#include "TypeSpec.h"
#include "Argument.h"
#include "DagNode.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "RbBoolean.h"
#include "RbFileManager.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlUtils.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "Workspace.h"
#include "WorkspaceToCoreWrapperObject.h"

using namespace RevLanguage;

using std::vector;
using std::set;
using std::string;

Model::Model() : WorkspaceToCoreWrapperObject<RevBayesCore::Model>()
{
 
    auto* dotArgRules = new ArgumentRules();
    dotArgRules->push_back( new ArgumentRule("file", RlString::getClassTypeSpec(), "The name of the file where to save the model graph.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    dotArgRules->push_back( new ArgumentRule("verbose", RlBoolean::getClassTypeSpec(), "Verbose output?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    dotArgRules->push_back( new ArgumentRule("bg", RlString::getClassTypeSpec(), "The background color.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("lavenderblush2") ) );
    methods.addFunction( new MemberProcedure("graph", RlUtils::Void, dotArgRules) );

    auto* ignoreDataRules = new ArgumentRules();
    ignoreDataRules->push_back( new Ellipsis( "Clamped variables to ignore.", RevObject::getClassTypeSpec() ) );
    methods.addFunction( new MemberProcedure("ignoreData", RlUtils::Void, ignoreDataRules) );

    auto ignoreAllDataRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure("ignoreAllData", RlUtils::Void, ignoreAllDataRules) );
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Model* Model::clone(void) const
{
    
	return new Model(*this);
}


void Model::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // update vector variables
    Workspace::userWorkspace().updateVectorVariables();
    
    // now allocate a model
    std::set<const RevBayesCore::DagNode*> s;
    for (std::set<RevPtr<const RevVariable> >::iterator it = sources.begin(); it != sources.end(); ++it)
    {
        RevBayesCore::DagNode* n = (*it)->getRevObject().getDagNode();
        s.insert( n );
    }
    value = new RevBayesCore::Model( s );
//    printModelDotGraph();
}

vector<RevPtr<const RevVariable>> getElementVariables(RevPtr<const RevVariable> var)
{
    if (not var->isVectorVariable())
        return {var};

    vector<RevPtr<const RevVariable>> elems;
    for(int i=0;i<var->getMaxElementIndex();i++)
    {
        auto elem = var->getElementVariable(i);
        auto sub_elems = getElementVariables(elem);
        elems.insert(elems.end(), sub_elems.begin(), sub_elems.end());
    }

    return elems;
}

/* Map calls to member methods */
RevPtr<RevVariable> Model::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "graph")
    {
        found = true;
        
        const std::string&   fn      = static_cast<const RlString &>( args[0].getVariable()->getRevObject() ).getValue();
        bool vb = static_cast< const RlBoolean &>( args[1].getVariable()->getRevObject() ).getValue();
        const std::string&   bg      = static_cast<const RlString &>( args[2].getVariable()->getRevObject() ).getValue();
        printModelDotGraph(fn, vb, bg);
        
        return nullptr;
    }
    else if (name == "ignoreData")
    {
        found = true;

        set<string> names;
        for(auto& arg: args)
        {
            for(auto& var: getElementVariables( arg.getVariable() ) )
                names.insert( var->getName() );
        }

        ignoreDataAtNodes( names );

        return nullptr;
    }
    else if (name == "ignoreAllData")
    {
        found = true;

        ignoreAllData();

        return nullptr;
    }
    
    return RevObject::executeMethod( name, args, found );
}

/** Get Rev type of object */
const std::string& Model::getClassType(void)
{
    
    static std::string rev_type = "Model";
    
	return rev_type; 
}

/** Get class type spec describing type of object */
const TypeSpec& Model::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::Model>::getClassTypeSpec() ) );
    
	return rev_type_spec; 
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string Model::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "model";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Model::getParameterRules(void) const
{
    
    static MemberRules modelMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        modelMemberRules.push_back( new ArgumentRule("x", RevObject::getClassTypeSpec(), "Any variable that is connected in the model graph.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        modelMemberRules.push_back( new Ellipsis( "Additional variables.", RevObject::getClassTypeSpec() ) );
        
        rules_set = true;
    }
    
    return modelMemberRules;
}


/** Get type spec */
const TypeSpec& Model::getTypeSpec( void ) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Print a simplified representation of the model for the user. */
void Model::printValue(std::ostream &o, bool user) const
{
    
    const std::vector<RevBayesCore::DagNode*>& theNodes = value->getDagNodes();
    std::vector<RevBayesCore::DagNode*>::const_iterator it;

    o << std::endl;
    std::stringstream s;
    
    // compute the number of nodes by only counting nodes that are not hidden
    size_t num_nodes = 0;
    for ( it=theNodes.begin(); it!=theNodes.end(); ++it )
    {
    
        if ( (*it)->isHidden() == false )
        {
            ++num_nodes;
        }
    
    }
    
    s << "Model with " << num_nodes << " nodes";
    o << s.str() << std::endl;
    for ( size_t i = 0; i < s.str().size(); ++i )
        o << "=";
    o << std::endl << std::endl;
    
    for ( it=theNodes.begin(); it!=theNodes.end(); ++it )
    {
        RevBayesCore::DagNode *the_node = *it;
        // skip hidden nodes
        if ( the_node->isHidden() == true )
        {
            continue;
        }
        
        if ( the_node->getName() != "" )
        {
            o << the_node->getName() <<  " :" << std::endl;
        }
        else
        {
            o << "<" << the_node << "> :" << std::endl;
        }
        
        o << "_value        = ";
        std::ostringstream o1;
        the_node->printValue( o1, ", ", true );
        o << StringUtilities::oneLiner( o1.str(), 54 ) << std::endl;

        the_node->printStructureInfo( o, false );

        o << std::endl;
    }
}


/** Set a member variable */
void Model::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "" || name == "x")
    {
        sources.insert( var );
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}


/* Write a file in DOT format for viewing the model DAG in graphviz */
//   This requires the user to have graphviz installed, or they can paste the file contents
//   into http://graphviz-dev.appspot.com/
void Model::printModelDotGraph(const RevBayesCore::path &fn, bool vb, const std::string &bgc)
{
    
    const std::vector<RevBayesCore::DagNode*>& theNodes = value->getDagNodes();
    std::vector<RevBayesCore::DagNode*>::const_iterator it;
    
    RevBayesCore::createDirectoryForFile( fn );
    
    std::ofstream o( fn.string() );

    o << "/* Graphical model description in DOT language                                    */\n";
    o << "/*    To view graph:                                                              */\n";
    o << "/*       open this file in the program Graphviz: http://www.graphviz.org          */\n";
    o << "/*       or paste contents into an online viewer: http://stamm-wilbrandt.de/GraphvizFiddle */\n\n";
	o << "digraph REVDAG {\n";
    std::string nrank = "   {rank=same";
    for ( it=theNodes.begin(); it!=theNodes.end(); ++it ){
        if ( !(*it)->isHidden() || vb){
            std::stringstream nname;
            if ( (*it)->getName() != "" )
                nname << (*it)->getName();
            else
                nname << (*it);
            std::string stname = nname.str();
            std::replace( stname.begin(), stname.end(), '[', '_');
			std::replace( stname.begin(), stname.end(), '.', '_');

            stname.erase(std::remove(stname.begin(), stname.end(), ']'), stname.end());  
                  
            std::stringstream rl;
			if ((*it)->getName() == "" && !vb){
				rl << "function" ;
			}
			else
				rl << nname.str() ;
           
            // only print values of constant nodes (only simple numeric values)
            if ( (*it)->getDagNodeType() == RevBayesCore::DagNode::CONSTANT )
            {
                std::stringstream trl;
                if ((*it)->isSimpleNumeric())  
                    (*it)->printValue(trl," ", true);
                else 
                    trl << " ... ";
                if (trl.str() != "" || vb){
                    std::string val = trl.str();
                    if (val == "")
                        val = "constant\\ndefault member";
                    if (!(*it)->isSimpleNumeric())
                        val = "...";
                    o << "   n_" << stname;
                    o << " [shape=";
                    if ( (*it)->getName() == "" )
                    {
                        o << "box, style=filled, fillcolor=white, ";
                        o << "label=\"" << val << "\"]\n";
                    }
                    else
                    {
                        o << "record, style=filled, fillcolor=white, ";
                        rl << "|" << val;
                        o << "label=\"{" << rl.str() << "}\"]\n";
                    }
                }
            }
            else if ((*it)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC){
                o << "   n_" << stname;
                o << " [shape=";
                std::stringstream strss;
                (*it)->printStructureInfo(strss);
                if (strss.str().find("function",0) < strss.str().npos)
                {
                    std::string w;
                    
                    while(strss >> w){
                        if (w == "_function")
                        {
                            strss >> w;
                            strss >> w;
//                            std::cout << w << std::endl;
//                            strss >> w;
                            rl << "\\n[ " << w << ") ]";
                        }
                    }
                }
                else{
                    std::string w;
                    while(strss >> w){
                        if (w == "_dagType")
                        {
                            strss >> w;
                            strss >> w;
                            rl << "\\n[ " << w;
                            strss >> w;
                            if (w != "DAG")
                                rl << " " << w;
                            rl << " ]";
                        }
                    }
                }
                o << "oval, style=\"dashed,filled\", fillcolor=white, label=\"" << rl.str() << "\"]\n";
            }
            else if ((*it)->getDagNodeType() == RevBayesCore::DagNode::STOCHASTIC){
                o << "   n_" << stname;
                o << " [shape=";
                o << "oval, ";
                if ((*it)->isIgnoredData()){
                    o << "style=filled, fillcolor=white, color=gray, fontcolor=gray,";
                    nrank += "; n_" + stname;
                }
                else if ((*it)->isClamped()){
                    o << "style=filled, fillcolor=gray, ";
                    nrank += "; n_" + stname;
                }
                else 
                    o << "style=filled, fillcolor=white, ";
                o << "label=\"" << rl.str() << "\"]\n";
            }
        }
    }
    for ( it=theNodes.begin(); it!=theNodes.end(); ++it )
    {
        if ( !(*it)->isHidden() || vb)
        {
            std::stringstream trl;
            (*it)->printValue(trl,",", true);
            if (trl.str() != "" || vb)
            {
                std::stringstream nname;
                if ( (*it)->getName() != "" )
                {
                    nname << (*it)->getName() ;
                }
                else
                {
                    nname << (*it);
                }
                std::string stname = nname.str();
                std::replace( stname.begin(), stname.end(), '[', '_');
				std::replace( stname.begin(), stname.end(), '.', '_');
                stname.erase(std::remove(stname.begin(), stname.end(), ']'), stname.end());  

                if ((*it)->getNumberOfChildren() > 0)
                {
                    const std::vector<RevBayesCore::DagNode*>& childeren = (*it)->getChildren();
                    std::vector<RevBayesCore::DagNode*>::const_iterator ci;
                    for ( ci=childeren.begin(); ci!=childeren.end(); ++ci )
                    {
                        if ( (*ci)->isHidden() && vb==false )
                        {
                            RevBayesCore::DagNode *ch = (*ci)->getFirstChild();
                            while (ch->isHidden())
                            {
                                ch = ch->getFirstChild();
                            }
                            std::stringstream cn;
                            if ( ch->getName() != "" )
                                cn << ch->getName();
                            else
                                cn << ch;
                            std::string stcn = cn.str();
                            std::replace( stcn.begin(), stcn.end(), '[', '_');
                            stcn.erase(std::remove(stcn.begin(), stcn.end(), ']'), stcn.end());  
                            
                            o << "   n_" << stname << " -> n_";
                            o << stcn;
                            if (ch->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC)
                                o << "[style=dashed]";
                            o << "\n";
                        }
                        else
                        {
                            std::stringstream cn;
                            if ( (*ci)->getName() != "" )
                            {
                                cn << (*ci)->getName();
                            }
                            else
                            {
                                cn << (*ci);
                            }
                            std::string stcn = cn.str();
                            std::replace( stcn.begin(), stcn.end(), '[', '_');
                            stcn.erase(std::remove(stcn.begin(), stcn.end(), ']'), stcn.end());  

                            o << "   n_" << stname << " -> n_";
                            o << stcn;
                            if ((*ci)->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC)
                                o << "[style=dashed]";
                            o << "\n";
                        }
                        
                    }
                
                }
                
            }
            
        }
        
    }
    if (nrank.size() > 13)
    {
        nrank += ";}\n";
        o << nrank;
    }
    o << "   graph [bgcolor=" << bgc << ", pad=0.25]\n";
    o << "}";
    o.close();
}

void Model::ignoreDataAtNodes(const set<string>& namesToIgnore)
{
    auto& graphNodes = value->getDagNodes();

    std::map<string,RevBayesCore::DagNode*> nodeForName;
    for(auto& node: graphNodes)
        nodeForName.insert({node->getName(), node});

    for(auto& name: namesToIgnore)
    {
        if (not nodeForName.count(name))
            RBOUT("Warning: can't ignore data at node '" + name + "' because it is not in the model.");
        else
            nodeForName.at(name)->setIgnoreData(true);
    }
}

void Model::ignoreAllData()
{
    auto& graphNodes = value->getDagNodes();

    for(auto& node: graphNodes)
    {
        if (node->isClamped())
            node->setIgnoreData(true);
    }
}
