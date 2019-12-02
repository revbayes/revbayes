#include "ArgumentRule.h"
#include "ConstantNode.h"
#include "Ellipsis.h"
#include "Func_readTrees.h"
#include "ModelVector.h"
#include "NclReader.h"
#include "NewickConverter.h"
#include "OptionRule.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlBranchLengthTree.h"
#include "RlString.h"
#include "RlTimeTree.h"
#include "RlUtils.h"
#include "StringUtilities.h"
#include "Tree.h"
#include "RlUserInterface.h"

#include <map>
#include <set>
#include <sstream>


using namespace RevLanguage;

/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_readTrees* Func_readTrees::clone( void ) const
{

    return new Func_readTrees( *this );
}


/** Execute function */
RevPtr<RevVariable> Func_readTrees::execute( void )
{

    // get the information from the arguments for reading the file
    const std::string&  fn = static_cast<const RlString&>( args[0].getVariable()->getRevObject() ).getValue();
    const std::string&  text = static_cast<const RlString&>( args[1].getVariable()->getRevObject() ).getValue();
    const std::string&  treetype = static_cast<const RlString&>( args[2].getVariable()->getRevObject() ).getValue();


    if (fn != "")
    {
        // get the global instance of the NCL reader and clear warnings from its warnings buffer
        RevBayesCore::NclReader reader = RevBayesCore::NclReader();

        if ( treetype == "clock" )
        {
            ModelVector<TimeTree> *trees = new ModelVector<TimeTree>();
            
            std::vector<RevBayesCore::Tree*> tmp = reader.readTimeTrees( fn );
            for (std::vector<RevBayesCore::Tree*>::iterator t = tmp.begin(); t != tmp.end(); ++t)
            {
                trees->push_back( TimeTree(*t) );
            }
            return new RevVariable( trees );
        }
        else if ( treetype == "non-clock" )
        {
            ModelVector<BranchLengthTree> *trees = new ModelVector<BranchLengthTree>();
            
            std::vector<RevBayesCore::Tree*>* tmp = reader.readBranchLengthTrees( fn );
            for (std::vector<RevBayesCore::Tree*>::iterator t = tmp->begin(); t != tmp->end(); ++t)
            {
                trees->push_back( BranchLengthTree(*t) );
            }
            delete tmp;
            return new RevVariable( trees );
        }
        
    }

    if (text != "")
    {

        string aux;
        std::istringstream iss(text);
        RevBayesCore::NewickConverter c;
        
        
        RevBayesCore::RbFileManager file_reader = RevBayesCore::RbFileManager();
        if ( treetype == "clock" )
        {
            ModelVector<TimeTree> *trees = new ModelVector<TimeTree>();
            while (file_reader.safeGetline(iss, aux))
            {
                RevBayesCore::Tree *blTree = c.convertFromNewick( aux );
                trees->push_back( TimeTree(*blTree) );
                
                delete blTree;
                
            }
            return new RevVariable( trees );
        }
        else if ( treetype == "non-clock" )
        {
            ModelVector<BranchLengthTree> *trees = new ModelVector<BranchLengthTree>();
            while (file_reader.safeGetline(iss, aux))
            {
                RevBayesCore::Tree *blTree = c.convertFromNewick( aux );
                trees->push_back( BranchLengthTree(*blTree) );
                
                delete blTree;
                
            }
            return new RevVariable( trees );
            
        }

        
        
    }


    return NULL;

}


/** Get argument rules */
const ArgumentRules& Func_readTrees::getArgumentRules( void ) const
{

    static ArgumentRules argumentRules = ArgumentRules();
    static bool rules_set = false;

    if ( rules_set == false )
    {

        argumentRules.push_back( new ArgumentRule( "file", RlString::getClassTypeSpec(), "The name of the file containing the trees.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
        argumentRules.push_back( new ArgumentRule( "text", RlString::getClassTypeSpec(), "A string containing one or several newick trees, separated by end of lines.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("") ) );
        
        std::vector<std::string> tree_options;
        tree_options.push_back( "clock" );
        tree_options.push_back( "non-clock" );
        argumentRules.push_back( new OptionRule( "treetype", new RlString("clock"), tree_options, "The type of trees." ) );
        
        rules_set = true;
    }

    return argumentRules;
}


/** Get Rev type of object */
const std::string& Func_readTrees::getClassType(void)
{

    static std::string rev_type = "Func_readTrees";

	return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Func_readTrees::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );

	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_readTrees::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "readTrees";

    return f_name;
}


/** Get type spec */
const TypeSpec& Func_readTrees::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get return type */
const TypeSpec& Func_readTrees::getReturnType( void ) const
{

    static TypeSpec return_typeSpec = ModelVector<TimeTree>::getClassTypeSpec();
    return return_typeSpec;
}



/**
 * Get the author(s) of this function so they can receive credit (and blame) for it.
 */
std::vector<std::string> Func_readTrees::getHelpAuthor(void) const
{
    // create a vector of authors for this function
    std::vector<std::string> authors;
    authors.push_back( "Bastien Boussau" );

    return authors;
}


/**
 * Get the (brief) description for this function
 */
std::string Func_readTrees::getHelpDescription(void) const
{
    // create a variable for the description of the function
    std::string description = "Reads trees from a file containing trees (Nexus, Phylip or Newick accepted), or from a string containing Newick representations of trees.";

    return description;
}


/**
 * Get the more detailed description of the function
 */
std::string Func_readTrees::getHelpDetails(void) const
{
    // create a variable for the description of the function
    std::string details = "Either a file name (with the file argument) or a string (with the text argument) must be provided as argument. If both are provided, trees will be read from both sources.";

    return details;
}


/**
 * Get an executable and instructive example.
 * These example should help the users to show how this function works but
 * are also used to test if this function still works.
 */
std::string Func_readTrees::getHelpExample(void) const
{
    // create an example as a single string variable.
    std::string example = "";
    example += "trees=readTrees(text=\"(a,(b,c));\\n(d:0.1,(e:0.1,f:0.1):0.1);\")\n";
    example += "print(trees)\n";
    return example;
}


/**
 * Get some references/citations for this function
 *
 */
std::vector<RevBayesCore::RbHelpReference> Func_readTrees::getHelpReferences(void) const
{
    // create an entry for each reference
    std::vector<RevBayesCore::RbHelpReference> references;

    return references;
}


/**
 * Get the names of similar and suggested other functions
 */
std::vector<std::string> Func_readTrees::getHelpSeeAlso(void) const
{
    // create an entry for each suggested function
    std::vector<std::string> see_also;

    see_also.push_back("write");
    see_also.push_back("readBranchLengthTrees");
    see_also.push_back("readCharacterData");
    see_also.push_back("readCharacterDataDelimited");
    see_also.push_back("readContinuousCharacterData");
    see_also.push_back("readDiscreteCharacterData");
    see_also.push_back("readDataDelimitedFile");
    see_also.push_back("readCharacterData");

    return see_also;
}


/**
 * Get the title of this help entry
 */
std::string Func_readTrees::getHelpTitle(void) const
{
    // create a title variable
    std::string title = "Function to read in trees.";

    return title;
}
