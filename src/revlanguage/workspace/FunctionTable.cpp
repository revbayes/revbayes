#include <cstddef>
#include <sstream>
#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "ArgumentRule.h"
#include "FunctionTable.h"
#include "RbException.h"
#include "RlFunction.h"
#include "Argument.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "StringUtilities.h"
#include "TypeSpec.h"

using namespace RevLanguage;

/** Basic constructor, empty table with or without parent */
FunctionTable::FunctionTable(FunctionTable* parent) : std::multimap<std::string, Function*>(),
    parentTable(parent)
{

}


/** Copy constructor */
FunctionTable::FunctionTable(const FunctionTable& x)
{
    
    for (std::multimap<std::string, Function *>::const_iterator it=x.begin(); it!=x.end(); ++it)
    {
        insert(std::pair<std::string, Function *>( it->first, ( it->second->clone() )));
    }
    
    parentTable = x.parentTable;
}


/** Destructor. We own the functions so we need to delete them. */
FunctionTable::~FunctionTable(void)
{

    clear();
    
}

/** Assignment operator */
FunctionTable& FunctionTable::operator=(const FunctionTable& x)
{

    if (this != &x)
    {

        clear();
        for (std::multimap<std::string, Function *>::const_iterator i=x.begin(); i!=x.end(); i++)
        {
            insert(std::pair<std::string, Function *>((*i).first, ( (*i).second->clone() ) ) );
        }
        
        parentTable = x.parentTable;
    }

    return (*this);
}


/**
 * Add function to table. We do various tests to ensure that the
 * function does not violate consistency rules, and we throw an
 * error if it does not have a distinct formal so that it cannot
 * overload existing functions, if its name is not unique.
 *
 * Note that we do not check parent frames, so the function can
 * hide (override if you wish) parent functions.
 *
 */
bool FunctionTable::addFunction( Function *func, bool avoid_except )
{
    std::string name = "";
    
    if ( func->isInternal() == true )
    {
        name = "_";
    }

    name += func->getFunctionName();

    // Test function compliance with basic rules
    testFunctionValidity( name, func );

    std::pair<std::multimap<std::string, Function *>::iterator,
              std::multimap<std::string, Function *>::iterator> ret_val;

    ret_val = equal_range(name);
    for (std::multimap<std::string, Function *>::iterator i=ret_val.first; i!=ret_val.second; i++)
    {
        if ( isDistinctFormal(i->second->getArgumentRules(), func->getArgumentRules()) == false )
        {
	    // Don't throw an exception if we are just going to ignore it.
	    if (avoid_except) return false;

            std::ostringstream msg;
            i->second->printValue(msg, true);
            msg << " cannot overload " << name << " = ";
            func->printValue(msg, true);
            msg << " : signatures are identical" << std::endl;
            
            // free memory
            delete func;
            
            // throw the error message
            throw RbException(msg.str());
        }
    }

    // Insert the function
    insert({name, func});
    
    std::vector<std::string> aliases = func->getFunctionNameAliases();
    for (size_t i=0; i < aliases.size(); ++i)
    {
        std::string a = aliases[i];
        // Insert the function
        insert({a, func->clone()});
    }
    return true;
}


/**
 * Clear table. We own the functions so we need
 * to delete them. When that is completed, we
 * call the base class clear function.
 */
void FunctionTable::clear(void)
{
    
    for ( std::multimap<std::string, Function *>::const_iterator i = begin(); i != end(); i++ )
    {
        Function *f = i->second;
        delete( f );
    }
    
    std::multimap<std::string, Function*>::clear();
    
}


/** Return a type-safe clone of the function table */
FunctionTable* FunctionTable::clone( void ) const
{
    
    return new FunctionTable( *this );
}


/** Erase function. @todo This does not work if there are several functions with the same name. Also memory leak. */
void FunctionTable::eraseFunction(const std::string& name)
{

    std::pair<std::multimap<std::string, Function *>::iterator,
              std::multimap<std::string, Function *>::iterator> ret_val;

    ret_val = equal_range(name);
    
    std::multimap<std::string, Function *>::iterator start = ret_val.first;
    std::multimap<std::string, Function *>::iterator end   = ret_val.second;
    for ( ; start != end; ++start )
    {
        Function *the_function = start->second;
        delete the_function;
    }
    
    erase(ret_val.first, ret_val.second);
    
}


/**
 * Does a function with the given name exist? We check this
 * function table and then delegate to parent table if we cannot
 * find the function here.
 */
bool FunctionTable::existsFunction(std::string const &name) const
{
    
    const std::map<std::string, Function *>::const_iterator& it = find( name );
    
    // if this table doesn't contain the function, then we ask the parent table
    if ( it == end() )
    {
        if ( parentTable != NULL )
        {
            return parentTable->existsFunction( name );
        }
        else
        {
            return false;
        }
    }
    
    return true;
}


/**
 * Does a function with the given name and signature exist in
 * this frame?
 */
bool FunctionTable::existsFunctionInFrame( std::string const &name, const ArgumentRules& r ) const
{
    std::map<std::string, Function *>::const_iterator it = find( name );
    
    // If the name does not exist, the answer is no
    if ( it == end() )
        return false;

    // The name exists, so we cycle through the functions and check whether the signature is distinct
    std::pair<std::multimap<std::string, Function *>::const_iterator,
              std::multimap<std::string, Function *>::const_iterator> range;
    range = equal_range( name );

    for ( it = range.first; it != range.second; ++it )
    {
        if ( !isDistinctFormal( it->second->getArgumentRules(), r ) )
            return true;
    }
    
    return false;
}


/**
 * Find functions matching name
 *
 * @todo Inherited functions are not returned if there
 *       are functions matching the name in the current
 *       workspace.
 */
std::vector<Function *> FunctionTable::findFunctions(const std::string& name) const
{

    std::vector<Function *>  the_functions;

    size_t hits = count(name);
    if (hits == 0)
    {
        if (parentTable != NULL)
        {
            return parentTable->findFunctions( name );
        }
        else
        {
            return the_functions;
        }
        
    }

    std::pair<std::multimap<std::string, Function *>::const_iterator,
              std::multimap<std::string, Function *>::const_iterator> ret_val;
    ret_val = equal_range( name );

    std::multimap<std::string, Function *>::const_iterator it;
    for ( it=ret_val.first; it!=ret_val.second; it++ )
    {
        the_functions.push_back( (*it).second->clone() );
    }
    
    return the_functions;
}


/** Find function (also processes arguments) */
const Function* FunctionTable::findFunction(const std::string& name, const std::vector<Argument>& args) const
{
    
    std::pair<std::multimap<std::string, Function *>::const_iterator,
              std::multimap<std::string, Function *>::const_iterator> ret_val;
    
    size_t hits = count(name);
    if (hits == 0)
    {
        
        if (parentTable != NULL)
        {
            // \TODO: We shouldn't allow const casts!!!
            FunctionTable* pt = const_cast<FunctionTable*>(parentTable);
            return pt->findFunction(name, args);
        }
        else
        {
            // Here we just return a nullptr to indicate failure.
            // The caller can produce better error messages, so it is the caller's job to throw an exception.
            return nullptr;
        }
        
    }
    ret_val = equal_range(name);
    if (hits == 1)
    {
        std::vector<bool> arg_mapped(args.size(), false);
        if (ret_val.first->second->checkArguments(args,NULL,arg_mapped) == false)
        {
            std::ostringstream msg;

            // get whitespace offset from function name (+2 for " (")
            std::string whitespace(name.size() + 2, ' ');
            
            msg << "Argument or label mismatch for function call.\n";
            msg << "Arguments which could not be matched (we stopped after the first mismatch):\n";
            for (size_t i=0; i<arg_mapped.size(); ++i)
            {
                if ( arg_mapped[i] == false )
                {
                    msg << "\t" << args[i].getLabel() << "\n";
                    break;
                }
            }
            msg << "Provided call:\n";
            msg << name << " (";

            // print the passed arguments
            for (std::vector<Argument>::const_iterator it = args.begin(); it != args.end(); it++)
            {
                // add a comma and whitespace before the every argument except the first
                if (it != args.begin())
                {
                    msg << ",\n" << whitespace;
                }
                
                // create the default type of the passed-in argument
                std::string type = "NULL";
                // get the type if the variable wasn't NULL
                if (it->getVariable() != NULL)
                {
                    type = it->getVariable()->getRevObject().getType();
                }
                msg << type;
                
                // create the default DAG type of the passed-in argument
                std::string dag_type = "";
                // get the type if the variable wasn't NULL
                if (it->getVariable() != NULL && it->getVariable()->getRevObject().isModelObject() == true && it->getVariable()->getRevObject().getDagNode() != NULL )
                {
                    if ( it->getVariable()->getRevObject().getDagNode()->getDagNodeType() == RevBayesCore::DagNode::DETERMINISTIC )
                    {
                        dag_type = "<deterministic>";
                    }
                    else if ( it->getVariable()->getRevObject().getDagNode()->getDagNodeType() == RevBayesCore::DagNode::STOCHASTIC )
                    {
                        dag_type = "<stochastic>";
                    }
                    else if ( it->getVariable()->getRevObject().getDagNode()->getDagNodeType() == RevBayesCore::DagNode::CONSTANT )
                    {
                        dag_type = "<constant>";
                    }
                    else
                    {
                        dag_type = "<?>";
                    }
                }
                msg << dag_type;
                
                if ( it->getLabel() != "" )
                {
                    msg << " '" << it->getLabel() << "'";
                }
            }
            msg << " )\n" << std::endl;
            msg << "Correct usage is:" << std::endl;
            ret_val.first->second->printValue( msg, true );
            msg << std::endl;
            throw RbException( msg.str() );
        }
        return ret_val.first->second;
    }
    else
    {
        std::vector<double>* match_score = new std::vector<double>();
        std::vector<double> best_score;
        Function* best_match = NULL;

        bool ambiguous = false;
        std::multimap<std::string, Function *>::const_iterator it;
        for (it=ret_val.first; it!=ret_val.second; it++)
        {
            match_score->clear();
            std::vector<bool> arg_mapped(args.size(), false);
            if ( (*it).second->checkArguments(args, match_score, arg_mapped) == true )
            {
                std::sort(match_score->begin(), match_score->end(), std::greater<double>());
                if ( best_match == NULL )
                {
                    best_score = *match_score;
                    best_match = it->second;
                    ambiguous = false;
                }
                else
                {
                    size_t j;
                    for (j=0; j<match_score->size() && j<best_score.size(); ++j)
                    {
                        
                        if ( (*match_score)[j] < best_score[j] )
                        {
                            best_score = *match_score;
                            best_match = it->second;
                            ambiguous = false;
                            break;
                        }
                        else if ((*match_score)[j] > best_score[j])
                        {
                            break;
                        }
                        
                    }
                    if (j==match_score->size() || j==best_score.size())
                    {
                        ambiguous = true;   // Continue checking, there might be better matches ahead
                    }
                    
                }
                
            }
            
        }
        
        // free the memory
        delete match_score;
        
        /* Delete all processed arguments except those of the best matching function, if it is ambiguous */
        for ( it = ret_val.first; it != ret_val.second; it++ )
        {
            if ( !( (*it).second == best_match && ambiguous == false ) )
            {
                (*it).second->clear();
            }
            
        }
        if ( best_match == NULL || ambiguous == true )
        {
            std::ostringstream msg;
            if ( best_match == NULL )
            {
                msg << "No overloaded function '" << name << "' matches for arguments (";
            }
            else
            {
                msg << "Ambiguous call to function '" << name << "' with arguments (";
            }
            // print the passed arguments
            for (std::vector<Argument>::const_iterator j = args.begin(); j != args.end(); j++)
            {
                if (j != args.begin())
                {
                    msg << ",";
                }
                const RevPtr<const RevVariable>& the_var = j->getVariable();
                msg << " " << the_var->getRevObject().getTypeSpec().getType();
                
            }
            msg << " )" << std::endl;
            
            msg << "Potentially matching functions are:" << std::endl;
            for ( it = ret_val.first; it != ret_val.second; it++ )
            {
                (*it).second->printValue( msg, true );
                msg << std::endl;
            }
            throw RbException( msg.str() );
        }
        else
        {
            return best_match;
        }
        
    }
    
}


/**
 * Get first function. This function will find the first function with a matching name without
 * throwing an error. Compare with the getFunction(name) function, which will throw an error
 * if the function name is overloaded.
 */
Function* FunctionTable::getFirstFunction( const std::string& name ) const
{
    // find the template function
    std::vector<Function *> the_functions = findFunctions(name);
    
    if ( the_functions.size() == 0 )
    {
        throw RbException() << "Could not find function with name '" << name << "'";
    }
    
    // free memory
    for (size_t i=1; i<the_functions.size(); ++i)
    {
        Function *the_function = the_functions[i];
        delete the_function;
        
        // just for savety
        the_functions[i] = NULL;
    }
    
    return the_functions[0];
}


/** Get function. This function will throw an error if the name is missing or if there are several matches (overloaded functions) */
Function* FunctionTable::getFunction( const std::string& name ) const
{
    
    // find the template function
    std::vector<Function *> the_functions = findFunctions(name);
    
    // free memory
    for (size_t i=1; i<the_functions.size(); ++i)
    {
        Function *the_function = the_functions[i];
        delete the_function;
        
        // just for savety
        the_functions[i] = NULL;
    }
    
    if ( the_functions.size() > 1 )
    {
        Function *the_function = the_functions[0];
        delete the_function;
        
        std::ostringstream o;
        o << "Found " << the_functions.size() << " functions with name \"" << name + "\". Identification not possible if arguments are not specified.";
        throw RbException( o.str() );
    }
    
    return the_functions[0];
}


/** Get function. This function will throw an error if the name and args do not match any named function. */
const Function& FunctionTable::getFunction(const std::string& name, const std::vector<Argument>& args) const
{
    
    // find the template function
    const Function* the_function = findFunction(name, args);

    if ( the_function == NULL )
    {
        throw RbException() << "No function named '" << name << "'"; 
    }
    
    return *the_function;
}

void FunctionTable::getFunctionNames(std::vector<std::string>& names) const
{
    for (std::multimap<std::string, Function *>::const_iterator i=begin(); i!=end(); i++)
    {
        std::string s = i->second->getFunctionName();
        names.push_back(s);
    }
    
    if ( parentTable != NULL)
    {
        parentTable->getFunctionNames(names);
    }
}


/** Check if two formals are unique */
bool FunctionTable::isDistinctFormal(const ArgumentRules& x, const ArgumentRules& y) const
{

    /* Check that all labels are unique in both sets of argument rules */
    for (size_t i=0; i<x.size(); i++)
    {
        for (size_t j=i+1; j < x.size(); j++)
        {
            if (x[i].getArgumentLabel().size() != 0 && x[j].getArgumentLabel().size() != 0)
            {
             
                if (x[i].getArgumentLabel() == x[j].getArgumentLabel())
                {
                    return false;
                }
                
            }
            
        }
        
    }
    for (size_t i=0; i<y.size(); i++)
    {
        for (size_t j=i+1; j<y.size(); j++)
        {
            
            if (y[i].getArgumentLabel().size() != 0 && y[j].getArgumentLabel().size() != 0)
            {
                if (y[i].getArgumentLabel() == y[j].getArgumentLabel())
                {
                    return false;
                }
                
            }
            
        }
        
    }

    /* Check that types are different for at least one argument without default values */
    size_t i;
    for (i=0; i<x.size() && i<y.size(); i++)
    {
        if ( !(x[i].hasDefault() == true && y[i].hasDefault() == true) &&
            !x[i].isEllipsis() && !y[i].isEllipsis() &&
            (x[i].getArgumentTypeSpec() != y[i].getArgumentTypeSpec() || x[i].getArgumentDagNodeType() != y[i].getArgumentDagNodeType() ))
        {
            return true;
        }
        
    }
    for (size_t j=i; j<x.size(); j++)
    {
        if (x[j].hasDefault() == false && !x[j].isEllipsis())
        {
            return true;
        }
        
    }
    for (size_t j=i; j<y.size(); j++)
    {
        if (y[j].hasDefault() == false && !y[j].isEllipsis())
        {
            return true;
        }
        
    }

    return false;
}


/**
 * Is the function with the given name a procedure? We check this
 * function table and then delegate to parent table if we cannot
 * find the function here.
 */
bool FunctionTable::isProcedure(const std::string& name) const
{
    const std::map<std::string, Function *>::const_iterator& it = find( name );
    
    // If we have the function, we know the answer
    if ( it != end() )
    {
        return it->second->isProcedure();
    }
    
    // If this table doesn't contain the function, then we ask the parent table
    if ( parentTable != NULL )
    {
        return parentTable->isProcedure( name );
    }
    else
    {
        throw RbException() << "No function or procedure '" << name << "'" ; 
    }
    
}


/** Print function table for user in pretty format */
void FunctionTable::printValue(std::ostream& o, bool env) const
{
    
    for (std::multimap<std::string, Function *>::const_iterator i=begin(); i!=end(); i++)
    {
        std::ostringstream s("");

        s << i->first << " = ";
        i->second->printValue( s, true );
        
        // nice printing with arguments aligned
        size_t offset = i->first.size() + 3;
        std::string temp = "";
        std::string repl = "";
        
        // remove the whitespace *before* the left parenthesis (look backward, while constructing the string)
        for (size_t i = 0; i < s.str().size(); ++i)
        {
            temp += s.str()[i];
            
            if (s.str()[i] == '(')
            {
                temp[i - 1] = *repl.c_str();
            }
        }
        
        // remove the newline and whitespace *after* a comma (look forward, after the string has been constructed)
        for (size_t i = 0; i < s.str().size(); ++i)
        {
            if (s.str()[i] == ',')
            {
                for (size_t j = 1; j < offset; ++j)
                {
                    temp[i + j] = *repl.c_str();
                }
            }
        }
        
        o << temp << std::endl;
    }
    
    // Print the parent table too
    if ( parentTable != NULL && env == true )
    {
        parentTable->printValue(o , env );
    }
}


/**
 * Replace function. Unlike the addFunction function, we do not throw an error if the
 * function exists. Instead we replace the existing function.
 *
 * Note that we make the same tests as in the addFunction() function to make sure that
 * the consistency of the function table is maintained.
 */
void FunctionTable::replaceFunction( const std::string& name, Function *func )
{
    // Test the function
    testFunctionValidity( name, func );
    
    // Find the function to be replaced
    std::pair<std::multimap<std::string, Function *>::iterator,
              std::multimap<std::string, Function *>::iterator> range;

    range = equal_range( name );
    for ( std::multimap<std::string, Function *>::iterator it = range.first; it != range.second; ++it )
    {
        if ( !isDistinctFormal( it->second->getArgumentRules(), func->getArgumentRules() ) )
        {
            delete it->second;
            it->second = func;
            return;
        }
    }
    
    // No match; simply insert the function
    insert(std::pair<std::string, Function* >( name, func ) );
    
    // Name the function so that it is aware of what it is called
    func->setName( name );
}


/**
 * Test whether the function can be added to the table. To be added, it needs to follow
 * these rules:
 *
 *    1. A function cannot overload a procedure, and vice versa
 *    2. A function or procedure overloading an existing function or procedure
 *       must have the same return type
 *
 * @todo We currently have a number of functions (from templated C++ code) with different
 *       return types registered as overloaded functions in RbRegister. Fix this.
 */
void FunctionTable::testFunctionValidity( const std::string& name, Function* func ) const
{
    // We only need to make these tests if the function name already exists
    if ( existsFunction( name ) )
    {
        Function* fxn = getFirstFunction( name );
        
        // Functions need to be of same type (procedure or function)
        if ( fxn->isProcedure() != func->isProcedure() )
        {
            // Construct an error message
            std::ostringstream msg;
            if ( func->isProcedure() )
            {
                msg << "Procedure ";
            }
            else
            {
                msg << "Function ";
            }
            
            msg << name << " =  ";
            func->printValue(msg, true);
            
            msg << " cannot overload ";
            if ( fxn->isProcedure() )
            {
                msg << " procedure ";
            }
            else
            {
                msg << " function ";
            }
            msg << name << " = ";
            fxn->printValue(msg, true);
            msg << " : procedure/function mismatch" << std::endl;
            
            // free function memory
            delete fxn;
            
            // throw the error message
            throw RbException(msg.str());
        }
        
        // Functions must have same return type
#if 0
        if ( fxn.getReturnType() != func->getReturnType() )
        {
            // Construct an error message
            std::ostringstream msg;
            if ( func->isProcedure() )
                msg << "Procedure ";
            else
                msg << "Function ";
            
            msg << name << " =  ";
            func->printValue(msg);
            
            msg << " cannot overload ";
            if ( fxn.isProcedure() )
                msg << " procedure ";
            else
                msg << " function ";
            msg << name << " = ";
            fxn.printValue(msg);
            msg << " : return types differ" << std::endl;
            
            // free function memory
            delete fxn;
            
            // throw the error message
            throw RbException(msg.str());
        }
#endif
        
        
        // free function memory
        delete fxn;
    }
    
}
