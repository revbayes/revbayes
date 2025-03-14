#include "Environment.h"

#include <sstream> // IWYU pragma: keep
#include <utility>

#include "RbException.h"
#include "RbHelpSystem.h"
#include "RlFunction.h"
#include "RevVariable.h"
#include "RbHelpFunction.h"

namespace RevLanguage { class Argument; }
namespace RevLanguage { class RevObject; }

using namespace RevLanguage;


/** Construct environment with NULL parent */
Environment::Environment(const std::string &n) :
    function_table(),
    numUnnamedVariables(0),
    parentEnvironment(NULL),
    variableTable(),
    children(),
    name( n )
{
    
}


/** Construct environment with parent */
Environment::Environment(const std::shared_ptr<Environment>& parentEnv, const std::string &n) :
    function_table(&parentEnv->getFunctionTable()),
    numUnnamedVariables(0),
    parentEnvironment(parentEnv),
    variableTable(),
    children(),
    name( n )
{
    
}


/** Copy constructor. We need to make sure we have a deep copy of the variable table, because a shallow copy is the default. */
Environment::Environment(const Environment &x) :
    function_table( x.function_table ),
    numUnnamedVariables( x.numUnnamedVariables ),
    parentEnvironment( x.parentEnvironment ),
    variableTable( x.variableTable ),
    children(),
    name( x.name )
{
    // Make a deep copy of the variable table by making sure we have clones of the variables
    for ( VariableTable::iterator i = variableTable.begin(); i != variableTable.end(); i++ )
    {
        i->second = i->second->clone();
    }
    
}


/** Destructor */
Environment::~Environment()
{
    // Clear the variable table and function table
    clear();
    
    // Child environments will be destroyed when they are completely unreferenced.
    children.clear();
}


/** Assignment operator. Ensure we have a deep copy of the variable table. */
Environment& Environment::operator=(const Environment &x)
{
    if (this != &x)
    {
        // Clean up variable and function table
        clear();

        // Copy parent environment pointer
        parentEnvironment = x.parentEnvironment;

        // Make a deep copy of function table by using assignment operator in FunctionTable
        function_table = x.function_table;
    }

    return *this;
}


/** Add alias to variable to frame. */
void Environment::addAlias( const std::string& name, const RevPtr<RevVariable>& the_var )
{
    
    /* Throw an error if the name string is empty. */
    if ( name == "" )
    {
        throw RbException("Invalid attempt to add unnamed reference variable to frame.");
    }
    
    /* Throw an error if the variable exists. Note that we cannot use the function
     existsVariable because that function looks recursively in parent frames, which
     would make it impossible to hide global variables. */
    if ( variableTable.find( name ) != variableTable.end() )
    {
        throw RbException() << "Variable " << name << " already exists in frame" ; 
    }
    
    /* Insert new alias to variable in variable table (we do not and should not name it) */
    variableTable.insert( std::pair<std::string, RevPtr<RevVariable> >( name, the_var ) );
    
}


/* Add function to frame. */
bool Environment::addFunction( Function* func )
{

    if ( existsVariable( func->getFunctionName() ) )
    {
        // free memory
        delete func;
        
        throw RbException() << "There is already a variable named '" << name << "' in the workspace";
    }
    
    function_table.addFunction( func );
    
    // add the help entry for this function to the global help system instance
    // but only if this is not an internal function
    if ( func->isInternal() == false )
    {
        RevBayesCore::RbHelpSystem::getHelpSystem().addHelpFunction( static_cast<RevBayesCore::RbHelpFunction*>(func->getHelpEntry()) );
    }

    return true;
}


/** Add an empty (NULL) variable to frame. */
void Environment::addNullVariable( const std::string& name )
{
    // Create a new RevVariable object
    RevPtr<RevVariable> var = RevPtr<RevVariable>( new RevVariable( NULL ) );
    
    // Add the variable to the table
    addVariable( name, var );
}


/** Add reference to variable to frame. */
void Environment::addReference( const std::string& name, const RevPtr<RevVariable>& the_var )
{
    
    /* Throw an error if the name string is empty. */
    if ( name == "" )
    {
        throw RbException("Invalid attempt to add unnamed reference variable to frame.");
    }
    
    /* Throw an error if the variable exists. Note that we cannot use the function
     existsVariable because that function looks recursively in parent frames, which
     would make it impossible to hide global variables. */
    if ( variableTable.find( name ) != variableTable.end() )
    {
        throw RbException() << "Variable " << name << " already exists in frame" ; 
    }
    
    /* Insert new reference variable in variable table */
    RevPtr<RevVariable> theRef = new RevVariable( the_var );
    variableTable.insert( std::pair<std::string, RevPtr<RevVariable> >( name, theRef ) );
    theRef->setName( name );
    
}


/** Add variable to frame. */
void Environment::addVariable( const std::string& name, const RevPtr<RevVariable>& the_var )
{
    
    /* Throw an error if the name string is empty. */
    if ( name == "" )
    {
        throw RbException("Invalid attempt to add unnamed variable to frame.");
    }
    
    /* Throw an error if the variable exists. Note that we cannot use the function
        existsVariable because that function looks recursively in parent frames, which
        would make it impossible to hide global variables. */
    if ( variableTable.find( name ) != variableTable.end() )
    {
        throw RbException() << "Variable " << name << " already exists in frame" ; 
    }
    
    /* Insert new RevVariable in variable table */
    variableTable.insert( std::pair<std::string, RevPtr<RevVariable> >( name, the_var ) );
    the_var->setName( name );

}


/** Add variable to frame from a "naked" Rev object. */
void Environment::addVariable(const std::string& name, RevObject* obj)
{
    // Create a new RevVariable object
    RevPtr<RevVariable> var = new RevVariable( obj ) ;

    // Add the variable to the table
    addVariable( name, var );
}


/** Clone the environment (frame) */
Environment* Environment::clone() const
{
    return new Environment(*this);
}


/**
 * Clear the variable table. The function table has its own deep destructor called by assignment or destruction.
 * However, we empty it here in case the clear function is called from some place other than the assignment
 * operator or the destructor.
 */
void Environment::clear(void)
{

    // Empty the variable table. It is as easy as this because we use smart pointers...
    variableTable.clear();

    // Empty the function table.
    function_table.clear();
}


/** Erase a variable by name given to it in the frame. */
void Environment::eraseVariable(const std::string& name)
{
    std::map<std::string, RevPtr<RevVariable> >::iterator it = variableTable.find(name);

    if ( it == variableTable.end() )
    {
        throw RbException(RbException::MISSING_VARIABLE, "Variable " + name + " does not exist in environment");
    }
    else
    {
        // Free the memory for the variable (smart pointer, so happens automatically) and
        // remove the variable from the map of variables
        variableTable.erase(it);
    }
}


/**
 * Erase a variable by its address. We just delegate the call to erase by name.
 */
void Environment::eraseVariable(const RevPtr<RevVariable>& var)
{
    
    // Delegate call
    eraseVariable( var->getName() );
}


/** Find whether a function name exists (current and enclosing frames) */
bool Environment::existsFunction(std::string const &name) const
{
    // We delegate the query to the function table
    return function_table.existsFunction(name);
}


/** Does variable exist in the environment (current frame and enclosing frames)? */
bool Environment::existsVariable(const std::string& name) const
{
    if (variableTable.find(name) == variableTable.end())
    {
        
        if (parentEnvironment != NULL)
        {
            return parentEnvironment->existsVariable(name);
        }
        else
        {
            return false;
        }
        
    }

    return true;
}


/** Does variable exist in the the current frame? */
bool Environment::existsVariableInFrame(const std::string& name) const
{
    if (variableTable.find(name) == variableTable.end())
        return false;
    else
        return true;
}


/** Generate a unique variable name */
std::string Environment::generateUniqueVariableName(void)
{
    std::ostringstream theName;

    theName << "var" << ++numUnnamedVariables;

    return theName.str();
}


std::shared_ptr<Environment> Environment::getParentEnvironment()
{
    return parentEnvironment;
}


std::shared_ptr<const Environment> Environment::getParentEnvironment() const
{
    return parentEnvironment;
}


std::shared_ptr<Environment> Environment::getChildEnvironment(const std::string &name)
{
    auto it = children.find(name);
    if ( it == children.end() )
    {
        std::shared_ptr<Environment> null;
        auto env = std::make_shared<Environment>(null, name);
        children.insert( {name, env} );
        return env;
    }
    else
    {
        return it->second;
    }
}


/** Get function. This call will throw an error if the name is missing or present multiple times. */
Function* Environment::getFunction(const std::string& name)
{
    return function_table.getFunction(name);
}


/* Get function. This call will throw an error if the function is missing. */
const Function& Environment::getFunction(const std::string& name, const std::vector<Argument>& args, bool once) const
{
    
    return function_table.getFunction(name, args, once);
}


/* Get function. This call will throw an error if the function is missing. */
const Function* Environment::findFunction(const std::string& name, const std::vector<Argument>& args, bool once) const
{
    return function_table.findFunction(name, args, once);
}


/** Return the function table (const) */
const FunctionTable& Environment::getFunctionTable(void) const
{

    return function_table;
}


/** Return the function table */
FunctionTable& Environment::getFunctionTable(void)
{

    return function_table;
}


/** Return the object inside a specific variable */
RevObject& Environment::getRevObject(const std::string& name)
{
    
    return getVariable( name )->getRevObject();
}


/** Return a specific variable (const version) */
const RevObject& Environment::getRevObject(const std::string& name) const
{
    
    return getVariable( name )->getRevObject();
}


/** Return a specific variable */
RevPtr<RevVariable>& Environment::getVariable(const std::string& name)
{
    std::map<std::string, RevPtr<RevVariable> >::iterator it = variableTable.find( name );
    
    if ( variableTable.find(name) == variableTable.end() )
    {
        
        if ( parentEnvironment != NULL )
        {
            return parentEnvironment->getVariable( name );
        }
        else
        {
            throw RbException(RbException::MISSING_VARIABLE, "Variable " + name + " does not exist");
        }
        
    }
    
    return it->second;
}


/** Return a specific variable (const version) */
const RevPtr<RevVariable>& Environment::getVariable(const std::string& name) const
{
    auto it = variableTable.find(name);
    
    if ( variableTable.find(name) == variableTable.end() )
    {
        if ( parentEnvironment )
        {
            getParentEnvironment()->getVariable( name );
        }
        else
        {
            throw RbException(RbException::MISSING_VARIABLE, "Variable " + name + " does not exist");
        }
        
    }
    
    return it->second;
}


/** Return variable table */
VariableTable& Environment::getVariableTable(void) {
    
    return variableTable;
}


/** Return variable table (const) */
const VariableTable& Environment::getVariableTable(void) const
{
    
    return variableTable;
}



bool Environment::hasChildEnvironment(const std::string &name)
{
    auto it = children.find(name);

    return it != children.end();
}


/**
 * Is the function with the name 'fxnName' a procedure? Simply ask the
 * function table.
 */
bool Environment::isProcedure(const std::string& fxnName) const
{
    return function_table.isProcedure( fxnName );
}

/**
 * Is Environment same or parent of otherEnvironment? We use this function
 * to decide when a reference from otherEnvironment to a variable in this
 * Environment is safe, and when it is not. The only time we know for sure
 * that it is safe is when this Environment is identical to, or a parent of,
 * otherEnvironment.
 */
bool Environment::isSameOrParentOf(const Environment& otherEnvironment) const
{

    if ( this == &otherEnvironment )
    {
        return true;
    }
    
    if ( otherEnvironment.parentEnvironment == NULL )
    {
        return false;
    }
    
    return isSameOrParentOf( *otherEnvironment.parentEnvironment );
}


/** Print workspace */
void Environment::printValue(std::ostream& o) const
{
    o << "Variable table:" << std::endl;

    VariableTable::const_iterator it;
    for ( it = variableTable.begin(); it != variableTable.end(); it++ )
    {
        o << (*it).first << " = ";
        (*it).second->printValue( o, true );
        o << std::endl;
    }

    o << std::endl;
}

