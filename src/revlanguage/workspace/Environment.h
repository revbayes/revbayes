#ifndef Environment_H
#define Environment_H

#include <cstddef>
#include <map>
#include <string>
#include <vector>
#include <iosfwd>

#include "FunctionTable.h"
#include "RevPtr.h"
#include "RevVariable.h"


namespace RevLanguage {
class Argument;
class Function;
class RevObject;

    typedef std::map<std::string, RevPtr<RevVariable> > VariableTable;                                                             //!< Typedef for convenience

    /**
     * @brief Environment: Base class for frames
     *
     * Instances of the Environment class are used to hold information about an
     * evaluation or execution environment, such as the stack frame of a function call
     * or a namespace.
     *
     * (a) Environment objects that represent stack frames only define local variables and
     *     contain a pointer to the enclosing (parent) environment to look up external variables.
     *     Such Environments contain a reference to the parent Environment so that the parent
     *     Environment will not be destroyed before the stack frame Environment.
     *
     * (b) Environment objects that represent a namespace do not contain a pointer to the parent
     *     Environment.  Instead, the parent Environment holds a references to the namespace.
     *     This ensures that the namespace Environment will not be destroyed before the
     *     enclosing Environment.
     *
     * The base environment is the global workspace. It is a
     * special type of environment, which is described in the class
     * Workspace, derived from Environment.
     *
     * An environment essentially consists of a table (map) or vector of
     * variables, which may or may not be labeled. The workspace
     * also contains a function table, and a class (type) table.
     *
     */

    class Environment {
        
    public:
        Environment(const std::string &n);                                                                                                              //!< Constructor of Environment with NULL parent
        Environment(const std::shared_ptr<Environment>& parentFr, const std::string &n);                                                                                             //!< Constructor of Environment with parent
        Environment(const Environment& x);                                                                                              //!< Copy Constructor
        virtual ~Environment(void);                                                                                                     //!< Destrcutor

        // Operators
        Environment&                        operator=(const Environment& x);                                                            //!< Assignment

        // Basic utility functions
        virtual Environment*                clone(void) const;                                                                          //!< Clone Environment
        virtual void                        printValue(std::ostream& o) const;                                                          //!< Print table for user

        // Regular functions
        void                                addAlias(const std::string& name, const RevPtr<RevVariable>& var);                          //!< Add alias of variable
        bool                                addFunction(Function* func);                                       //!< Add function
        void                                addNullVariable(const std::string& name);                                                   //!< Add variable with null object
        void                                addReference(const std::string& name, const RevPtr<RevVariable>& var);                      //!< Add reference to variable
        void                                addVariable(const std::string& name, const RevPtr<RevVariable>& var);                       //!< Add variable
        void                                addVariable(const std::string& name, RevObject* object);                                    //!< Add variable from object
        void                                clear(void);                                                                                //!< Clear the environment (variables and functions)
        void                                eraseVariable(const std::string& name);                                                     //!< Erase a variable by name
        void                                eraseVariable(const RevPtr<RevVariable>& var);                                              //!< Erase a variable by variable address
//        RevPtr<RevVariable>                    executeFunction(const std::string& name, const std::vector<Argument>& args);                //!< Execute function
        bool                                existsFunction(const std::string &name) const;                                              //!< Does a function exists with given name?
        bool                                existsVariable(const std::string& name) const;                                              //!< Does variable exist?
        bool                                existsVariableInFrame(const std::string& name) const;                                       //!< Does variable exist in this frame?
        std::string                         generateUniqueVariableName(void);                                                           //!< Automatically generate a unique variable name

        std::shared_ptr<Environment>        getParentEnvironment();                                                                     //!< Get parent environment
        std::shared_ptr<const Environment>  getParentEnvironment() const;                                                              //!< Get parent environment
        std::shared_ptr<Environment>        getChildEnvironment(const std::string &name);                                               //!< Get child environment by name
        Function*                           getFunction(const std::string& name);                                                       //!< Get function reference
        const Function&                     getFunction(const std::string& name, const std::vector<Argument>& args, bool once) const;   //!< Get function reference
        const Function*                     findFunction(const std::string& name, const std::vector<Argument>& args, bool once) const;  //!< Get function reference
        const FunctionTable&                getFunctionTable(void) const;                                                               //!< Get function table (const)
        FunctionTable&                      getFunctionTable(void);                                                                     //!< Get function table (non-const)
        const RevObject&                    getRevObject(const std::string& name) const;                                                //!< Convenient alternative for [name]->getValue()
        RevObject&                          getRevObject(const std::string& name);                                                      //!< Convenient alternative for [name]->getValue() (non-const to return non-const value)
        RevPtr<RevVariable>&                getVariable(const std::string& name);                                                       //!< Get variable
        const RevPtr<RevVariable>&          getVariable(const std::string& name) const;                                                 //!< Get variable (const)
        const VariableTable&                getVariableTable(void) const;                                                               //!< Get the table with the variables (const)
        VariableTable&                      getVariableTable(void);                                                                     //!< Get the table with the variables (non-const)
        bool                                hasChildEnvironment(const std::string &name);                                               //!< Has a child environment with the name
        bool                                isProcedure(const std::string& fxnName) const;                                              //!< Is 'fxnName' a procedure?
        virtual bool                        isSameOrParentOf(const Environment& otherEnvironment) const;                                //!< Is the Environment same or parent of other Environment?
        void                                setParentEnvironment(const std::shared_ptr<Environment>& newEnvironment);                   //!< Set parent Environment
        size_t                              size(void) const;                                                                           //!< Get size of variable table
        

    protected:

        FunctionTable                       function_table;                                                                              //!< Table holding functions
        int                                 numUnnamedVariables;                                                                        //!< Current number of unnamed variables
        std::shared_ptr<Environment>        parentEnvironment;                                                                          //!< Pointer to enclosing Environment
        VariableTable                       variableTable;                                                                              //!< Variable table
    
        std::map<std::string, std::shared_ptr<Environment>> children;
        std::string                         name;
    };

}

#endif
