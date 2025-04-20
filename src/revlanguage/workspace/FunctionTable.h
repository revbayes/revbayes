#ifndef FunctionTable_H
#define FunctionTable_H

#include "RevPtr.h"

#include <map>
#include <ostream>
#include <string>
#include <vector>

namespace RevLanguage {
    
    class Argument;
    class ArgumentRule;
    class ArgumentRules;
    class Function;
    class RevObject;
    class RevVariable;

    /**
     * @brief FunctionTable: A multimap from function names to functions
     *
     * FunctionTable is used to hold functions in Workspace and Environment (frame)
     * objects. It holds the functions, which it owns, in a std::multimap, which it
     * is derived from. Function tables can be nested; each table defers
     * calls to its parent(s) when the task cannot be solved locally.
     *
     */
    class FunctionTable : public std::multimap<std::string, Function*> {
        
    public:

        FunctionTable(FunctionTable* parent = NULL);                                                                                        //!< Empty table
        FunctionTable(const FunctionTable& x);                                                                                              //!< Copy constructor
        virtual                                 ~FunctionTable();                                                                           //!< Delete functions

        // Assignment operator
        FunctionTable&                          operator=(const FunctionTable& x);                                                          //!< Assignment operator

        // Basic utility functions
        virtual FunctionTable*                  clone(void) const;                                                                          //!< Clone object
        void                                    printValue(std::ostream& o, bool env) const;                                                //!< Print table for user

        // FunctionTable functions
        bool                                    addFunction(Function *func, bool no_except=false);                                          //!< Add function
        void                                    clear(void);                                                                                //!< Clear table
//        RevPtr<RevVariable>                        executeFunction(const std::string&           name,
//                                                                const std::vector<Argument>& args);                                         //!< Evaluate function (once)
        bool                                    existsFunction(const std::string &name) const;                                              //!< Does this table contain a function with given name?
        bool                                    existsFunctionInFrame(const std::string &name, const ArgumentRules& r) const;               //!< Does this table contain a function with given name?
        void                                    eraseFunction(const std::string& name);                                                     //!< Erase a function (all versions)
        std::vector<Function*>                  findFunctions(const std::string& name) const;                                               //!< Return functions matching name
        void                                    getFunctionNames(std::vector<std::string>& names) const;
        Function*                               getFirstFunction(const std::string& name) const;                                            //!< Get first function with given name
        Function*                               getFunction(const std::string& name) const;                                                 //!< Get function, throw an error if overloaded
        const Function*                         findFunction(const std::string& name, const std::vector<Argument>& args) const;             //!< Get function
        const Function&                         getFunction(const std::string& name, const std::vector<Argument>& args) const;              //!< Get function
        bool                                    isDistinctFormal(const ArgumentRules& x, const ArgumentRules& y) const;                     //!< Are formals unique?
        bool                                    isProcedure(const std::string& fxnName) const;                                              //!< Is 'fxnName' a procedure?
        void                                    replaceFunction(const std::string &name, Function* func);                                   //!< Replace existing function
        void                                    setParentTable(const FunctionTable* ft) { parentTable = ft; }                               //!< Set parent table

    protected:
        
        void                                    testFunctionValidity(const std::string& name, Function* func) const;                        //!< Test whether function can be added
        
        // Member variables
        const FunctionTable*                    parentTable;                                                                                //!< Enclosing table

};
    
}

#endif

