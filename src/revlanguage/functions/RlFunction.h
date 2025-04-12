#ifndef RlFunction_H
#define RlFunction_H

#include <string>
#include <vector>
#include <memory>

#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Environment.h"
#include "RbHelpFunction.h"
#include "RevPtr.h"
#include "TypeSpec.h"


namespace RevLanguage {

    /**
     * \brief Function: Abstract base class for Rev functions
     *
     * This is the interface and abstract base class for functions in
     * Rev. Function instances are put in the function table in the
     * relevant Environment (user workspace or base environment) if they are
     * global. If they are member functions of user-defined type, they
     * are instead associated with the function table of the approprioate
     * class in the class table of the user workspace.
     *
     * A function instance knows its argument rules and can process a
     * vector of labeled argument values according to these rules to pro-
     * duce a ready-to-use list of arguments. The processing involves label
     * matching as well as filling of missing values with default values.
     * The processing of labeled argument values is done in the function
     * processArguments, which will throw an error if the provided
     * arguments do not match.
     */
    class Function : public RevObject {
    public:
        virtual                                        ~Function(void);                                                                     //!< Destructor
                                                        Function(const Function &x);                                                        //!< Copy constuctor

        // Basic utility functions you have to override
        virtual Function*                               clone(void) const = 0;                                                              //!< Clone object
        static const std::string&                       getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                             //!< Get class type spec

        // Basic utility functions you may want to override
        virtual std::string                             callSignature(void) const;                                                          //!< Return call signature

        // Basic utility functions you should not have to override
        std::string                                     getHelpUsage(void) const;
        std::string                                     getRevDeclaration(void) const;                                                      //!< Get Rev declaration of the function
        void                                            printValue(std::ostream& o, bool user) const;                                       //!< Print the general information on the function ('usage')
        void                                            setExecutionEnviroment(const std::shared_ptr<Environment>& e);                                             //!< Set the environment from which the function was called.
    
        // Functions you have to override
        virtual RevPtr<RevVariable>                     execute(void) = 0;                                                                  //!< Execute function
        virtual const ArgumentRules&                    getArgumentRules(void) const = 0;                                                   //!< Get argument rules
        virtual std::string                             getConstructorFunctionName(void) const;
        virtual std::vector<std::string>                getConstructorFunctionAliases(void) const;
        virtual std::string                             getFunctionName(void) const = 0;                                                    //!< Get the name of the function
        virtual std::vector<std::string>                getFunctionNameAliases(void) const;                                                 //!< Get aliases for the Rev name of the function
        virtual const TypeSpec&                         getReturnType(void) const = 0;                                                      //!< Get type of return value

        
        
        // Functions you may want to override
        virtual bool                                    checkArguments(const std::vector<Argument>& passed_args,
                                                                       std::vector<double>*         matchScore,
                                                                       std::vector<bool>&           arg_mapped,
                                                                       bool                         once);                                  //!< Process args, return a match score if pointer is not null
        virtual bool                                    isInternal(void) const { return false; }                                            //!< Is the function a procedure?
        virtual bool                                    isProcedure(void) const { return false; }                                           //!< Is the function a procedure?
        virtual void                                    processArguments(const std::vector<Argument>& passed_args);                         //!< Process args, return a match score if pointer is not null


        // Function functions you should not override
        void                                            clear(void);                                                                        //!< Clear argument frame "args"
        const std::vector<Argument>&                    getArguments(void) const;                                                           //!< Get processed arguments in argument Environment "args"
        std::vector<Argument>&                          getArguments(void);                                                                 //!< Get processed arguments in argument Environment "args"
        std::shared_ptr<Environment>                    getEnvironment(void) const;                                                         //!< Get the execution environment
        
        
        
	protected:
                                                        Function(void);                                                                     //!< Basic constructor
        
        // virtual method that may be overwritten
        virtual RevBayesCore::RbHelpFunction*           constructTypeSpecificHelp(void) const;
        virtual void                                    addSpecificHelpFields(RevBayesCore::RbHelpEntry* e) const;

        // Function you may want to override
        virtual void                                    clearArguments(void);                                                               //!< Clear arguments

        // Member variables
        bool                                            argsProcessed;                                                                      //!< Are arguments processed?
        std::vector<Argument>                           args;                                                                               //!< Vector of arguments
        std::shared_ptr<Environment>                    env;                                                                                //!< Evaluation environment

    private:
        double                                          computeMatchScore(const RevVariable* arg, const ArgumentRule& rule);
};
    
}

#endif
