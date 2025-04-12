#ifndef ArgumentRule_H
#define ArgumentRule_H

#include "Environment.h"
#include "RevPtr.h"
#include "TypeSpec.h"

#include <string>

namespace RevLanguage {

    /**
     * @brief Argument rule class
     *
     * The argument rule class describes the rules for fitting an
     * argument (a variable) into an argument or a member variable
     * slot.
     */
    class ArgumentRule {

    public:
        
        enum EvaluationType { BY_VALUE, BY_REFERENCE, BY_CONSTANT_REFERENCE };
        enum DagNodeType { STOCHASTIC, DETERMINISTIC, ANY, DYNAMIC, CONSTANT };
        
        ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::string& argName, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::string& argName, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::string& argName, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::vector<std::string>& argNames, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::vector<std::string>& argNames, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::vector<std::string>& argNames, const TypeSpec& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal);                                                                               //!< Constructor of rule without default value
        ArgumentRule(const std::vector<std::string>& argNames, const std::vector<TypeSpec>& argTypeSp, const std::string& argDesc, EvaluationType et, DagNodeType dt, RevObject *defVal);                                                                               //!< Constructor of rule without default value

        virtual                            ~ArgumentRule(void) {}
        
        // Basic utility functions
        virtual ArgumentRule*               clone(void) const;                                                              //!< Clone object
        virtual void                        printValue(std::ostream& o) const;                                              //!< Print value for user

        // ArgumentRule functions
        Argument                            fitArgument( Argument& arg, bool once ) const;                                  //!< Fit a passed argument into a slot using the rule
        const std::vector<std::string>&     getArgumentAliases(void) const;                                                 //!< Get aliases of argument
        DagNodeType                         getArgumentDagNodeType(void) const;                                             //!< Get the DAG node type
        const std::string&                  getArgumentLabel(void) const;                                                   //!< Get label of argument
        const std::vector<TypeSpec>&        getArgumentTypeSpec(void) const;                                                //!< Get argument type spec
        virtual const RevVariable&          getDefaultVariable(void) const;                                                 //!< Get default argument
        const std::string&                  getArgumentDescription(void) const;                                             //!< Get label of argument
        EvaluationType                      getEvaluationType(void) const;                                                  //!< Get the evaluation type
        bool                                hasDefault(void) const;                                                         //!< Does the rule have a default value?
        virtual double                      isArgumentValid(Argument &arg) const;                                           //!< Is var a valid argument?
        virtual bool                        isEllipsis(void) const;                                                         //!< Is this an ellipsis?
    

    protected:
        
        std::vector<TypeSpec>               argTypeSpecs;                                                                   //!< Type specs
        RevPtr<RevVariable>                 defaultVar;                                                                     //!< Default value
        EvaluationType                      evalType;                                                                       //!< Is rule const?
        DagNodeType                         nodeType;
        std::vector<std::string>            aliases;                                                                        //!< Aliases of argument
        std::string                         label;                                                                          //!< Selected label of argument
        std::string                         description;                                                                    //!< Description of argument
        bool                                hasDefaultVal;                                                                  //!< Has default

    };
    
}

#endif
