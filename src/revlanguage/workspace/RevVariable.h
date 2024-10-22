#ifndef RevVariable_H
#define RevVariable_H

#include "RevObject.h"
#include "RevPtr.h"
#include "TypeSpec.h"

#include <string>
#include <set>
#include <boost/optional.hpp>

namespace RevLanguage {
    
    /**
     * @brief Variable: named (in frames) or unnamed (temporary) Rev variables
     *
     * A Variable instance corresponds to a variable in the Rev language.
     * For instance "a <- 1" will create a variable holding a pointer to
     * a Rev object containing a constant DAG node with value 1. A variable
     * is either unnamed (it is a temporary variable), or it is inserted into
     * a frame ("Environment"), in which case it will be named.
     *
     * Rev variables can either be standard variables, control variables or
     * reference variables. A standard variable contains a RevObject instance,
     * which in turn may have a stochastic, deterministic or stochastic DAG
     * node inside it. There are also other types of Rev objects that are
     * not associated with model DAG nodes: workspace objects and factor objects.
     *
     * A control variable is a standard variable flagged as being a control variable,
     * which causes it to behave differently in equation and tilde assignments. The
     * typical case is a loop variable. Finally, a reference variable is an alias to
     * another variable, and keeps a smart pointer to the variable it references.
     *
     * The type of variable is controlled by different assignment operators. The
     * standard operators ('<-', ':=' and '~') all create standard variables. A
     * control variable is created by the control assignment operator ('<<-') and
     * a reference variable by the reference assignment operator ('<-&').
     *
     * @author RevBayes core development team
     */
    class RevVariable {
    public:
        RevVariable(const TypeSpec& ts, const std::string& n="");                       //!< Constructor of empty variable with a type specification
        RevVariable(RevObject *revObj, const std::string &n = "");                      //!< Constructor of filled variable
        RevVariable(const RevPtr<RevVariable>& refVar, const std::string &n = "");      //!< Constructor of reference variable
        RevVariable(const RevVariable &v);                                              //!< Copy constructor
        
        virtual                 ~RevVariable(void);                                     //!< Virtual destructor
        
        RevVariable&            operator=(const RevVariable &v);                        //!< Assignment operator

        // Regular functions
        void                    addIndex(int idx, const RevPtr<RevVariable>& elem);     //!< Resize the vector to include this index.
        RevVariable*            clone(void) const;                                      //!< Clone variable
        size_t                  getMaxElementIndex(void) const;                         //!< Get the set of element indices for this vector variable.
        const std::string&      getName(void) const;                                    //!< Get the name of the variable
        RevObject&              getRevObject(void) const;                               //!< Get the value of the variable (non-const to return non-const value)
        const TypeSpec&         getRequiredTypeSpec(void) const;                        //!< Get the required Rev object type spec
        bool                    isAssignable(void) const;                               //!< Is variable or upstream members assignable?
        bool                    isHiddenVariable(void) const;                           //!< Is this a hidden variable?
        bool                    isElementVariable(void) const;                          //!< Is this a variable an element of a vector?
        bool                    isReferenceVariable(void) const;                        //!< Is this a reference variable?
        bool                    isVectorVariable(void) const;                           //!< Is this a vector variable?
        bool                    isWorkspaceVariable(void) const;                        //!< Is this a workspace variable?
        void                    makeReference(const RevPtr<RevVariable>& refVar);       //!< Make this a reference variable
        void                    printValue(std::ostream& o, bool toScreen) const;       //!< Print value of variable
        void                    setElementVariableState(bool flag = true);              //!< Set (or unset) element variable status
        void                    setHiddenVariableState(bool flag = true);               //!< Set (or unset) hidden variable status
        void                    setToVectorVariable();                                  //!< Set (or unset) vector variable status
        void                    setWorkspaceVariableState(bool flag = true);            //!< Set (or unset) control variable status
        void                    setName(const std::string &n);                          //!< Set the name of this variable
        void                    replaceRevObject(RevObject *newObj);                    //!< Replace the Rev object of this variable
        void                    setRequiredTypeSpec(const TypeSpec& ts);                //!< set the required Rev object type spec
        
        // Smart pointer functions
        size_t                  decrementReferenceCount(void) const;                    //!< Decrement the reference count for reference counting in smart pointers
        size_t                  getReferenceCount(void) const;                          //!< Get the reference count for reference counting in smart pointers
        void                    incrementReferenceCount(void) const;                    //!< Increment the reference count for reference counting in smart pointers
        
    private:
        // Member variables

        // variables related to building the components of a vector
        boost::optional<std::vector<RevPtr<RevVariable> > > vector_var_elements;          //!< Elements if this is a vector variable.
        mutable bool            needs_building = false;                                 //!< Do we need to construct the revobject?

        bool                    is_element_var = false;                                 //!< Is this variable an element of a vector?
        bool                    is_hidden_var = false;                                  //!< Is this a hidden variable?
        bool                    is_workspace_var = false;                               //!< Is this a workspace variable?
        std::string             name;                                                   //!< Name of variable
        mutable size_t          ref_count = 0;                                          //!< Reference count used by RevPtr
        RevPtr<RevVariable>     referenced_variable = nullptr;                          //!< Smart pointer to referenced variable
        RevObject*              rev_object = nullptr;                                   //!< Pointer to the Rev object inside the variable
        TypeSpec                required_type_spec;                                     //!< Required type of the object
    };
    
}


#endif
