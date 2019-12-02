#ifndef MemberProcedure_H
#define MemberProcedure_H

#include "Procedure.h"
#include "RlMemberMethod.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace RevLanguage {
    
    class ArgumentRule;
    class RevObject;
    
    
    /**
     * \brief Member procedures which call a specific procedure of the object.
     *
     * This class simply manages the function call from Rev to a member procedure,
     * for example, "a.methods()". Internally we simply call executeMethod
     * for the instance of the object and pass it the method name.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     *
     */
    class MemberProcedure : public Procedure, public MemberMethod {
        
    public:
        MemberProcedure(const std::string &name, const TypeSpec retType, ArgumentRules* argRules);                                               //!< Constructor
        MemberProcedure(const MemberProcedure &mp);
        virtual                                            ~MemberProcedure();
        
        MemberProcedure&                                    operator=(const MemberProcedure &mp);
        
        // Basic utility functions
        MemberProcedure*                                    clone(void) const;                                                          //!< Clone the object
        static const std::string&                           getClassType(void);                                                         //!< Get class name
        static const TypeSpec&                              getClassTypeSpec(void);                                                     //!< Get class type spec
        const TypeSpec&                                     getTypeSpec(void) const;                                                    //!< Get language type of the object
        
        // Regular functions
        virtual RevPtr<RevVariable>                         execute(void);                                                              //!< Execute function
        const ArgumentRules&                                getArgumentRules(void) const;                                               //!< Get argument rules
        std::string                                         getFunctionName(void) const;                                                //!< Get the name of the function
        const TypeSpec&                                     getReturnType(void) const;                                                  //!< Get type of return value
        bool                                                isProcedure(void) const;                                                    //!< Is this a procedure?
        void                                                setMemberObject(const RevPtr<RevVariable> &obj);                            //!< Set the member object to which this function belongs
        
        
    protected:

        
        ArgumentRules*                                      argument_rules;                                                              //!< Argument rules (different for different member functions)
        std::string                                         proc_name;
        RevPtr<RevVariable>                                 object;                                                                     //!< The member object to which this function belongs (we do not own the member object because of cyclic ownership)
        TypeSpec                                            return_type;                                                                 //!< Return type (different for different member functions)
        
    };
    
}

#endif

