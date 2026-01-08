#ifndef Move_EmpiricalTree_H
#define Move_EmpiricalTree_H

#include <ostream>
#include <vector>

#include "RlMove.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
class TypeSpec;
    
    /**
     * @brief Rev Wrapper of an empirical tree distribution move.
     *
     * This class is the RevLanguage wrapper of EmpiricalTreeMove.
     *
     * @author The RevBayes Development Core Team 
     * @copyright GPL version 3
     */
    class Move_EmpiricalTree : public Move {
        
    public:
        
        Move_EmpiricalTree(void);                                                                                                                   //!< Default constructor
        
        // Basic utility functions
        virtual Move_EmpiricalTree*                 clone(void) const;                                                                              //!< Clone the object
        void                                        constructInternalObject(void);                                                                  //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                        //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                                  //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        RevPtr<const RevVariable>                   tree;
    };
    
}


#endif
