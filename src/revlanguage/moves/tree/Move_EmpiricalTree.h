#ifndef Move_EmpiricalTree_H
#define Move_EmpiricalTree_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
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
        std::vector<std::string>                    getHelpAuthor(void) const;                                  //!< Get the author(s) of this function
        std::string                                 getHelpDescription(void) const;                             //!< Get the description for this function
        std::string                                 getHelpExample(void) const;                                 //!< Get an executable and instructive example
        std::vector<std::string>                    getHelpSeeAlso(void) const;                                 //!< Get suggested other functions
        std::string                                 getMoveName(void) const;                                                                        //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                                  //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        RevPtr<const RevVariable>                   tree;
        RevPtr<const RevVariable>                   metropolisHastings;
    };
    
}


#endif
