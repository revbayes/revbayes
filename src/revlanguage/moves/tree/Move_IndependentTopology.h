#ifndef Move_IndependentTopology_H
#define Move_IndependentTopology_H

#include <ostream>

#include "RlMove.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
class TypeSpec;
    
    /**
     * @brief Rev Wrapper of a constrained tree distribution move.
     *
     * This class is the RevLanguage wrapper of EmpiricalTreeMove.
     *
     * @author The RevBayes Development Core Team 
     * @copyright GPL version 3
     */
    class Move_IndependentTopology : public Move {
        
    public:
        
        Move_IndependentTopology(void);                                                                                                                   //!< Default constructor
        
        // Basic utility functions
        virtual Move_IndependentTopology*           clone(void) const;                                                                              //!< Clone the object
        void                                        constructInternalObject(void);                                                                  //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                        //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                                  //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        RevPtr<const RevVariable>                   outgroup;
        RevPtr<const RevVariable>                   proposal_distribution;
        RevPtr<const RevVariable>                   substitution_rates;
        RevPtr<const RevVariable>                   tree;
    };
    
}


#endif
