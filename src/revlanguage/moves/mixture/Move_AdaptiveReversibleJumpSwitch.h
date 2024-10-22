#ifndef Move_AdaptiveReversibleJumpSwitch_H
#define Move_AdaptiveReversibleJumpSwitch_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of an adaptive reversible jump move between a constant value and
     *  a random distributedly value.
     *
     * This class is the RevLanguage wrapper of ReversibleJumpMixtureProposal.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2014-08-18, version 1.0
     */
    class Move_AdaptiveReversibleJumpSwitch : public Move {
        
    public:
        
        Move_AdaptiveReversibleJumpSwitch(void);                                                                                                            //!< Default constructor
        
        // Basic utility functions
        virtual Move_AdaptiveReversibleJumpSwitch*  clone(void) const;                                                                              //!< Clone the object
        void                                        constructInternalObject(void);                                                                  //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                        //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                                  //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        RevPtr<const RevVariable>                   x;                                                                                              //!< The variable holding the real valued vector.
        RevPtr<const RevVariable>                   wait_before_learning;                                                                                              //!< The variable holding the real valued vector.
        RevPtr<const RevVariable>                   wait_before_using;                                                                                              //!< The variable holding the real valued vector.
        RevPtr<const RevVariable>                   update_every;                                                                                              //!< The variable holding the real valued vector.

    };
    
}


#endif
