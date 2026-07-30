#ifndef Move_RandomCategoryWalk_H
#define Move_RandomCategoryWalk_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a random category walk move.
     *
     * This class is the RevLanguage wrapper of RandomCategoryWalkProposal.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2022-08-04, version 1.2
     */
    class Move_RandomCategoryWalk : public Move {
        
    public:
        
        Move_RandomCategoryWalk(void);                                                                                                                     //!< Default constructor
        
        // Basic utility functions
        virtual Move_RandomCategoryWalk*            clone(void) const;                                                                              //!< Clone the object
        void                                        constructInternalObject(void);                                                                  //!< We construct the a new internal move.
        static const std::string&                   getClassType(void);                                                                             //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                         //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                        //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                                  //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                        //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                              //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);               //!< Set member variable
        
        RevPtr<const RevVariable>                   x;                                                                                              //!< The variable holding the vector.
        
    };
    
}


#endif /* defined(__revbayes_proj__Move_RandomCategoryWalk__) */
