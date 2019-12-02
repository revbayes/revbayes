#ifndef Move_NodeTimeSlidePathTruncatedNormal_H
#define Move_NodeTimeSlidePathTruncatedNormal_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    
    /**
     * @brief Rev wrapper class for the NodeTimeSlidePathTruncatedNormal move.
     *
     * This class is the Rev wrapper class for the NodeTimeSlidePathTruncatedNormal move,
     * a move that proposes a node age.
     *
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     */
    class Move_NodeTimeSlidePathTruncatedNormal : public Move {
        
    public:
        
        Move_NodeTimeSlidePathTruncatedNormal(void);                                                                                                    //!< Default constructor
        
        // Basic utility functions
        virtual Move_NodeTimeSlidePathTruncatedNormal*          clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal Move.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   tree;
        RevPtr<const RevVariable>                   tune;
        RevPtr<const RevVariable>                   sigma;
        RevPtr<const RevVariable>                   scaleByAge;
        
    };
    
}

#endif
