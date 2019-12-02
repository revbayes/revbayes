#ifndef Move_SynchronizedVectorFixedSingleElementSlide_H
#define Move_SynchronizedVectorFixedSingleElementSlide_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    class Move_SynchronizedVectorFixedSingleElementSlide : public Move {
        
    public:
        
        Move_SynchronizedVectorFixedSingleElementSlide(void);                                                                                           //!< Default constructor
        
        // Basic utility functions
        virtual Move_SynchronizedVectorFixedSingleElementSlide*     clone(void) const;                                                                      //!< Clone object
        void                                                        constructInternalObject(void);                                                          //!< We construct the a new internal SlidingMove.
        static const std::string&                                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                                   v;
        RevPtr<const RevVariable>                                   lambda;
        RevPtr<const RevVariable>                                   weight;
        RevPtr<const RevVariable>                                   which_element;
        RevPtr<const RevVariable>                                   tune;
        
    };
    
}

#endif
