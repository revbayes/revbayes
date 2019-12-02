//
//  Move_MatrixRealSymmetricSlide.h
//  revbayes
//
//  Created by Nicolas Lartillot on 2014-03-28.
//  Copyright (c) 2014 revbayes team. All rights reserved.
//

#ifndef Move_MatrixRealSymmetricSlide_H
#define Move_MatrixRealSymmetricSlide_H

#include <iostream>

#include "RlMove.h"
#include "TypedDagNode.h"


namespace RevLanguage {
    
    class Move_MatrixRealSymmetricSlide : public Move {
        
    public:
        
        Move_MatrixRealSymmetricSlide(void);                                                                                                //!< Default constructor
        
        // Basic utility functions
        virtual Move_MatrixRealSymmetricSlide*      clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal SlidingMove.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   mat;
        RevPtr<const RevVariable>                   delta;
        RevPtr<const RevVariable>                   weight;
        RevPtr<const RevVariable>                   tune;
        
    };
    
}


#endif /* defined(__revbayes__Move_MatrixRealSymmetricSlide__) */
