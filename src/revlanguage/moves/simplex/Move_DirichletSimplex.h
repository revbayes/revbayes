#ifndef Move_DirichletSimplex_H
#define Move_DirichletSimplex_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev Wrapper of a Dirichlet-simplex move on all elements of a simplex.
     *
     * This class is the RevLanguage wrapper of DirichletSimplexMove.
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     * @since 2015-05-21, version 1.0
     */
    class Move_DirichletSimplex : public Move {
        
    public:
        
        Move_DirichletSimplex(void);                                                                                                        //!< Default constructor
        
        // Basic utility functions
        virtual Move_DirichletSimplex*              clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal SlidingMove.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                    getMoveAliases(void) const;                                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   x;
        RevPtr<const RevVariable>                   alpha;
        RevPtr<const RevVariable>                   numCats;
        RevPtr<const RevVariable>                   offset;
        RevPtr<const RevVariable>                   kappa;
        RevPtr<const RevVariable>                   tune;
        RevPtr<const RevVariable>                   alphaUpper;
        
    };
    
}

#endif
