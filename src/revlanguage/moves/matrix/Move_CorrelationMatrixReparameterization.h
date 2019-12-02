#ifndef Move_CorrelationMatrixReparameterization_H
#define Move_CorrelationMatrixReparameterization_H

#include "CorrelationMatrixReparameterizationMove.h"
#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    /**
     * @brief Rev wrapper class for the rate-age-beta-shift move.
     *
     * This class is the Rev wrapper class for the rate-age-beta-shift move,
     * a move that proposes a new node age and branch rates so that the actual
     * branch lengths remain constant. The new node age is proposed from a scaled
     * beta distribution between its parent age and its oldest child age.
     *
     *
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @copyright GPL version 3
     */
    class Move_CorrelationMatrixReparameterization : public Move {
        
    public:
        
        Move_CorrelationMatrixReparameterization(void);                                                                                                        //!< Default constructor
        
        // Basic utility functions
        virtual Move_CorrelationMatrixReparameterization*              clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal Move.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   correlation_matrix;
        RevPtr<const RevVariable>                   mvbm;
        RevPtr<const RevVariable>                   variance;

        
    };
    
}

#endif
