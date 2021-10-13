#ifndef Dist_FBDPSpeciesMatrix_H
#define Dist_FBDPSpeciesMatrix_H

#include "FossilizedBirthDeathMatrixProcess.h"
#include "RlFossilizedBirthDeathSpeciesProcess.h"
#include "ModelVector.h"
#include "RlMatrixReal.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the Fossilized-Birth-Death Range Matrix Process
     *
     * The RevLanguage wrapper of the fossilized-birth-death range matrix process connects
     * the variables/parameters of the process and creates the internal FossilizedBirthDeathRangeMatrixProcess object.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *c
     */
    class Dist_FBDPSpeciesMatrix : public FossilizedBirthDeathSpeciesProcess<MatrixReal> {
        
    public:
        Dist_FBDPSpeciesMatrix( void );
        
        // Basic utility functions
        Dist_FBDPSpeciesMatrix*                                 clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::FossilizedBirthDeathMatrixProcess*        createDistribution(void) const;

    };
    
}

#endif
