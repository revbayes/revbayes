#ifndef Dist_FBDRangeMatrix_H
#define Dist_FBDRangeMatrix_H

#include "ModelVector.h"
#include "FossilizedBirthDeathMatrixProcess.h"
#include "RlFossilizedBirthDeathRangeProcess.h"
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
    class Dist_FBDRangeMatrix : public FossilizedBirthDeathRangeProcess<MatrixReal> {
        
    public:
        Dist_FBDRangeMatrix( void );
        
        // Basic utility functions
        Dist_FBDRangeMatrix*                                    clone(void) const;                                                                      //!< Clone the object
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
