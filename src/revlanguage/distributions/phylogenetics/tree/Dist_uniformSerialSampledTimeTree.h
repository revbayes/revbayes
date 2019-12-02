#ifndef Dist_UniformSerialSampledTimeTree_H
#define Dist_UniformSerialSampledTimeTree_H

#include "UniformSerialSampledTimeTreeDistribution.h"
#include "RlTimeTree.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * @file
     * This file contains the declaration of the Rev uniform time tree distribution
     *
     * @brief Declaration of the Rev uniform time tree distribution
     *
     * @author Fredrik Ronquist
     */
 class Dist_uniformSerialSampledTimeTree :  public TypedDistribution<TimeTree> {
        
    public:
                Dist_uniformSerialSampledTimeTree( void );        //!< Constructor
        virtual ~Dist_uniformSerialSampledTimeTree();             //!< Virtual destructor
        
        // Basic utility functions
        Dist_uniformSerialSampledTimeTree*              clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)
        
        
        // Distribution functions you have to override
        RevBayesCore::UniformSerialSampledTimeTreeDistribution* createDistribution(void) const;                                                         //!< Create the core object corresponding to this wrapper
        
    protected:
        
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        
    private:
        
        RevPtr<const RevVariable>                       taxa;
        RevPtr<const RevVariable>                       start_age;
        std::string                                     start_condition;                                                                        //!< The start condition of the process (rootAge/originAge)

    };

}

#endif
