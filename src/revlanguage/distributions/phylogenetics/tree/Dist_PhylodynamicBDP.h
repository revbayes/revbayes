#ifndef Dist_PhylodynamicBDP_H
#define Dist_PhylodynamicBDP_H

#include "Dist_BDSTP.h"

namespace RevLanguage {

    /**
     * RevLanguage wrapper of the phylodynamic Birth-Death Sampling Treatment Process
     * 
     * Helper distribution of the BDSTP for phylodynamic applications, with removal probability = 1 (i.e. sampling = death) and extant sampling = 0 by default, 
     * and no bursts or mass extinctions
     * 
     */
    class Dist_PhylodynamicBDP : public Dist_BDSTP {

    public:
        // Basic utility functions
        Dist_PhylodynamicBDP*                                   clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const)

    protected:

        void                                                    setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
    };

}

#endif
