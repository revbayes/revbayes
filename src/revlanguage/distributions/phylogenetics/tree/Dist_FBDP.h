#ifndef Dist_FBDP_H
#define Dist_FBDP_H

#include "Dist_BDSTP.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the Fossilized-Birth-Death Process
     *
     * Helper distribution of the BDSTP for fossil applications, with fixed removal probability = 0 (i.e. no death after sampling)
     */
    class Dist_FBDP : public Dist_BDSTP {

    public:

        // Basic utility functions
        Dist_FBDP*                                              clone(void) const override;                                                             //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        std::vector<std::string>                                getDistributionFunctionAliases(void) const override;                                    //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const override;                                       //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const override;                                                       //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const override;                                                 //!< Get member rules (const)

    protected:

        virtual void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) override; //!< Set member variable
        virtual RevBayesCore::DagNode*                          getRemovalProbability( void ) const override;
    };

}

#endif
