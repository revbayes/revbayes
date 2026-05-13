#ifndef Dist_PseudoData_h
#define Dist_PseudoData_h

#include "PseudoDist.h"
#include "RlTypedDistribution.h"
#include "RlPseudoDataLikelihood.h"
#include "RlPseudoObservation.h"

namespace RevLanguage
{

    class Dist_Pseudo: public TypedDistribution<PseudoObservation>
    {
    public:
        Dist_Pseudo() {}

        Dist_Pseudo* clone() const override
        {
            return new Dist_Pseudo(*this);
        }

        //!< Get Rev type
        static const std::string&                       getClassType(void)
        {
            static std::string rev_type = "PseudoDist";

            return rev_type;
        }

        //!< Get class type spec
        static const TypeSpec&                          getClassTypeSpec(void)
        {
            static TypeSpec rev_type_spec ( TypedDistribution<PseudoObservation>::getClassTypeSpec() );

            return rev_type_spec;
        }

        //!< Get the type spec of the instance
        const TypeSpec&                                 getTypeSpec(void) const
        {
            static TypeSpec ts = getClassTypeSpec();

            return ts;
        }

        /**
         * Get the Rev name for the distribution.
         * This name is used for the constructor and the distribution functions,
         * such as the density and random value function
         *
         * \return Rev name of constructor function.
         */
        //!< Get the Rev-name for this distribution.
        std::string                                     getDistributionFunctionName(void) const override
        {
            return "Pseudo";
        }

        const MemberRules&                              getParameterRules(void) const override;


        // Distribution functions you have to override
        RevBayesCore::PseudoDist*                       createDistribution() const override;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var) override;       //!< Set member variable

    private:
        RevPtr<const RevVariable> pseudoDataLikelihood;
    };

}


#endif
