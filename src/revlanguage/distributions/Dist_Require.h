#ifndef Dist_RequireData_h
#define Dist_RequireData_h

#include "RequireDist.h"
#include "RlTypedDistribution.h"
#include "RlPseudoObservation.h"

namespace RevLanguage
{

    class Dist_Require: public TypedDistribution<PseudoObservation>
    {
    public:
        Dist_Require() {}

        Dist_Require* clone() const
        {
            return new Dist_Require(*this);
        }

        //!< Get Rev type
        static const std::string&                       getClassType(void)
        {
            static std::string rev_type = "RequireDist";

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
        std::string                                     getDistributionFunctionName(void) const
        {
            return "Require";
        }

        const MemberRules&                              getParameterRules(void) const;


        // Distribution functions you have to override
        RevBayesCore::RequireDist*                      createDistribution() const override;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

    private:
        RevPtr<const RevVariable> predicate;
        RevPtr<const RevVariable> weight;
    };

}


#endif
