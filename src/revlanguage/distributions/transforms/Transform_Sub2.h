#ifndef Transform_Sub2_H
#define Transform_Sub2_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"
#include "TransformedDistribution.h"

namespace RevLanguage {

    /* This class creates a distribution that implements
     *
     *    value - distribution
     *
     * This is equivalent to distribution |> tnScale(-1) |> tnShift(value).
     *
     * Note that Transform_Sub1 and Transform_Sub2 implement different argument
     * orders and therefore different transformations.
     */

    class Transform_Sub2 : public TypedDistribution<Real> {

    public:
        Transform_Sub2( void );
        virtual ~Transform_Sub2();

        // Basic utility functions
        Transform_Sub2*                                 clone(void) const;                                                                      //!< Clone the object
        static const std::string&                       getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                          getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                     getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                 getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::TransformedDistribution*          createDistribution(void) const;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:
	RevPtr<const RevVariable>                       first = nullptr;
        RevPtr<const RevVariable>                       second_distribution = nullptr;
    };

}

#endif // Transform_Sub2_H
