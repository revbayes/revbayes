#ifndef Transform_Shift_Pos_H
#define Transform_Shift_Pos_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"
#include "TransformedDistribution.h"

namespace RevLanguage {

    class Transform_Shift_Pos : public TypedDistribution<RealPos> {

    public:
        Transform_Shift_Pos( void );
        virtual ~Transform_Shift_Pos();

        // Basic utility functions
        Transform_Shift_Pos*                                clone(void) const;                                                                      //!< Clone the object
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
        RevPtr<const RevVariable>                       base_distribution;
	RevPtr<const RevVariable>                       delta;
    };

}

#endif // Transform_Shift_Pos_H