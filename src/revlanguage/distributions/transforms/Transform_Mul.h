#ifndef Transform_Mul_H
#define Transform_Mul_H

#include "ModelVector.h"
#include "RealPos.h"
#include "RlTypedDistribution.h"
#include "TypeSpec.h"
#include "TransformedDistribution.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "Real.h"
#include "RlSimplex.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "Transforms.h"

namespace RevLanguage {

    /* This class creates a distribution that implements either
     *
     *    tnShift(distribution, value)      if isOp == false
     *             or
     *    distribution + value              if isOp == true.
     *
     */

    /*
     * The class is templated to depend on type variable T
     * where T is either Real or RealPos.
     *
     * The signatures are then
     *
     *    distribution<Real>    + Real    -> distribution<Real>
     *    distribution<RealPos> + RealPos -> distribution<RealPos>
     *
     */

    template <class T, bool isOp=true>
    class Transform_Mul : public TypedDistribution<T> {

    public:
        Transform_Mul( void )
        {
            this->markAsTransform();
        }

        virtual ~Transform_Mul() {}

        //!< Clone the object
        Transform_Mul* clone(void) const
        {
            return new Transform_Mul(*this);
        }

        //!< Get Rev type
        static const std::string&                       getClassType(void)
        {
            static std::string rev_type = "Transform_Mul";

            return rev_type;
        }

        //!< Get class type spec
        static const TypeSpec&                          getClassTypeSpec(void)
        {
            static TypeSpec rev_type_spec ( TypedDistribution< T >::getClassTypeSpec() );

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
            // create a distribution name variable that is the same for all instance of this class
            if constexpr (isOp)
            {
                return "_mul";
            }
            else
            {
                return "scale";
            }
        }

        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::TransformedDistribution*          createDistribution(void) const;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:
        RevPtr<const RevVariable>                       base_distribution = nullptr;
	RevPtr<const RevVariable>                       lambda = nullptr;
    };


    template <class T, bool isOp>
    RevBayesCore::TransformedDistribution* Transform_Mul<T,isOp>::createDistribution( void ) const
    {
        using namespace Transforms;
        namespace Core = RevBayesCore;

        // get the parameters
        const Distribution& rl_vp               = static_cast<const Distribution &>( base_distribution->getRevObject() );
        Core::TypedDistribution<double>* vp     = static_cast<Core::TypedDistribution<double>* >( rl_vp.createDistribution() );

        Core::TypedDagNode<double>* l           = static_cast<const T &>( lambda->getRevObject() ).getDagNode();

        Core::TransformedDistribution* dist = new Core::TransformedDistribution(*vp, mul_transform, mul_inverse, log_mul_prime, {l});

        delete vp;

        return dist;
    }


    /** Return member rules (no members) */
    template <class T, bool isOp>
    const MemberRules& Transform_Mul<T, isOp>::getParameterRules(void) const
    {
        static MemberRules dist_member_rules;
        static bool rules_set = false;

        if ( rules_set == false )
        {
            /* This should be simple but is a bit complicated because
             * the subtype relationships for T don't work on Distribution<T>
             *
             * That is:  RevBayes knows that Probability < RealPos < Real.
             *           RevBayes does NOT know that Distribution__Probability < Distribution__RealPos < Distribution__Real.
             *
             * Therefore, depending on T, we need to list the types of the distributions that are allowed.
             * - for T=Real, the distribution can be on Real, RealPos, or Probability.
             * - for T=RealPos, the distribution can be on RealPos or Probability.
             * - for T=Probability, the distribution can be on Probability.
             * This encodes the fact that:
             * - a RealPos times a RealPos yields a RealPos.
             * - a Probability times a Probability yields a Probability.
             */
             
            std::vector<TypeSpec> distTypes;

            if (Real::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<Real>::getClassTypeSpec());

            if (RealPos::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<RealPos>::getClassTypeSpec());

            if (Probability::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<Probability>::getClassTypeSpec());

            dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
            dist_member_rules.push_back( new ArgumentRule( "lambda", T::getClassTypeSpec()   , "The amount muled to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

            rules_set = true;
        }

        return dist_member_rules;
    }



    /** Set a member variable */
    template <class T, bool isOp>
    void Transform_Mul<T, isOp>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
    {
        if ( name == "baseDistribution" )
        {
            base_distribution = var;
        }
        else if ( name == "lambda" )
        {
            lambda = var;
        }
        else
        {
            TypedDistribution< T >::setConstParameter(name, var);
        }
    }


}

#endif // Transform_Mul_H
