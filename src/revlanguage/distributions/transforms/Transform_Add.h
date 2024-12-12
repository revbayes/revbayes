#ifndef Transform_Add_H
#define Transform_Add_H

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
     *    tnScale(distribution, value)      if isOp == false
     *             or
     *    distribution * value              if isOp == true.
     *
     */

    /*
     * The class is templated to depend on type variable T
     * where T is Real, RealPos, or Probability.
     *
     * The signatures are then
     *
     *    distribution<Real>        * Real      -> distribution<Real>
     *    distribution<RealPos>     * RealPos   -> distribution<RealPos>
     *    distribution<Probability> * Probility -> distribution<Probability>
     *
     */

    template <class T, bool isOp=true>
    class Transform_Add : public TypedDistribution<T> {

    public:
        Transform_Add( void )
        {
            this->markAsTransform();
        }

        virtual ~Transform_Add() {}

        //!< Clone the object
        Transform_Add*                                  clone(void) const
        {
            return new Transform_Add(*this);
        }

        //!< Get Rev type
        static const std::string&                       getClassType(void)
        {
            static std::string rev_type = "Transform_Add";

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
            // create a distribution name variable that is the same for all instance of this class
            if constexpr (isOp)
            {
                return "_add";
            }
            else
            {
                return "shift";
            }
        }

        const MemberRules&                              getParameterRules(void) const;                                                          //!< Get member rules (const)


        // Distribution functions you have to override
        RevBayesCore::TransformedDistribution*          createDistribution(void) const;

    protected:

        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable


    private:
        RevPtr<const RevVariable>                       base_distribution;
	RevPtr<const RevVariable>                       delta;
    };

    template <class T, bool isOp>
    RevBayesCore::TransformedDistribution* Transform_Add<T,isOp>::createDistribution( void ) const
    {
        using namespace Transforms;

        // get the parameters
        const Distribution& rl_vp                      = static_cast<const Distribution &>( base_distribution->getRevObject() );
        RevBayesCore::TypedDistribution<double>* vp    = static_cast<RevBayesCore::TypedDistribution<double>* >( rl_vp.createDistribution() );

        RevBayesCore::TypedDagNode<double>* d           = static_cast<const T &>( delta->getRevObject() ).getDagNode();

        RevBayesCore::TransformedDistribution* dist = new RevBayesCore::TransformedDistribution(*vp, add_transform, add_inverse, log_add_prime, {d});

        delete vp;

        return dist;
    }


     /** Return member rules (no members) */
    template <class T, bool isOp>
    const MemberRules& Transform_Add<T,isOp>::getParameterRules(void) const
    {
        static MemberRules dist_member_rules;
        static bool rules_set = false;

        if ( rules_set == false )
        {
            /* This should be simple but is a bit complicated because
             * the subtype relationships for T don't work on Distribution<T>
             *
             * That is:  RevBayes knows that RealPos < Real.
             *           RevBayes does NOT know that Distribution__RealPos < Distribution__Real.
             *
             * Therefore, depending on T, we need to list the types of the distributions that are allowed.
             * - for T=Real, the distribution can be on Real, RealPos, or Probability.
             * - for T=RealPos, the distribution can be on RealPos or Probability.
             * This encodes the fact that a RealPos plus a RealPos yields a RealPos.
             */
             
            std::vector<TypeSpec> distTypes;

            if (Real::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<Real>::getClassTypeSpec());

            if (RealPos::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<RealPos>::getClassTypeSpec());

            if (Probability::getClassTypeSpec().isDerivedOf(T::getClassTypeSpec()))
                distTypes.push_back(TypedDistribution<Probability>::getClassTypeSpec());

            dist_member_rules.push_back( new ArgumentRule( "baseDistribution", distTypes, "The distribution to be transformed.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
            dist_member_rules.push_back( new ArgumentRule( "delta", T::getClassTypeSpec()   , "The amount added to base random variable.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

            rules_set = true;
        }

        return dist_member_rules;
    }


    /** Set a member variable */
    template <class T, bool isOp>
    void Transform_Add<T,isOp>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
    {
        if ( name == "baseDistribution" )
        {
            base_distribution = var;
        }
        else if ( name == "delta" )
        {
            delta = var;
        }
        else
        {
            TypedDistribution< T >::setConstParameter(name, var);
        }
    }
}

#endif // Transform_Add_H
