#ifndef RlFossilizedBirthDeathRangeProcess_H
#define RlFossilizedBirthDeathRangeProcess_H

#include "ModelVector.h"
#include "OptionRule.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the FossilizedBirthDeathRangeProcess
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    template <typename rlType>
    class FossilizedBirthDeathRangeProcess : public TypedDistribution<rlType> {
        
    public:
        virtual                                             ~FossilizedBirthDeathRangeProcess(void);                                          //!< Destructor
        
        // Basic utility functions you have to overwrite
        virtual FossilizedBirthDeathRangeProcess<rlType>*   clone(void) const = 0;                                                              //!< Clone the object

        // Basic utility functions you may want to overwrite
        const MemberRules&                                  getParameterRules(void) const;                                                      //!< Get member rules (const): skeleton + reporting args (the fused process)

        // Basic utility functions
        static const std::string&                           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                             //!< Get class type spec

    protected:
        FossilizedBirthDeathRangeProcess( void );

        const MemberRules&                                  getSkeletonParameterRules(void) const;                                              //!< Member rules WITHOUT the reporting args; the b/d skeleton (dnFBDRP) takes its reporting model from a downstream dnFossilRecord node instead.
        static void                                         appendParameterRules(MemberRules &rules, bool include_reporting);                   //!< Build the rule list; the reporting args sit at their fused-process position so the facade's argument order is unchanged.

        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable
    
        // members        
        RevPtr<const RevVariable>                           lambda;                                                                             //!< The speciation rate(s)
        RevPtr<const RevVariable>                           mu;                                                                                 //!< The extinction rate(s)
        RevPtr<const RevVariable>                           psi;                                                                                //!< The fossilization rate(s)
        RevPtr<const RevVariable>                           rho;                                                                                //!< The extant sampling proportion
        RevPtr<const RevVariable>                           timeline;                                                                           //!< The interval times
        RevPtr<const RevVariable>                           taxa;                                                                               //!< The taxa
        RevPtr<const RevVariable>                           condition;                                                                          //!< The condition of the process
        RevPtr<const RevVariable>                           complete;                                                                           //!< Is the fossil record complete?
        RevPtr<const RevVariable>                           reporting;                                                                          //!< Reporting model when the record is incomplete
        RevPtr<const RevVariable>                           resample;

    };
    
}


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
template <typename rlType>
RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::FossilizedBirthDeathRangeProcess() : TypedDistribution<rlType>()
{

}


/**
 * Default destructor.
 *
 * The default destructor does nothing.
 */
template <typename rlType>
RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::~FossilizedBirthDeathRangeProcess()
{

}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
template <typename rlType>
const std::string& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getClassType(void)
{

    static std::string rev_type = "FossilizedBirthDeathRangeProcess";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
template <typename rlType>
const TypeSpec& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<rlType>::getClassTypeSpec() ) );

    return rev_type_spec;
}



/**
 * Build the member rules shared by the matrix and tree range processes.
 *
 * The reporting args (complete, reporting) describe how sampled occurrences make it into the
 * observed record. They belong to the fossil-record term, so the bare skeleton (dnFBDRP) omits
 * them and picks its reporting model up from a downstream dnFossilRecord node. They are appended
 * in place rather than at the end so that including them reproduces the historical argument
 * order of the fused process exactly.
 *
 * \param[in,out] rules              The rule list to append to.
 * \param[in]     include_reporting  Include the reporting args (fused process) or omit them (skeleton).
 */
template <typename rlType>
void RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::appendParameterRules(MemberRules &rules, bool include_reporting)
{
    std::vector<TypeSpec> paramTypes;
    paramTypes.push_back( RealPos::getClassTypeSpec() );
    paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
    rules.push_back( new ArgumentRule( "lambda",  paramTypes, "The speciation rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    rules.push_back( new ArgumentRule( "mu",      paramTypes, "The extinction rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
    rules.push_back( new ArgumentRule( "psi",     paramTypes, "The fossil sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
    rules.push_back( new ArgumentRule( "rho",     Probability::getClassTypeSpec(), "The extant sampling fraction.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

    rules.push_back( new ArgumentRule( "timeline",   ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

    std::vector<std::string> optionsCondition;
    optionsCondition.push_back( "time" );
    optionsCondition.push_back( "sampling" );
    optionsCondition.push_back( "survival" );
    rules.push_back( new OptionRule( "condition", new RlString("time"), optionsCondition, "The condition of the process." ) );
    rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The taxa with fossil occurrence information.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

    if ( include_reporting )
    {
        rules.push_back( new ArgumentRule( "complete", RlBoolean::getClassTypeSpec(), "Is the fossil record complete (every sampled occurrence reported)?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        std::vector<std::string> optionsReporting;
        optionsReporting.push_back( "firstlast" );
        optionsReporting.push_back( "uniform" );
        rules.push_back( new OptionRule( "reporting", new RlString("firstlast"), optionsReporting, "Reporting model for an incomplete record (used when complete=FALSE): firstlast (oldest and youngest occurrence) or uniform (exchangeable, capped at the max observed count)." ) );
    }

    rules.push_back( new ArgumentRule( "resample", RlBoolean::getClassTypeSpec(), "Resample augmented ages?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(true) ) );
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the fossilized birth-death range process are:
 * (1) time of the process since the origin.
 * (2) time of the process since the rootAge.
 * (3) the sampling probability.
 * (4) the sampling strategy.
 * (5) the condition.
 * (6) the number of taxa.
 * (7) the taxon names.
 * (8) the clade constraints.
 *
 * \return The member rules.
 */
template <typename rlType>
const MemberRules& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        appendParameterRules( memberRules, true );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Get the member rules of the bare b/d skeleton: everything except the reporting args.
 *
 * \return The member rules.
 */
template <typename rlType>
const MemberRules& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getSkeletonParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        appendParameterRules( memberRules, false );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
template <typename rlType>
void RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "psi" )
    {
        psi = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "timeline" )
    {
        timeline = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "complete" )
    {
        complete = var;
    }
    else if ( name == "reporting" )
    {
        reporting = var;
    }
    else if ( name == "resample" )
    {
        resample = var;
    }
    else
    {
        TypedDistribution<rlType>::setConstParameter(name,var);
    }

}

#endif
