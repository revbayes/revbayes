#include "Dist_FBDRMatrix.h"

#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "FossilizedBirthDeathRangeProcess.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MethodTable.h"
#include "StochasticNode.h"
#include "RlStochasticNode.h"
#include "DistributionMemberFunction.h"
#include "RlDistributionMemberFunction.h"
#include "ModelVector.h"
#include "Natural.h"
#include "OptionRule.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RbException.h"
#include "RlUserInterface.h"
#include "DagMemberFunction.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "MatrixReal.h"
#include "ModelObject.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlAbstractHomologousDiscreteCharacterData.h"
#include "RlBoolean.h"
#include "RlDagMemberFunction.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "RlStochasticNode.h"
#include "RlDistribution.h"
#include "RlTypedDistribution.h"
#include "RlTypedFunction.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "TypedFunction.h"

namespace RevBayesCore { class DagNode; }

using namespace RevLanguage;

/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
Dist_FBDRMatrix::Dist_FBDRMatrix() : FossilizedBirthDeathRangeProcess<MatrixReal>()
{

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Dist_FBDRMatrix* Dist_FBDRMatrix::clone( void ) const
{
    return new Dist_FBDRMatrix(*this);
}


/**
 * Create a new internal distribution object.
 *
 * This function simply dynamically allocates a new internal distribution object that can be
 * associated with the variable. The internal distribution object is created by calling its
 * constructor and passing the distribution-parameters (other DAG nodes) as arguments of the
 * constructor. The distribution constructor takes care of the proper hook-ups.
 *
 * \return A new internal distribution object.
 */
RevBayesCore::FossilizedBirthDeathRangeProcess* Dist_FBDRMatrix::createDistribution( void ) const
{
    static bool warned = false;
    if ( warned == false )
    {
        RBOUT("\nWarning! `dnFBDRMatrix` is deprecated. It fuses the birth-death range process with the");
        RBOUT("         fossil record, and takes the occurrences as an argument instead of as clamped data.");
        RBOUT("         Use `dnFBDRP` for the range process and `dnFossilRecord` for the record instead.");
        RBOUT("         See `?dnFossilRecord` for an example.\n");
        warned = true;
    }

    // get the parameters

    // sampling condition
    const std::string& cond  = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // get the taxa to simulate either from a vector of rev taxon objects or a vector of names
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // speciation rate
    RevBayesCore::DagNode* l = lambda->getRevObject().getDagNode();
    // extinction rate
    RevBayesCore::DagNode* m = mu->getRevObject().getDagNode();
    // fossilization rate
    RevBayesCore::DagNode* p = psi->getRevObject().getDagNode();

    // sampling probability
    RevBayesCore::TypedDagNode<double>* r = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();

    // rate change times
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
        rt = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    }

    // complete=TRUE reports every occurrence; FALSE is first/last, or the truncated model when truncated is given
    bool comp = static_cast<const RlBoolean &>( complete->getRevObject() ).getValue();

    // a cap selects the truncated model: a taxon reporting K may have unreported occurrences
    size_t K = 0;
    if ( truncated->getRevObject() != RevNullObject::getInstance() )
    {
        if ( comp == true )
        {
            RBOUT( "Warning: \"truncated\" is ignored when complete=TRUE." );
        }
        else
        {
            K = size_t( static_cast<const Natural &>( truncated->getRevObject() ).getValue() );

            if ( K == 0 )
            {
                throw(RbException("The truncated (exchangeable occurrence) reporting model requires a reporting cap of at least 1."));
            }
        }
    }

    bool re = static_cast<const RlBoolean &>( resample->getRevObject() ).getValue();

    // optional origin time of the process
    RevBayesCore::TypedDistribution<double>* op = createOriginPrior();

    // report_internally = true: the fused facade adds the fossil-record term inline
    RevBayesCore::FossilizedBirthDeathRangeProcess* d = new RevBayesCore::FossilizedBirthDeathRangeProcess(l, m, p, r, rt, cond, t, comp, K, re, NULL, op, true);

    return d;
}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
const std::string& Dist_FBDRMatrix::getClassType( void )
{

    static std::string rev_type = "Dist_FBDRMatrix";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
const TypeSpec& Dist_FBDRMatrix::getClassTypeSpec( void )
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<ModelVector<ModelVector<RealPos> > >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution.
 * This name is used for the constructor and the distribution functions,
 * such as the density and random value function
 *
 * The canonical name FossilizedBirthDeathRange now belongs to the range process (dnFBDRP), so the
 * deprecated fused form registers only under dnFBDRMatrix -- the name the fbd_range tutorials
 * use on this branch. It is not given a long name of its own, since it is on its way out.
 *
 * \return Rev name of constructor function.
 */
std::string Dist_FBDRMatrix::getDistributionFunctionName( void ) const
{
    // create a distribution name variable that is the same for all instance of this class
    std::string d_name = "FBDRMatrix";

    return d_name;
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * The member rules of the fossilized birth-death process are:
 * (1) the speciation rate lambda which must be a positive real.
 * (2) the extinction rate mu that must be a positive real.
 * (3) the fossil sampling rate psi that must be a positive real.
 * (4) the extant sampling rate rho that must be a positive real.
 *
 * \return The member rules.
 */
const MemberRules& Dist_FBDRMatrix::getParameterRules(void) const
{

    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        dist_member_rules.push_back( originPriorRule() );

        // add the rules from the base class, including the reporting args
        const MemberRules &parentRules = FossilizedBirthDeathRangeProcess<MatrixReal>::getParameterRules();
        dist_member_rules.insert(dist_member_rules.end(), parentRules.begin(), parentRules.end());

        rules_set = true;
    }

    return dist_member_rules;
}


/**
 * Get type-specification on this object (non-static).
 *
 * \return The type spec of this object.
 */
const TypeSpec& Dist_FBDRMatrix::getTypeSpec( void ) const
{

    static TypeSpec ts = getClassTypeSpec();

    return ts;
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
void Dist_FBDRMatrix::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

           FossilizedBirthDeathRangeProcess<MatrixReal>::setConstParameter(name,var);

}


/**
 * The augmented first (tau_1) and last (tau_K) occurrence ages. They are internal to the
 * distribution, so a deterministic node is the only way a monitor can reach them.
 */
RevLanguage::MethodTable Dist_FBDRMatrix::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<MatrixReal>::getDistributionMethods();

    ArgumentRules* first_ages_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDRMatrix, ModelVector<RealPos> >( "getAugmentedFirstAges", variable, first_ages_arg_rules, true ) );

    ArgumentRules* last_ages_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDRMatrix, ModelVector<RealPos> >( "getAugmentedLastAges", variable, last_ages_arg_rules, true ) );

    ArgumentRules* origin_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDRMatrix, RealPos >( "getOrigin", variable, origin_arg_rules, true ) );

    ArgumentRules* birth_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDRMatrix, ModelVector<RealPos> >( "getBirthAges", variable, birth_arg_rules, true ) );

    ArgumentRules* death_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_FBDRMatrix, ModelVector<RealPos> >( "getDeathAges", variable, death_arg_rules, true ) );

    return methods;
}
