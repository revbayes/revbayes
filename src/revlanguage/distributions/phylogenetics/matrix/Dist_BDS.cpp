#include "Dist_BDS.h"

#include <cstddef>
#include <string>
#include <vector>

#include "FossilizedBirthDeathRangeProcess.h"

#include "ModelVector.h"
#include "Probability.h"
#include "RealPos.h"
#include "RlString.h"
#include "RlBoolean.h"
#include "RlTaxon.h"
#include "RlMatrixReal.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevVariable.h"
#include "TypeSpec.h"
#include "Taxon.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"
#include "RbVector.h"

namespace RevBayesCore { class DagNode; }

using namespace RevLanguage;

/**
 * Default constructor.
 */
Dist_BDS::Dist_BDS() : FossilizedBirthDeathRangeProcess<MatrixReal>()
{

}


/**
 * Clone the object.
 */
Dist_BDS* Dist_BDS::clone( void ) const
{
    return new Dist_BDS(*this);
}


/**
 * Create a new internal distribution object, with the BDS likelihood always enabled.
 * This is identical to Dist_FBDRP::createDistribution except that use_bds is fixed to true.
 */
RevBayesCore::FossilizedBirthDeathRangeProcess* Dist_BDS::createDistribution( void ) const
{
    // sampling condition
    const std::string& cond  = static_cast<const RlString &>( condition->getRevObject() ).getValue();

    // taxa
    std::vector<RevBayesCore::Taxon> t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    // rates
    RevBayesCore::DagNode* l = lambda->getRevObject().getDagNode();
    RevBayesCore::DagNode* m = mu->getRevObject().getDagNode();
    RevBayesCore::DagNode* p = psi->getRevObject().getDagNode();

    // sampling probability
    RevBayesCore::TypedDagNode<double>* r = static_cast<const Probability &>( rho->getRevObject() ).getDagNode();

    // rate change times
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rt = NULL;
    if ( timeline->getRevObject() != RevNullObject::getInstance() )
    {
        rt = static_cast<const ModelVector<RealPos> &>( timeline->getRevObject() ).getDagNode();
    }

    bool comp = static_cast<const RlBoolean &>( complete->getRevObject() ).getValue();
    std::string c  = comp ? "complete" : static_cast<const RlString &>( reporting->getRevObject() ).getValue();
    bool re = static_cast<const RlBoolean &>( resample->getRevObject() ).getValue();

    // BDS model: use_bds is always true
    RevBayesCore::FossilizedBirthDeathRangeProcess* d = new RevBayesCore::FossilizedBirthDeathRangeProcess(l, m, p, r, rt, cond, t, c, re, true);

    return d;
}


/**
 * Get Rev type of object.
 */
const std::string& Dist_BDS::getClassType( void )
{
    static std::string rev_type = "Dist_BDS";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 */
const TypeSpec& Dist_BDS::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<ModelVector<ModelVector<RealPos> > >::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the distribution (constructor is dn<name>, i.e. dnBDS).
 */
std::string Dist_BDS::getDistributionFunctionName( void ) const
{
    std::string d_name = "BDS";

    return d_name;
}


/**
 * Get type-specification on this object (non-static).
 */
const TypeSpec& Dist_BDS::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}
