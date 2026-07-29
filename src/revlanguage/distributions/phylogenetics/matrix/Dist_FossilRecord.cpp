#include "Dist_FossilRecord.h"

#include <string>
#include <vector>

#include "FossilRecordProcess.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlMatrixReal.h"
#include "RlTimeTree.h"
#include "RlTaxon.h"
#include "Taxon.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Dist_FossilRecord::Dist_FossilRecord() : TypedDistribution< ModelVector<Taxon> >()
{
}


Dist_FossilRecord* Dist_FossilRecord::clone( void ) const
{
    return new Dist_FossilRecord(*this);
}


RevBayesCore::FossilRecordProcess* Dist_FossilRecord::createDistribution( void ) const
{
    // the range process stochastic node (a dnFBDRP or dnFBDSP node)
    RevBayesCore::DagNode* rn = ranges->getRevObject().getDagNode();

    // complete=TRUE reports every occurrence, FALSE is first/last; the truncated (exchangeable occurrence) model is reachable only
    // through the deprecated dnFBDRMatrix
    bool comp = static_cast<const RlBoolean &>( complete->getRevObject() ).getValue();

    // the occurrences are read from the range process, so no taxa argument is needed here
    RevBayesCore::FossilRecordProcess* d = new RevBayesCore::FossilRecordProcess( rn, comp );

    return d;
}


const std::string& Dist_FossilRecord::getClassType( void )
{
    static std::string rev_type = "Dist_FossilRecord";
    return rev_type;
}


const TypeSpec& Dist_FossilRecord::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Taxon> >::getClassTypeSpec() ) );
    return rev_type_spec;
}


std::string Dist_FossilRecord::getDistributionFunctionName( void ) const
{
    std::string d_name = "FossilRecord";
    return d_name;
}


const MemberRules& Dist_FossilRecord::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        std::vector<TypeSpec> rangesTypes;
        rangesTypes.push_back( MatrixReal::getClassTypeSpec() );  // dnFBDRP
        rangesTypes.push_back( TimeTree::getClassTypeSpec() );    // dnFBDSP
        dist_member_rules.push_back( new ArgumentRule( "ranges",  rangesTypes, "The FBD range process (a dnFBDRP or dnFBDSP node) supplying b/d, tau, psi and the timeline.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        dist_member_rules.push_back( new ArgumentRule( "complete",  RlBoolean::getClassTypeSpec(),         "Is the fossil record complete (every sampled occurrence reported)? FALSE is first/last.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_FossilRecord::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


void Dist_FossilRecord::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "ranges" )
    {
        ranges = var;
    }
    else if ( name == "complete" )
    {
        complete = var;
    }
    else
    {
        TypedDistribution< ModelVector<Taxon> >::setConstParameter(name, var);
    }
}
