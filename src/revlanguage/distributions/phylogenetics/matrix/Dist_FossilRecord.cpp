#include "Dist_FossilRecord.h"

#include <string>
#include <vector>

#include "FossilRecordProcess.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "RlString.h"
#include "RlMatrixReal.h"
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
    // the skeleton stochastic node (a MatrixReal-valued dnFBDRP node)
    RevBayesCore::DagNode* sk = skeleton->getRevObject().getDagNode();

    // the reporting model (complete | firstlast | uniform), pushed onto the skeleton
    const std::string& rep = static_cast<const RlString &>( reporting->getRevObject() ).getValue();

    // the observed occurrences (data)
    const std::vector<RevBayesCore::Taxon>& t = static_cast<const ModelVector<Taxon> &>( taxa->getRevObject() ).getValue();

    RevBayesCore::FossilRecordProcess* d = new RevBayesCore::FossilRecordProcess( sk, rep, t );

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
        dist_member_rules.push_back( new ArgumentRule( "skeleton",  MatrixReal::getClassTypeSpec(),        "The FBD-range skeleton (a dnFBDRP node) supplying b/d, tau, psi and the timeline.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        dist_member_rules.push_back( new ArgumentRule( "reporting", RlString::getClassTypeSpec(),          "Reporting model: complete | firstlast | uniform.",                                 ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("uniform") ) );
        dist_member_rules.push_back( new ArgumentRule( "taxa",      ModelVector<Taxon>::getClassTypeSpec(),"The taxa with their fossil occurrences (the observed record).",                     ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );

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
    if ( name == "skeleton" )
    {
        skeleton = var;
    }
    else if ( name == "reporting" )
    {
        reporting = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else
    {
        TypedDistribution< ModelVector<Taxon> >::setConstParameter(name, var);
    }
}
