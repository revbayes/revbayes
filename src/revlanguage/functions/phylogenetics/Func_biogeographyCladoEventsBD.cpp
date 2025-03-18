#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "BiogeographyCladogeneticBirthDeathFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "Func_biogeographyCladoEventsBD.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbException.h"
#include "RbVector.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_biogeographyCladoEventsBD::Func_biogeographyCladoEventsBD( void ) : TypedFunction<CladogeneticSpeciationRateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_biogeographyCladoEventsBD* Func_biogeographyCladoEventsBD::clone( void ) const {
    
    return new Func_biogeographyCladoEventsBD( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix >* Func_biogeographyCladoEventsBD::createFunction( void ) const
{
    
    // speciation rate vector
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* sr = static_cast<const ModelVector<RealPos> &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // parameterized geography for within-region features (shape == (N))
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* wf = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    // parameterized geography for between-region features (shape == (N,N))
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* bf = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    // support for hidden rate classes
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rm = NULL;
    if (this->args[3].getVariable()->getRevObject() != RevNullObject::getInstance())
    {
        RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rm = static_cast<const ModelVector<RealPos> &>( this->args[3].getVariable()->getRevObject() ).getDagNode();
        if (rm->getValue().size() > 1)
        {
            throw RbException("Only one hidden rate category currently implemented.");
        }
    }
    
    // maximum range size (default == N)
    int n_regions = (int)wf->getValue().size();
    int mrs = 0;
    if (this->args[4].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        mrs = (int)static_cast<const Natural &>( this->args[4].getVariable()->getRevObject() ).getValue();
        if (mrs <= 1 || mrs > n_regions) mrs = n_regions;
        assert( mrs > 1 );
    }
    
    // maximum subrange split size (default == max_range_size)
    int msss = 0;
    if (this->args[5].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        msss = (int)static_cast<const Natural &>( this->args[5].getVariable()->getRevObject() ).getValue();
        if (msss <= 0) msss = mrs;
        assert( msss > 0 );
    }
    
    
    RevBayesCore::BiogeographyCladogeneticBirthDeathFunction* f = new RevBayesCore::BiogeographyCladogeneticBirthDeathFunction( sr, wf, bf, mrs, msss );
    f->setRateMultipliers(rm);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_biogeographyCladoEventsBD::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "speciation_rates", ModelVector<RealPos>::getClassTypeSpec() , "The vector of speciation rates (1: within-region, 2: between-region).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "within_region_features", ModelVector<RealPos>::getClassTypeSpec(), "The within-region feature vector.", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT, NULL ) );
        argumentRules.push_back( new ArgumentRule( "between_region_features", ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "The between-region feature matrix.", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT, NULL ) );
        argumentRules.push_back( new ArgumentRule( "rate_multipliers", ModelVector<RealPos>::getClassTypeSpec() , "The rate multipliers for hidden rate classes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        argumentRules.push_back( new ArgumentRule( "max_range_size",          Natural::getClassTypeSpec(), "The maximum range size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0LL) ) );
        argumentRules.push_back( new ArgumentRule( "max_subrange_split_size", Natural::getClassTypeSpec(), "The maximum size of a daughter subrange following a between-region speciation event.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0LL) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_biogeographyCladoEventsBD::getClassType(void)
{
    
    static std::string rev_type = "Func_biogeographyCladoEventsBD";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_biogeographyCladoEventsBD::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_biogeographyCladoEventsBD::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnBiogeographyCladoEventsBD";
    
    return f_name;
}


const TypeSpec& Func_biogeographyCladoEventsBD::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
