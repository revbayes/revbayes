//
//  Func_traitBiogeographyCladoEventsBD.cpp
//  revbayes_traitfig
//
//  Created by Michael Landis on 9/27/24.
//

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "TraitBiogeographyCladogeneticBirthDeathFunction.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "Func_traitBiogeographyCladoEventsBD.h"
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
Func_traitBiogeographyCladoEventsBD::Func_traitBiogeographyCladoEventsBD( void ) : TypedFunction<CladogeneticSpeciationRateMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_traitBiogeographyCladoEventsBD* Func_traitBiogeographyCladoEventsBD::clone( void ) const {
    
    return new Func_traitBiogeographyCladoEventsBD( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticSpeciationRateMatrix >* Func_traitBiogeographyCladoEventsBD::createFunction( void ) const
{
    
    // base rate for within-region speciation
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rho_w = static_cast<const ModelVector<RealPos>& >( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // base rate for between-region speciation
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* rho_b = static_cast<const ModelVector<RealPos>& >( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    // parameterized geography for within-region features (shape == (N))
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* m_w = static_cast<const ModelVector<ModelVector<RealPos> >&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
    
    // parameterized geography for between-region features (shape == (N,N))
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > > >* m_b = static_cast<const ModelVector<ModelVector<ModelVector<RealPos> > >& >( this->args[3].getVariable()->getRevObject() ).getDagNode();
    
    // get dimension information
    
    int num_traits = (int)m_w->getValue().size();
    int num_regions = (int)m_w->getValue()[0].size();
    
    // maximum range size (default == N)
    int mrs = 0;
    if (this->args[4].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        mrs = (int)static_cast<const Natural &>( this->args[4].getVariable()->getRevObject() ).getValue();
        if (mrs <= 1 || mrs > num_regions) mrs = num_regions;
        assert( mrs > 1 );
    }
    
    // maximum subrange split size (default == max_range_size)
    int msss = 0;
    if (this->args[5].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        msss = (int)static_cast<const Natural &>( this->args[5].getVariable()->getRevObject() ).getValue();
        if (msss <= 0) msss = mrs;
        assert( msss > 0 );
    }
    
    // normalize split scores?
    bool nss = false;
    if (this->args[6].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        nss = (bool)static_cast<const RlBoolean &>( this->args[6].getVariable()->getRevObject() ).getValue();
    }
    
    RevBayesCore::TraitBiogeographyCladogeneticBirthDeathFunction* f = new RevBayesCore::TraitBiogeographyCladogeneticBirthDeathFunction( rho_w, rho_b, m_w, m_b, mrs, msss, nss );
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_traitBiogeographyCladoEventsBD::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        argumentRules.push_back( new ArgumentRule( "rho_w", ModelVector<RealPos>::getClassTypeSpec() , "Base speciation rates for within-region speciation events (one per trait).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rho_b", ModelVector<RealPos>::getClassTypeSpec() , "Base speciation rates for between-region speciation events (one per trait).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "m_w", ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "The within-region feature vector (traits x regions).", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT, NULL ) );
        argumentRules.push_back( new ArgumentRule( "m_b", ModelVector<ModelVector<ModelVector<RealPos> > >::getClassTypeSpec(), "The between-region feature matrix (traits x regions x regions).", ArgumentRule::BY_VALUE, ArgumentRule::CONSTANT, NULL ) );
        argumentRules.push_back( new ArgumentRule( "max_range_size",          Natural::getClassTypeSpec(), "The maximum range size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0L) ) );
        argumentRules.push_back( new ArgumentRule( "max_subrange_split_size", Natural::getClassTypeSpec(), "The maximum size of a daughter subrange following a between-region speciation event.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0L) ) );
        argumentRules.push_back( new ArgumentRule( "normalize_split_score", RlBoolean::getClassTypeSpec(), "Normalize the split scores to have geometric mean of 1?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_traitBiogeographyCladoEventsBD::getClassType(void)
{
    
    static std::string rev_type = "Func_traitBiogeographyCladoEventsBD";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_traitBiogeographyCladoEventsBD::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_traitBiogeographyCladoEventsBD::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnTraitBiogeographySpeciationRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_traitBiogeographyCladoEventsBD::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
