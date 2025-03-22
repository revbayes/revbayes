#include <iosfwd>
#include <vector>

#include "Argument.h"
#include "ArgumentRule.h"
#include "Func_PhylogeneticIndependentContrastsMultiSample.h"
#include "ModelVector.h"
#include "PhylogeneticIndependentContrastsMultiSampleFunction.h"
#include "RlContinuousCharacterData.h"
#include "RlTaxon.h"
#include "RlTimeTree.h"
#include "TypeSpec.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbBoolean.h"
#include "RbVector.h"
#include "Real.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlFunction.h"
#include "RlTypedFunction.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevBayesCore { class ContinuousCharacterData; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevLanguage;

/** Default constructor */
Func_PhylogeneticIndependentContrastsMultiSample::Func_PhylogeneticIndependentContrastsMultiSample( void ) : TypedFunction<ModelVector<Real> >()
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_PhylogeneticIndependentContrastsMultiSample* Func_PhylogeneticIndependentContrastsMultiSample::clone( void ) const
{
    
    return new Func_PhylogeneticIndependentContrastsMultiSample( *this );
}


RevBayesCore::TypedFunction<RevBayesCore::RbVector<double> >* Func_PhylogeneticIndependentContrastsMultiSample::createFunction( void ) const
{
    const RevBayesCore::TypedDagNode<RevBayesCore::Tree>* tau = static_cast<const TimeTree &>( args[0].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<RevBayesCore::ContinuousCharacterData>* data = static_cast<const ContinuousCharacterData &>( args[1].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::TypedDagNode<std::int64_t>* site = static_cast<const Natural &>( args[2].getVariable()->getRevObject() ).getDagNode();
    const RevBayesCore::RbVector<RevBayesCore::Taxon>& taxa = static_cast<const ModelVector<Taxon> &>( args[3].getVariable()->getRevObject() ).getValue();
    const bool& normalized = static_cast<const RlBoolean &>( args[4].getVariable()->getRevObject() ).getValue();
    
    RevBayesCore::PhylogeneticIndependentContrastsMultiSampleFunction* f = new RevBayesCore::PhylogeneticIndependentContrastsMultiSampleFunction( tau, data, site, taxa, normalized );
    
    return f;
}


/** Get argument rules */
const ArgumentRules& Func_PhylogeneticIndependentContrastsMultiSample::getArgumentRules( void ) const
{
    
    static ArgumentRules argument_rules = ArgumentRules();
    static bool rules_set = false;
    
    if ( rules_set == false )
    {
        
        argument_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "data", ContinuousCharacterData::getClassTypeSpec(), "The character data object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "site", Natural::getClassTypeSpec(), "The site for which we compute the contrasts.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The vector of taxa which have species and individual names.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argument_rules.push_back( new ArgumentRule( "normalized", RlBoolean::getClassTypeSpec(), "If we should normalize the contrasts by their standard deviations.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        
        rules_set = true;
    }
    
    return argument_rules;
}


/** Get Rev type of object */
const std::string& Func_PhylogeneticIndependentContrastsMultiSample::getClassType(void)
{
    
    static std::string rev_type = "Func_PhylogeneticIndependentContrastsMultiSample";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& Func_PhylogeneticIndependentContrastsMultiSample::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_PhylogeneticIndependentContrastsMultiSample::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnPhylogeneticIndependentContrastsMultiSample";
    
    return f_name;
}

std::vector<std::string> Func_PhylogeneticIndependentContrastsMultiSample::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnPIC" );
    
    return a_names;
}

/** Get type spec */
const TypeSpec& Func_PhylogeneticIndependentContrastsMultiSample::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}



