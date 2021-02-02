#ifndef Func_fnNormalizedQuantile_H
#define Func_fnNormalizedQuantile_H

#include "ModelVector.h"
#include "RlTypedFunction.h"

#include <string>
#include "TypedDagNode.h"

namespace RevLanguage {
    
    /** @brief Function to create a normalized vector from a continuous distribution
     *
     * @see RevBayesCore::NormalizeVectorFunction for the created core object
     **/
    template <typename valType>
    class Func_fnNormalizedQuantile : public TypedFunction< ModelVector<valType> > {
        
    public:
        Func_fnNormalizedQuantile();
        
        // Basic utility functions
        Func_fnNormalizedQuantile*                                                  clone(void) const;
        static const std::string&                                                   getClassType(void);
        static const TypeSpec&                                                      getClassTypeSpec(void);
        std::string                                                                 getFunctionName(void) const;
        const TypeSpec&                                                             getTypeSpec(void) const;
        
        // Regular functions
        RevBayesCore::TypedFunction< RevBayesCore::RbVector< double > >*            createFunction(void) const;  //!< Create a new NormalizeVectorFunction object
        const ArgumentRules&                                                        getArgumentRules(void) const;
        
    };
    
}


#include "ArgumentRule.h"
#include "Ellipsis.h"
#include "NormalizeVectorFunction.h"
#include "RlDeterministicNode.h"
#include "QuantileFunction.h"
#include "TypeSpec.h"
#include "VectorFunction.h"
#include "TypedDistribution.h"
#include "ContinuousDistribution.h"
#include "RlContinuousDistribution.h"
#include "RlTypedDistribution.h"

using namespace RevLanguage;

template <typename valType>
Func_fnNormalizedQuantile<valType>::Func_fnNormalizedQuantile() : TypedFunction< ModelVector<valType> >()
{}


template <typename valType>
Func_fnNormalizedQuantile<valType>* Func_fnNormalizedQuantile<valType>::clone( void ) const
{
    return new Func_fnNormalizedQuantile( *this );
}


template <typename valType>
RevBayesCore::TypedFunction< RevBayesCore::RbVector<double> >* Func_fnNormalizedQuantile<valType>::createFunction( void ) const
{
    
    const ContinuousDistribution& rlDistribution    = static_cast<const ContinuousDistribution &>( this->args[0].getVariable()->getRevObject() );
    RevBayesCore::ContinuousDistribution* dist     = static_cast<RevBayesCore::ContinuousDistribution* >( rlDistribution.createDistribution() );    
    
    int nc = (int)static_cast<const Integer &>( this->args[1].getVariable()->getRevObject() ).getValue();
    
    double binWidth = 1.0 / ((double)nc);
    double binMid = binWidth * 0.5;
    std::vector<const RevBayesCore::TypedDagNode<double> *> gRates = std::vector<const RevBayesCore::TypedDagNode<double>* >();
    for (int i=1; i<=nc; i++)
    {
        RevBayesCore::ConstantNode<double> *quant = new RevBayesCore::ConstantNode<double>("", new double((double)i * binWidth - binMid));
        RevBayesCore::DeterministicNode<double> *tmp = new RevBayesCore::DeterministicNode<double>("", new RevBayesCore::QuantileFunction(quant,  dist));
		gRates.push_back( tmp );
    }
    RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> > *discRates = new RevBayesCore::DeterministicNode< RevBayesCore::RbVector<double> >( "", new RevBayesCore::VectorFunction<double>(gRates) );
    RevBayesCore::NormalizeVectorFunction *func = new RevBayesCore::NormalizeVectorFunction( discRates, new RevBayesCore::ConstantNode<double>("", new double(1.0) ) );
    
    return func;
}


template <typename valType>
const ArgumentRules& Func_fnNormalizedQuantile<valType>::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "contDistribution", TypedDistribution<valType>::getClassTypeSpec(), "The distribution which we discretize.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCategories", Integer::getClassTypeSpec()                      , "How many discrete categories?", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        rules_set = true;
    }
    
    return argumentRules;
}


template <typename valType>
const std::string& Func_fnNormalizedQuantile<valType>::getClassType(void)
{
    
    static std::string revClassType = "Func_fnNormalizedQuantile";
    
	return revClassType;
}


template <typename valType>
const RevLanguage::TypeSpec& Func_fnNormalizedQuantile<valType>::getClassTypeSpec(void)
{
    
    static TypeSpec revClassType = TypeSpec( Func_fnNormalizedQuantile<valType>::getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return revClassType;
}


template <typename valType>
std::string RevLanguage::Func_fnNormalizedQuantile<valType>::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnNormalizedQuantile";
    
    return f_name;
}


template <typename valType>
const RevLanguage::TypeSpec& Func_fnNormalizedQuantile<valType>::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


#endif


