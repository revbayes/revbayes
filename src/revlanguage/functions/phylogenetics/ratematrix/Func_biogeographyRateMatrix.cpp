//
//  Func_biogeographyRateMatrix.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 3/16/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//


#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

#include "ConstantNode.h"
#include "BiogeographyRateMatrixFunction.h"
#include "Func_biogeographyRateMatrix.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlDeterministicNode.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RateGenerator.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevNullObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlFunction.h"
#include "RlString.h"
#include "RlTypedFunction.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "TypedFunction.h"

using namespace RevLanguage;

/** default constructor */
Func_biogeographyRateMatrix::Func_biogeographyRateMatrix( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_biogeographyRateMatrix* Func_biogeographyRateMatrix::clone( void ) const
{
    
    return new Func_biogeographyRateMatrix( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_biogeographyRateMatrix::createFunction( void ) const
{
    
    // dispersal rates
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* dr;
    dr = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // extirpation rates
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* er;
    er = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    size_t num_statesEr = er->getValue().size();
    std::vector<size_t> num_statesDr;
    RevBayesCore::RbVector<RevBayesCore::RbVector<double> > dr_tmp = dr->getValue();
    for (size_t i = 0; i < dr_tmp.size(); i++)
    {
        num_statesDr.push_back( dr_tmp[i].size() );
        if (num_statesDr[i] != num_statesEr)
        {
            throw RbException("The dimension between dispersal and extirpation rates does not match.");
        }
        if (i > 0)
        {
            if (num_statesDr[i] != num_statesDr[i-1])
            {
                throw RbException("The dispersal matrix is not square.");
            }
        }
    }
    if (dr->getValue().size() != num_statesEr)
    {
        throw RbException("The dimensions of the dispersal and extirpation rates does not match.");
    }

    // maximum range size
    size_t mrs = static_cast<const Natural&>(this->args[2].getVariable()->getRevObject() ).getValue();
    
    if (mrs < 1 || mrs > er->getValue().size())
    {
        mrs = er->getValue().size();
    }
    RevBayesCore::BiogeographyRateMatrixFunction* f = new RevBayesCore::BiogeographyRateMatrixFunction( dr, er, mrs );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_biogeographyRateMatrix::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "dispersalRates"  , ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "Matrix of dispersal rates between areas.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "extirpationRates", ModelVector<RealPos>::getClassTypeSpec(), "Matrix of extirpation rates.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "maxRangeSize", Natural::getClassTypeSpec(), "Maximum range size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0LL) ));
    
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_biogeographyRateMatrix::getClassType(void)
{
    
    static std::string rev_type = "Func_biogeographyRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_biogeographyRateMatrix::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_biogeographyRateMatrix::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnBiogeographyRateMatrix";
    
    return f_name;
}


const TypeSpec& Func_biogeographyRateMatrix::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
