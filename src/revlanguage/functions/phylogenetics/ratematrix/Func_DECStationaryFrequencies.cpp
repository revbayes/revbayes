//
//  Func_DECStationaryFrequencies.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/5/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "ConstantNode.h"
#include "DECStationaryFrequenciesFunction.h"
#include "DECRateMatrixFunction.h"
#include "Func_DECStationaryFrequencies.h"
#include "Integer.h"
#include "Natural.h"
#include "OptionRule.h"
#include "RbConstants.h"
#include "Real.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlDeterministicNode.h"
#include "RlMatrixReal.h"
#include "RlRateMatrix.h"
#include "RlSimplex.h"
#include "TypedDagNode.h"

using namespace RevLanguage;

/** default constructor */
Func_DECStationaryFrequencies::Func_DECStationaryFrequencies( void ) : TypedFunction<RateMatrix>( )
{
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_DECStationaryFrequencies* Func_DECStationaryFrequencies::clone( void ) const
{
    
    return new Func_DECStationaryFrequencies( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::RateGenerator >* Func_DECStationaryFrequencies::createFunction( void ) const
{
    
    // dispersal rates
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* dr;
    dr = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    
    // extirpation rates
    RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >* er;
    er = static_cast<const ModelVector<RealPos> &>( this->args[1].getVariable()->getRevObject() ).getDagNode();
    
    size_t num_statesEr = er->getValue().size();
    std::vector<size_t> num_statesDr;
    for (size_t i = 0; i < dr->getValue().size(); i++)
    {
        num_statesDr.push_back( dr->getValue()[i].size() );
        if (num_statesDr[i] != num_statesEr)
        {
            throw RbException("The dimension between dispersal and extirpation rates does not match.");
        }
        for (size_t j = 0; j < i; j++)
        {
            if (num_statesDr[i] != num_statesDr[j])
            {
                throw RbException("The dispersal matrix is not square.");
            }
        }
    }
    if (dr->getValue().size() != num_statesEr)
    {
        throw RbException("The dimension between dispersal and extirpation rates does not match.");
    }
    
    // range size probabilities
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex >* rs = NULL;
    if ( this->args[2].getVariable() != NULL && this->args[2].getVariable()->getRevObject() != RevNullObject::getInstance()) {
        
        rs = static_cast<const Simplex&>( this->args[2].getVariable()->getRevObject() ).getDagNode();
        if (rs->getValue().size() != num_statesEr && rs->getValue().size() > 0) {
            throw RbException("The probabilities of range sizes must equal the number of areas.");
        }
        else {
            size_t n = num_statesEr;
            double p = 1.0 / n;
            rs = new RevBayesCore::ConstantNode<RevBayesCore::Simplex >("", new RevBayesCore::Simplex(n,p));
        }
    }
    else {
        size_t n = num_statesEr+1;
        double p = 1.0 / n;
        rs = new RevBayesCore::ConstantNode<RevBayesCore::Simplex >("", new RevBayesCore::Simplex(n,p));
    }
    
    std::string nullRangeStr = static_cast<const RlString &>( this->args[3].getVariable()->getRevObject() ).getValue();
//    bool cs = nullRangeStr=="CondSurv";
//    bool ex = nullRangeStr=="Exclude";
    //    std::cout << nullRangeStr << " " << cs << " " << ex << "\n";
    
//    bool os = static_cast<const RlBoolean&>(this->args[4].getVariable()->getRevObject() ).getValue();
    
//    bool uc = false;
    RevBayesCore::DECRateMatrixFunction* f = NULL;// = new RevBayesCore::DECRateMatrixFunction( dr, er, rs, cs, ex, os, uc );
//    RevBayesCore::DECStationaryFrequenciesFunction* f; // = new RevBayesCore::DECStationaryFrequenciesFunction(dr, er, rs, cs, ex, os, uc);
    
    return f;
}


/* Get argument rules */
const ArgumentRules& Func_DECStationaryFrequencies::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "dispersalRates"  , ModelVector<ModelVector<RealPos> >::getClassTypeSpec(), "Matrix of dispersal rates between areas.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "extirpationRates", ModelVector<RealPos>::getClassTypeSpec(), "The per area extinction rates.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "rangeSize",        Simplex::getClassTypeSpec(), "Relative proportions of range sizes.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new Simplex( RevBayesCore::RbVector<double>() ) ) );
        argumentRules.push_back(  new ArgumentRule( "birthRate",        RealPos::getClassTypeSpec(), "Birth rate.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0)) );
        argumentRules.push_back( new ArgumentRule( "cladoProbs",  MatrixReal::getClassTypeSpec(), "Cladogenetic probabilities.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::DETERMINISTIC ) );
        
        std::vector<std::string> options;
        options.push_back( "CondSurv" );
        options.push_back( "Exclude" );
//        options.push_back( "Include" );
        argumentRules.push_back( new OptionRule( "nullRange", new RlString("CondSurv"), options, "How should DEC handle the null range?" ) );
        
        argumentRules.push_back( new ArgumentRule( "orderStatesBySize", RlBoolean::getClassTypeSpec(), "Order states by size?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ));
        
        
        rules_set = true;
    }
    
    return argumentRules;
}


const std::string& Func_DECStationaryFrequencies::getClassType(void)
{
    
    static std::string rev_type = "Func_DECStationaryFrequencies";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_DECStationaryFrequencies::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_DECStationaryFrequencies::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDECStationaryFrequencies";
    
    return f_name;
}


const TypeSpec& Func_DECStationaryFrequencies::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
