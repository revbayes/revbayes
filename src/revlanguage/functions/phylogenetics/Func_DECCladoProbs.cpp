//
//  Func_DECCladoProbs.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 1/19/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#include "Func_DECCladoProbs.h"

#include <cstddef>
#include <map>
#include <utility>

#include "CladogeneticProbabilityMatrix.h"
#include "ConstantNode.h"
#include "DECCladogeneticStateFunction.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlDeterministicNode.h"
#include "RlSimplex.h"
#include "RlString.h"
#include "TypedDagNode.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DagNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "Natural.h"
#include "RbBoolean.h"
#include "RbException.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevVariable.h"
#include "RlBoolean.h"
#include "RlConstantNode.h"
#include "RlFunction.h"
#include "Simplex.h"
#include "StringUtilities.h"
#include "TypeSpec.h"
#include "UserFunctionNode.h"

using namespace RevLanguage;

/** default constructor */
Func_DECCladoProbs::Func_DECCladoProbs( void ) : TypedFunction<CladogeneticProbabilityMatrix>( ) {
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
Func_DECCladoProbs* Func_DECCladoProbs::clone( void ) const {
    
    return new Func_DECCladoProbs( *this );
}


RevBayesCore::TypedFunction< RevBayesCore::CladogeneticProbabilityMatrix >* Func_DECCladoProbs::createFunction( void ) const
{
    
    // supplied arguments
    RevBayesCore::TypedDagNode<RevBayesCore::Simplex>* ep = static_cast<const Simplex &>( this->args[0].getVariable()->getRevObject() ).getDagNode();
    unsigned nc = (unsigned)static_cast<const Natural &>( this->args[1].getVariable()->getRevObject() ).getValue();
    unsigned max_num_areas = RevBayesCore::DECCladogeneticStateFunction::MAX_NUM_AREAS;
    if (nc > max_num_areas) {
        std::stringstream ss;
        ss << "fnDECCladoProbs only supports up to " << max_num_areas << " areas\n";
        throw RbException(ss.str());
    }
    
    unsigned mrs = (unsigned)static_cast<const Natural &>( this->args[2].getVariable()->getRevObject() ).getValue();
    if (mrs <= 1) mrs = nc;
    std::string pt = static_cast<const RlString &> ( this->args[3].getVariable()->getRevObject() ).getValue();
    bool ept = (pt == "pattern");
    bool wa = static_cast<const RlBoolean &>( this->args[4].getVariable()->getRevObject() ).getValue();
    ModelVector<RlString> et_tmp = static_cast<const ModelVector<RlString>& >(this->args[5].getVariable()->getRevObject()).getValue();
    
    // default arguments
    unsigned ns = 2;
    
    // check event type vector
    std::map<std::string, std::string> valid_event_names;
    valid_event_names["s"]               = "s";
    valid_event_names["subset_sympatry"] = "s";
    valid_event_names["a"]               = "a";
    valid_event_names["allopatry"]       = "a";
    valid_event_names["f"]               = "f";
    valid_event_names["full_sympatry"]   = "f";
    valid_event_names["j"]               = "j";
    valid_event_names["jump_dispersal"]  = "j";
    
    RevBayesCore::RbVector<std::string> et;
    
    for (size_t i = 0; i < et_tmp.size(); i++)
    {
        std::string s = et_tmp[i];
        std::map<std::string, std::string>::iterator it = valid_event_names.find(s);
        if (it != valid_event_names.end()) {
            et.push_back( it->second );
        }
        else {
            throw RbException() << "\"" << s << "\" is not a valid element for eventTypes." ; 
        }
    }
    
    if (et.size() != ep->getValue().size())
    {
        throw RbException("eventProbs and eventTypes must have equal sizes.");
    }
    
    // connectivity matrix (possible ranges)
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* cg = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[6].getVariable()->getRevObject() ).getDagNode();
    if (cg->getValue().size() > 0 && cg->getValue().size() != nc) {
        throw RbException("Size of connectivity graph does not match number of characters.");
    }
    else if (cg->getValue().size() == 0)
    {
        for (size_t i = 0; i < nc; i++)
        {
            cg->getValue().push_back( RevBayesCore::RbVector<double>(nc, 1.0) );
        }
    }
    
    // vicariance matrix (range-split outcomes)
    RevBayesCore::TypedDagNode< RevBayesCore::RbVector<RevBayesCore::RbVector<double> > >* vg = static_cast<const ModelVector<ModelVector<RealPos> > &>( this->args[7].getVariable()->getRevObject() ).getDagNode();
    bool useVicariance = true;
    if (vg->getValue().size() > 0 && vg->getValue().size() != nc) {
        throw RbException("Size of connectivity graph does not match number of characters.");
    }
    else if (vg->getValue().size() == 0)
    {
        useVicariance = false;
        for (size_t i = 0; i < nc; i++)
        {
            vg->getValue().push_back( RevBayesCore::RbVector<double>(nc, 1.0) );
        }
    }

    // create P matrix
    RevBayesCore::DECCladogeneticStateFunction* f;
    f = new RevBayesCore::DECCladogeneticStateFunction( ep, cg, vg, nc, ns, et, ept, wa, useVicariance, mrs );

    return f;
}


/* Get argument rules */
const ArgumentRules& Func_DECCladoProbs::getArgumentRules( void ) const
{
    
    static ArgumentRules argumentRules = ArgumentRules();
    static bool          rules_set = false;
    
    if ( !rules_set )
    {
        
        argumentRules.push_back( new ArgumentRule( "eventProbs", Simplex::getClassTypeSpec(), "The probabilities of the different event types.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "numCharacters", Natural::getClassTypeSpec(), "The number of characters.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        argumentRules.push_back( new ArgumentRule( "maxRangeSize", Natural::getClassTypeSpec(), "The maximum range size.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(0) ) );

        std::vector<std::string> options;
        options.push_back( "pattern" );
        options.push_back( "class" );
        argumentRules.push_back( new OptionRule( "probType", new RlString("pattern"), options, "Assign event weights over classes of patterns or over specific patterns" ) );
        
        argumentRules.push_back( new ArgumentRule( "widespreadAllopatry", RlBoolean::getClassTypeSpec(), "Allopatry may result in both daughter ranges being larger than size 1.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        
        argumentRules.push_back( new ArgumentRule( "eventTypes",
                                 ModelVector<RlString>::getClassTypeSpec(),
                                 "Vector of cladogenetic event types.",
                                 ArgumentRule::BY_VALUE,
                                 ArgumentRule::ANY) );
        
        argumentRules.push_back( new ArgumentRule( "connectivityGraph",
                                                  ModelVector<ModelVector<RealPos> >::getClassTypeSpec(),
                                                  "Connectivity graph of allowed ranges.",
                                                  ArgumentRule::BY_VALUE,
                                                  ArgumentRule::CONSTANT,
                                                  new ModelVector<ModelVector<RealPos> >() ));
        
        argumentRules.push_back( new ArgumentRule( "vicarianceGraph",
                                                  ModelVector<ModelVector<RealPos> >::getClassTypeSpec(),
                                                  "Graph to model vicariance events.",
                                                  ArgumentRule::BY_VALUE,
                                                  ArgumentRule::CONSTANT,
                                                  new ModelVector<ModelVector<RealPos> >() ));
//                                 options, "Assign event weights over classes of patterns or over specific patterns" ) );
        
        
        rules_set = true;
    }
    
    return argumentRules;
}



const std::string& Func_DECCladoProbs::getClassType(void)
{
    
    static std::string rev_type = "Func_DECCladoProbs";
    
	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& Func_DECCladoProbs::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Function::getClassTypeSpec() ) );
    
	return rev_type_spec;
}


/**
 * Get the primary Rev name for this function.
 */
std::string Func_DECCladoProbs::getFunctionName( void ) const
{
    // create a name variable that is the same for all instance of this class
    std::string f_name = "fnDECCladoProbs";
    
    return f_name;
}

std::vector<std::string> Func_DECCladoProbs::getFunctionNameAliases( void ) const
{
    // create alternative constructor function names variable that is the same for all instance of this class
    std::vector<std::string> a_names;
    a_names.push_back( "fnCladoProbs" );
    
    return a_names;
}


const TypeSpec& Func_DECCladoProbs::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}
