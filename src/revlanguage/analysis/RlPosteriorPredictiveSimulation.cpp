#include "RlPosteriorPredictiveSimulation.h"

#include <cstddef>
#include <string>

#include "ArgumentRules.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "Natural.h"
#include "RlModel.h"
#include "RlString.h"
#include "RlAncestralStateTrace.h"
#include "RlModelTrace.h"
#include "RlUtils.h"
#include "WorkspaceVector.h"
#include "Argument.h"
#include "ArgumentRule.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevNullObject.h"
#include "StringUtilities.h"
#include "Trace.h"
#include "TypeSpec.h"

namespace RevBayesCore { class Model; }

using namespace RevLanguage;

PosteriorPredictiveSimulation::PosteriorPredictiveSimulation() : WorkspaceToCoreWrapperObject<RevBayesCore::PosteriorPredictiveSimulation>()
{
    
    ArgumentRules* runArgRules = new ArgumentRules();
    runArgRules->push_back( new ArgumentRule("thinning", Natural::getClassTypeSpec(), "The number of samples to jump over.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1)) );
    this->methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );
    
}


PosteriorPredictiveSimulation::PosteriorPredictiveSimulation(const RevBayesCore::PosteriorPredictiveSimulation &m) : WorkspaceToCoreWrapperObject<RevBayesCore::PosteriorPredictiveSimulation>( new RevBayesCore::PosteriorPredictiveSimulation( m ) )
{
    
    ArgumentRules* runArgRules = new ArgumentRules();
    runArgRules->push_back( new ArgumentRule("thinning", Natural::getClassTypeSpec(), "The number of samples to jump over.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Natural(1)) );
    this->methods.addFunction( new MemberProcedure( "run", RlUtils::Void, runArgRules) );
    
    
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
PosteriorPredictiveSimulation* RevLanguage::PosteriorPredictiveSimulation::clone(void) const
{
    
    return new PosteriorPredictiveSimulation(*this);
}


void PosteriorPredictiveSimulation::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a new PosteriorPredictiveSimulation object
    const RevBayesCore::Model&                              mdl     = static_cast<const Model &>( model->getRevObject() ).getValue();
    
    RevBayesCore::RbVector<RevBayesCore::ModelTrace>      pt;
    const WorkspaceVector<ModelTrace> &      tmp_pt     = static_cast<const WorkspaceVector<ModelTrace> &>( trace->getRevObject() );
    for ( size_t i=0; i<tmp_pt.size(); ++i)
    {
        const RevBayesCore::ModelTrace &tr = tmp_pt.getElement( i )->getValue();
        pt.push_back( tr );
    }
    
    const std::string &    dir   = static_cast<const RlString &>( directory->getRevObject() ).getValue();
    
    // get vector of ancestral state traces
    if ( ancestral_state_trace->getRevObject() != RevNullObject::getInstance() )
    {
        const WorkspaceVector<AncestralStateTrace>& ast = static_cast<const WorkspaceVector<AncestralStateTrace> &>( ancestral_state_trace->getRevObject() );
        std::vector<RevBayesCore::AncestralStateTrace> ancestral_state_traces;
        for (int i = 0; i < ast.size(); ++i)
        {
            ancestral_state_traces.push_back( ast[i].getValue() );
        }
        value = new RevBayesCore::PosteriorPredictiveSimulation(mdl, dir, pt, ancestral_state_traces);
    }
    else
    {
        value = new RevBayesCore::PosteriorPredictiveSimulation(mdl, dir, pt);
    }
}


/* Map calls to member methods */
RevPtr<RevLanguage::RevVariable> RevLanguage::PosteriorPredictiveSimulation::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "run")
    {
        found = true;
        
        double t            = static_cast<const Natural &>( args[0].getVariable()->getRevObject() ).getValue();
        this->value->run( t );
        
        return NULL;
    }
    
    return RevObject::executeMethod( name, args, found );
}


/** Get Rev type of object */
const std::string& RevLanguage::PosteriorPredictiveSimulation::getClassType(void)
{
    
    //    static std::string rev_type = "PosteriorPredictiveSimulation<" + treeType::getClassType() + ">";
    static std::string rev_type = "PosteriorPredictiveSimulation";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const RevLanguage::TypeSpec& RevLanguage::PosteriorPredictiveSimulation::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::PosteriorPredictiveSimulation>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string PosteriorPredictiveSimulation::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "posteriorPredictiveSimulation";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& PosteriorPredictiveSimulation::getParameterRules(void) const
{
    
    static MemberRules memberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        memberRules.push_back( new ArgumentRule("model", Model::getClassTypeSpec(), "The reference model instance.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("directory", RlString::getClassTypeSpec(), "The name of the directory where we store the simulations.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("trace", WorkspaceVector<ModelTrace>::getClassTypeSpec(), "The sample trace object.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
        memberRules.push_back( new ArgumentRule("ancestralStateTrace", WorkspaceVector<AncestralStateTrace>::getClassTypeSpec(), "The ancestral state trace object. Used only for simulating CDBDP when conditioning on sampled tip states.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        
        rules_set = true;
    }
    
    return memberRules;
}


/** Get type spec */
const TypeSpec& PosteriorPredictiveSimulation::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Get type spec */
void PosteriorPredictiveSimulation::printValue(std::ostream &o) const
{
    
    o << "PosteriorPredictiveSimulation";
}


/** Set a member variable */
void PosteriorPredictiveSimulation::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "model")
    {
        model = var;
    }
    else if ( name == "filetype" )
    {
        filetype = var;
    }
    else if ( name == "directory" )
    {
        directory = var;
    }
    else if ( name == "trace" )
    {
        trace = var;
    }
    else if ( name == "ancestralStateTrace" )
    {
        ancestral_state_trace = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
}
