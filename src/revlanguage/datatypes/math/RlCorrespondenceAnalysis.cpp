
#include <cstddef>
#include <vector>
#include <string>
#include <ostream>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "CorrespondenceAnalysis.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlCorrespondenceAnalysis.h"
#include "RlMatrixReal.h"
#include "TypeSpec.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MatrixReal.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RbVector.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"
#include "WorkspaceToCoreWrapperObject.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

CorrespondenceAnalysis::CorrespondenceAnalysis() : WorkspaceToCoreWrapperObject<RevBayesCore::CorrespondenceAnalysis>()
{
    
    // Add methods
    ArgumentRules* columnWeightsArgRules    = new ArgumentRules();
    ArgumentRules* principalAxesArgRules    = new ArgumentRules();
    ArgumentRules* rowCoordArgRules         = new ArgumentRules();
    ArgumentRules* colCoordArgRules         = new ArgumentRules();
    ArgumentRules* rankArgRules             = new ArgumentRules();
    
    methods.addFunction( new MemberProcedure( "columnWeights", ModelVector<RealPos>::getClassTypeSpec(), columnWeightsArgRules ) );
    methods.addFunction( new MemberProcedure( "principalAxes", MatrixReal::getClassTypeSpec(), principalAxesArgRules ) );
    methods.addFunction( new MemberProcedure( "rowCoordinates", MatrixReal::getClassTypeSpec(), rowCoordArgRules ) );
    methods.addFunction( new MemberProcedure( "columnCoordinates", MatrixReal::getClassTypeSpec(), colCoordArgRules ) );
    methods.addFunction( new MemberProcedure( "rank", Natural::getClassTypeSpec(), rankArgRules ) );

}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'b'.
 *
 * \return A new copy of the process.
 */
CorrespondenceAnalysis* CorrespondenceAnalysis::clone(void) const
{
    
    return new CorrespondenceAnalysis(*this);
}


void CorrespondenceAnalysis::constructInternalObject( void )
{
    // we free the memory first
    delete value;
    
    // now allocate a CorrespondenceAnalysis
    RevBayesCore::MatrixReal d  = static_cast<const MatrixReal &>( data->getRevObject() ).getValue();
    unsigned int na             = (unsigned int)static_cast<const Natural &>( num_axes->getRevObject() ).getValue();
    double tol                  = static_cast<const RealPos &>( tolerance->getRevObject() ).getValue();
    
    value = new RevBayesCore::CorrespondenceAnalysis( d, na, tol );
    
}

/* Map calls to member methods */
RevPtr<RevVariable> CorrespondenceAnalysis::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if ( name == "columnWeights" )
    {
        found = true;
        
        const std::vector<double> &w = value->getColumnWeights();
        
        return new RevVariable( new ModelVector<RealPos>( w ) );
    }
    else if ( name == "principalAxes" )
    {
        found = true;
        
        const RevBayesCore::MatrixReal &m = value->getPrincipalAxes();
        
        return new RevVariable( new MatrixReal( m ) );
    }
    else if ( name == "rowCoordinates" )
    {
        found = true;
        
        const RevBayesCore::MatrixReal &m = value->getRowCoordinates();
        
        return new RevVariable( new MatrixReal( m ) );
    }
    else if ( name == "columnCoordinates" )
    {
        found = true;
        
        const RevBayesCore::MatrixReal &m = value->getColCoordinates();
        
        return new RevVariable( new MatrixReal( m ) );
    }
    else if ( name == "rank" )
    {
        found = true;
        
        size_t n = value->getNbOfKeptAxes();
        
        return new RevVariable( new Natural( n ) );
    }
    
    return RevObject::executeMethod( name, args, found );
}

/** Get Rev type of object */
const std::string& CorrespondenceAnalysis::getClassType(void)
{
    
    static std::string rev_type = "CorrespondenceAnalysis";
    
    return rev_type;
}

/** Get class type spec describing type of object */
const TypeSpec& CorrespondenceAnalysis::getClassTypeSpec(void)
{
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( WorkspaceToCoreWrapperObject<RevBayesCore::CorrespondenceAnalysis>::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 *
 * \return Rev name of constructor function.
 */
std::string CorrespondenceAnalysis::getConstructorFunctionName( void ) const
{
    // create a constructor function name variable that is the same for all instance of this class
    std::string c_name = "CorrespondenceAnalysis";
    
    return c_name;
}


/** Return member rules (no members) */
const MemberRules& CorrespondenceAnalysis::getParameterRules(void) const
{
    
    static MemberRules CorrespondenceAnalysisMemberRules;
    static bool rules_set = false;
    
    if ( !rules_set )
    {
        
        CorrespondenceAnalysisMemberRules.push_back( new ArgumentRule("data",       MatrixReal::getClassTypeSpec(), "The matrix of numerical values.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        CorrespondenceAnalysisMemberRules.push_back( new ArgumentRule("numAxes",    Natural::getClassTypeSpec(),    "The number of principle components.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
        CorrespondenceAnalysisMemberRules.push_back( new ArgumentRule("tolerance",  RealPos::getClassTypeSpec(),    "The allowed machine tolerance.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0000001) ) );
        
        rules_set = true;
    }
    
    return CorrespondenceAnalysisMemberRules;
}


/** Get type spec */
const TypeSpec& CorrespondenceAnalysis::getTypeSpec( void ) const
{
    
    static TypeSpec type_spec = getClassTypeSpec();
    
    return type_spec;
}


/** Print a simplified representation of the CorrespondenceAnalysis for the user. */
void CorrespondenceAnalysis::printValue(std::ostream &o) const
{
    
    o << "CorrespondenceAnalysis";
    
}


/** Set a member variable */
void CorrespondenceAnalysis::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    
    if ( name == "data" )
    {
        data = var;
    }
    else if ( name == "numAxes" )
    {
        num_axes = var;
    }
    else if ( name == "tolerance" )
    {
        tolerance = var;
    }
    else
    {
        RevObject::setConstParameter(name, var);
    }
    
}

