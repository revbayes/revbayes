#include "RlRateMatrix.h"

#include <string>
#include <cstdint>

#include "ArgumentRule.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RlMemberFunction.h"
#include "RlSimplex.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "MemberProcedure.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RateGenerator.h"
#include "RbVector.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlTypedFunction.h"
#include "Simplex.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

RateMatrix::RateMatrix(void) : RateGenerator()
{
    initMethods();
}


RateMatrix::RateMatrix( const RevBayesCore::RateMatrix &v) : RateGenerator( v.clone() )
{
    initMethods();
}


RateMatrix::RateMatrix( RevBayesCore::RateMatrix *v) : RateGenerator( v )
{
    initMethods();
}


RateMatrix::RateMatrix( RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator> *m) : RateGenerator( m )
{
    initMethods();
}


RateMatrix* RateMatrix::clone() const
{
    return new RateMatrix( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RateMatrix::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{

    if (name == "size")
    {
        found = true;

        std::int64_t n = this->dag_node->getValue().getNumberOfStates();
        return new RevVariable( new Natural(n) );
    }

    return RateGenerator::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& RateMatrix::getClassType(void) {

    static std::string rev_type = "RateMatrix";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& RateMatrix::getClassTypeSpec(void) {

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RateGenerator::getClassTypeSpec() ) );

	return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RateMatrix::getTypeSpec(void) const {

    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}

void RateMatrix::initMethods(void) {

    // member procedures
    ArgumentRules* sizeArgRules = new ArgumentRules();
    methods.addFunction( new MemberProcedure( "size", Natural::getClassTypeSpec(), sizeArgRules) );

    // member functions
    ArgumentRules* squareBracketArgRules = new ArgumentRules();
    squareBracketArgRules->push_back( new ArgumentRule( "index" , Natural::getClassTypeSpec(), "The index of the row.", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    methods.addFunction( new MemberFunction<RateMatrix, ModelVector<Real> >( "[]", this, squareBracketArgRules   ) );
    
    // member functions
    ArgumentRules* stationaryFrequenciesArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<RateMatrix, Simplex>( "getStationaryFrequencies", this, stationaryFrequenciesArgRules  ) );
    
}
