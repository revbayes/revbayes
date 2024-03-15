#include <iosfwd>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlMemberFunction.h"
#include "RlRateGenerator.h"
#include "ArgumentRules.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "MethodTable.h"
#include "ModelObject.h"
#include "RateGenerator.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlDeterministicNode.h"
#include "RlTypedFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

RateGenerator::RateGenerator(void) : ModelObject<RevBayesCore::RateGenerator>()
{
    initMethods();
}


RateGenerator::RateGenerator( const RevBayesCore::RateGenerator &v) : ModelObject<RevBayesCore::RateGenerator>( v.clone() )
{
    initMethods();
}


RateGenerator::RateGenerator( RevBayesCore::RateGenerator *v) : ModelObject<RevBayesCore::RateGenerator>( v )
{
    initMethods();
}


RateGenerator::RateGenerator( RevBayesCore::TypedDagNode<RevBayesCore::RateGenerator> *m) : ModelObject<RevBayesCore::RateGenerator>( m )
{
    initMethods();
}


RateGenerator* RateGenerator::clone() const
{
    return new RateGenerator( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> RateGenerator::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    ; // do nothing for now
    return ModelObject<RevBayesCore::RateGenerator>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& RateGenerator::getClassType(void) {

    static std::string rev_type = "RateGenerator";

	return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& RateGenerator::getClassTypeSpec(void) {

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );

	return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& RateGenerator::getTypeSpec(void) const {

    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}

void RateGenerator::initMethods(void)
{
    
    
    // add the DAG node member methods
    // note that this is a sage case because all DAG nodes are member objects
    if ( dag_node != NULL )
    {
        const MethodTable &dagMethods = dynamic_cast<RevMemberObject*>( dag_node )->getMethods();
        methods.insertInheritedMethods( dagMethods );
    }
    
    // member functions
    ArgumentRules* transitionProbabilityArgRules = new ArgumentRules();
    transitionProbabilityArgRules->push_back( new ArgumentRule( "rate", RealPos::getClassTypeSpec(), "The rate of the process (or duration of the process assuming rate=1).", ArgumentRule::BY_VALUE, ArgumentRule::ANY ) );
    transitionProbabilityArgRules->push_back( new ArgumentRule( "startAge", RealPos::getClassTypeSpec(), "The start age of the process.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
    transitionProbabilityArgRules->push_back( new ArgumentRule( "endAge", RealPos::getClassTypeSpec(), "The end age of the process.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.0) ) );
    methods.addFunction( new MemberFunction<RateGenerator, ModelVector<ModelVector<RealPos> > >( "getTransitionProbabilities", this, transitionProbabilityArgRules   ) );
    
}

/**
 * Print the value of the matrix, respecting the formatting of the
 * model object elements. We do this by retrieving the elements, one
 * by one, and printing them using their own printValue implementation.
 * Among other things, this takes care of proper formatting.
 *
 * We make a perfectly safe const cast here, since we only utilize the
 * const printValue function of the element.
 */
void RateGenerator::printValue( std::ostream& o, bool user ) const
{

    this->getDagNode()->printValue( o, ",", -1, true, user, false );

}
