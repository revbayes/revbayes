//
//  RlCladogeneticSpeciationRateMatrix.cpp
//
//  Created by Will Freyman on 8/1/2017.
//


#include <iosfwd>
#include <string>
#include <vector>

#include "CladogeneticSpeciationRateMatrix.h"
#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "MemberFunction.h"
#include "ModelObject.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
#include "RlMemberFunction.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"
#include "UserFunctionNode.h"

namespace RevLanguage { class Argument; }

using namespace RevLanguage;

CladogeneticSpeciationRateMatrix::CladogeneticSpeciationRateMatrix(void) : ModelObject<RevBayesCore::CladogeneticSpeciationRateMatrix>()
{
    initMethods();
}


CladogeneticSpeciationRateMatrix::CladogeneticSpeciationRateMatrix( const RevBayesCore::CladogeneticSpeciationRateMatrix &v) : ModelObject<RevBayesCore::CladogeneticSpeciationRateMatrix>( v.clone() )
{
    initMethods();
}


CladogeneticSpeciationRateMatrix::CladogeneticSpeciationRateMatrix( RevBayesCore::CladogeneticSpeciationRateMatrix *v) : ModelObject<RevBayesCore::CladogeneticSpeciationRateMatrix>( v )
{
    initMethods();
}


CladogeneticSpeciationRateMatrix::CladogeneticSpeciationRateMatrix( RevBayesCore::TypedDagNode<RevBayesCore::CladogeneticSpeciationRateMatrix> *m) : ModelObject<RevBayesCore::CladogeneticSpeciationRateMatrix>( m )
{
    initMethods();
}


CladogeneticSpeciationRateMatrix* CladogeneticSpeciationRateMatrix::clone() const
{
    return new CladogeneticSpeciationRateMatrix( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> CladogeneticSpeciationRateMatrix::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    
    if (name == "getSpeciationRateSumsPerState")
    {
        found = true;
        std::vector<double> v = this->dag_node->getValue().getSpeciationRateSumPerState();
        return new RevVariable( new ModelVector<RealPos>( v ) );
    }
    else if (name == "getCladogeneticProbabilityMatrix")
    {
        found = true;
        RevBayesCore::CladogeneticProbabilityMatrix p = this->dag_node->getValue().getCladogeneticProbabilityMatrix();
        return new RevVariable( new CladogeneticProbabilityMatrix( p ) );
    }
    else if (name == "getRate")
    {
        found = true;

        // create state for anc -> ch1, ch2
        std::vector<unsigned> state;
        unsigned anc_state = unsigned( static_cast<const Natural&>( args[0].getVariable()->getRevObject() ).getValue() );
        unsigned ch1_state = unsigned( static_cast<const Natural&>( args[1].getVariable()->getRevObject() ).getValue() );
        unsigned ch2_state = unsigned( static_cast<const Natural&>( args[2].getVariable()->getRevObject() ).getValue() );
        state.push_back(anc_state);
        state.push_back(ch1_state);
        state.push_back(ch2_state);

        // get rate event map
        const std::map< std::vector<unsigned>, double >& rates = this->dag_node->getValue().getEventMap();
        
        double v = 0.0;
        auto it = rates.find(state);
        if (it != rates.end()) {
            v = it->second;
        }

        return new RevVariable( new RealPos( v ) );
    }
    
    return ModelObject<RevBayesCore::CladogeneticSpeciationRateMatrix>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& CladogeneticSpeciationRateMatrix::getClassType(void) {
    
    static std::string rev_type = "CladogeneticSpeciationRateMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& CladogeneticSpeciationRateMatrix::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& CladogeneticSpeciationRateMatrix::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}

void CladogeneticSpeciationRateMatrix::initMethods(void) {
    
    ArgumentRules* speciationRateSumPerStateArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<CladogeneticSpeciationRateMatrix, ModelVector<RealPos> >( "getSpeciationRateSumPerState", this, speciationRateSumPerStateArgRules ) );

    ArgumentRules* cladogeneticProbabilityMatrixArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<CladogeneticSpeciationRateMatrix, CladogeneticProbabilityMatrix>( "getCladogeneticProbabilityMatrix", this, cladogeneticProbabilityMatrixArgRules ) );

    ArgumentRules* getEventsArgRules = new ArgumentRules();
    methods.addFunction( new MemberFunction<CladogeneticSpeciationRateMatrix, ModelVector<ModelVector<Natural> > >( "getEvents", this, getEventsArgRules ) );
    
    ArgumentRules* getRateArgRules = new ArgumentRules();
    getRateArgRules->push_back(new ArgumentRule( "ancestorState", Natural::getClassTypeSpec(), "The state of the ancestral lineage before cladogenesis.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    getRateArgRules->push_back(new ArgumentRule( "childState1", Natural::getClassTypeSpec(), "The state of the first child lineage after cladogenesis.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    getRateArgRules->push_back(new ArgumentRule( "childState2", Natural::getClassTypeSpec(), "The state of the second child lineage after cladogenesis.", ArgumentRule::BY_VALUE, ArgumentRule::ANY  ) );
    methods.addFunction( new MemberFunction<CladogeneticSpeciationRateMatrix, RealPos>( "getRate", this, getRateArgRules ) );
}
