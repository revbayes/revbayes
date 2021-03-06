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
#include "ConstantNode.h"
#include "DagNode.h"
#include "DeterministicNode.h"
#include "DynamicNode.h"
#include "IndirectReferenceFunction.h"
#include "ModelObject.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlConstantNode.h"
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
    ; // do nothing for now
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
    ; // do nothing
}
