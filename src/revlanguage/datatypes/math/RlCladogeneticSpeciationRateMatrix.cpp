//
//  RlCladogeneticSpeciationRateMatrix.cpp
//
//  Created by Will Freyman on 8/1/2017.
//


#include "ArgumentRule.h"
#include "CladogeneticSpeciationRateMatrix.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlCladogeneticSpeciationRateMatrix.h"
#include "RlRateMatrix.h"

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
