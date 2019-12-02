
//
//  RlCladogeneticProbabilityMatrix.cpp
//  revbayes-proj
//
//  Created by Michael Landis on 10/23/16.
//  Copyright © 2016 Michael Landis. All rights reserved.
//

#include "ArgumentRule.h"
#include "CladogeneticProbabilityMatrix.h"
#include "ModelVector.h"
#include "Natural.h"
#include "Real.h"
#include "RealPos.h"
#include "RlBoolean.h"
#include "RlCladogeneticProbabilityMatrix.h"
#include "RlRateMatrix.h"

using namespace RevLanguage;

CladogeneticProbabilityMatrix::CladogeneticProbabilityMatrix(void) : ModelObject<RevBayesCore::CladogeneticProbabilityMatrix>()
{
    initMethods();
}


CladogeneticProbabilityMatrix::CladogeneticProbabilityMatrix( const RevBayesCore::CladogeneticProbabilityMatrix &v) : ModelObject<RevBayesCore::CladogeneticProbabilityMatrix>( v.clone() )
{
    initMethods();
}


CladogeneticProbabilityMatrix::CladogeneticProbabilityMatrix( RevBayesCore::CladogeneticProbabilityMatrix *v) : ModelObject<RevBayesCore::CladogeneticProbabilityMatrix>( v )
{
    initMethods();
}


CladogeneticProbabilityMatrix::CladogeneticProbabilityMatrix( RevBayesCore::TypedDagNode<RevBayesCore::CladogeneticProbabilityMatrix> *m) : ModelObject<RevBayesCore::CladogeneticProbabilityMatrix>( m )
{
    initMethods();
}


CladogeneticProbabilityMatrix* CladogeneticProbabilityMatrix::clone() const
{
    return new CladogeneticProbabilityMatrix( *this );
}


/* Map calls to member methods */
RevPtr<RevVariable> CladogeneticProbabilityMatrix::executeMethod(std::string const &name, const std::vector<Argument> &args, bool &found)
{
    ; // do nothing for now
    return ModelObject<RevBayesCore::CladogeneticProbabilityMatrix>::executeMethod( name, args, found );
}


/* Get Rev type of object */
const std::string& CladogeneticProbabilityMatrix::getClassType(void) {
    
    static std::string rev_type = "CladogeneticProbabilityMatrix";
    
    return rev_type;
}

/* Get class type spec describing type of object */
const TypeSpec& CladogeneticProbabilityMatrix::getClassTypeSpec(void) {
    
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( RevObject::getClassTypeSpec() ) );
    
    return rev_type_spec;
}


/** Get the type spec of this class. We return a member variable because instances might have different element types. */
const TypeSpec& CladogeneticProbabilityMatrix::getTypeSpec(void) const {
    
    static TypeSpec type_spec = getClassTypeSpec();
    return type_spec;
}

void CladogeneticProbabilityMatrix::initMethods(void) {
    ; // do nothing
}
