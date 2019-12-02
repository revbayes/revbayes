#include "MixtureCladogeneticStateFunction.h"
#include "RbException.h"

#include <math.h>

using namespace RevBayesCore;

MixtureCladogeneticStateFunction::MixtureCladogeneticStateFunction(const TypedDagNode< Simplex > *mw, const TypedDagNode< RbVector<CladogeneticProbabilityMatrix> > *cp, unsigned nc, unsigned ns) : TypedFunction<CladogeneticProbabilityMatrix>( new CladogeneticProbabilityMatrix(0) ),
mixtureWeights( mw ),
cladoProbs( cp ),
numCharacters(nc),
numStates(ns),
numComponents( mw->getValue().size() )
{
    // add the lambda parameter as a parent
    addParameter( mixtureWeights );
    addParameter( cladoProbs );
    
    buildEventMap();
    
    update();
}

MixtureCladogeneticStateFunction::~MixtureCladogeneticStateFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


void MixtureCladogeneticStateFunction::buildEventMap( void ) {
    
    std::map<std::vector<unsigned>, double> eventMapProbsMixture = this->getValue().getEventMap();
    
//    this->getValue()[it->first] = 0.0;
    for (size_t i = 0; i < numComponents; i++)
    {
        const CladogeneticProbabilityMatrix& cp = cladoProbs->getValue()[i];
        std::map<std::vector<unsigned>, double> emp = cp.getEventMap();
        std::map<std::vector<unsigned>, double>::iterator it;
        for (it = emp.begin(); it != emp.end(); it++)
        {
            eventMapProbsMixture[it->first] = 0.0;
        }
    }
    
    this->getValue().setEventMap(eventMapProbsMixture);
}



MixtureCladogeneticStateFunction* MixtureCladogeneticStateFunction::clone( void ) const
{
    return new MixtureCladogeneticStateFunction( *this );
}

std::map< std::vector<unsigned>, double >  MixtureCladogeneticStateFunction::getEventMap(double t)
{
    return this->getValue().getEventMap();
}

const std::map< std::vector<unsigned>, double >&  MixtureCladogeneticStateFunction::getEventMap(double t) const
{
    return this->getValue().getEventMap();
}

void MixtureCladogeneticStateFunction::update( void )
{
    
    // get the information from the arguments for reading the file
    const std::vector<double>& mp = mixtureWeights->getValue();

    std::map<std::vector<unsigned>, double> eventMapProbsMixture = this->getValue().getEventMap();
    
    //    this->getValue()[it->first] = 0.0;
    for (size_t i = 0; i < numComponents; i++)
    {
        const CladogeneticProbabilityMatrix& cp = cladoProbs->getValue()[i];
        std::map<std::vector<unsigned>, double> emp = cp.getEventMap();
        std::map<std::vector<unsigned>, double>::iterator it;
        for (it = emp.begin(); it != emp.end(); it++)
        {
            eventMapProbsMixture[it->first] = mp[i] * it->second;
        }
    }
    
    this->getValue().setEventMap(eventMapProbsMixture);
    
    return;
}



void MixtureCladogeneticStateFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == mixtureWeights)
    {
        mixtureWeights = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    else if (oldP == cladoProbs)
    {
        cladoProbs = static_cast<const TypedDagNode< RbVector<CladogeneticProbabilityMatrix> >* >( newP );
    }
    
}
