#include "FlowT2PopulationsRateMatrixFunction.h"
#include "RbException.h"

using namespace RevBayesCore;


FlowT2PopulationsRateMatrixFunction::FlowT2PopulationsRateMatrixFunction(   const TypedDagNode< RbVector<double> > *n,
                                                                            const TypedDagNode< RbVector<double> > *m,
                                                                            const TypedDagNode< RbVector<double> > *l,
                                                                            const TypedDagNode< RbVector<double> > *s ) : TypedFunction<RateGenerator>( new RateMatrix_FlowT2Populations( computeNumStates(n->getValue()), n->getValue() ) ),
nu( n ),
mu( m ),
lambda( l ),
sigma( s )

{
    // add the lambda parameter as a parent
    addParameter( nu );
    addParameter( mu );
    addParameter( lambda );
    addParameter( sigma );

    update();
}


FlowT2PopulationsRateMatrixFunction::~FlowT2PopulationsRateMatrixFunction( void )
{
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}


size_t FlowT2PopulationsRateMatrixFunction::computeNumStates( const std::vector<double> &n )
{

    size_t numStates = (n[0]+1)*(n[1]+1);
    std::cout << numStates;
    
    return numStates;

}



FlowT2PopulationsRateMatrixFunction* FlowT2PopulationsRateMatrixFunction::clone( void ) const
{
    
    return new FlowT2PopulationsRateMatrixFunction( *this );
}


void FlowT2PopulationsRateMatrixFunction::update( void )
{
    // get the information from the arguments for reading the file
    const std::vector<double>& n = nu->getValue();
    const std::vector<double>& m = mu->getValue();
    const std::vector<double>& l = lambda->getValue();
    const std::vector<double>& s = sigma->getValue();

    // set the base frequencies
    static_cast< RateMatrix_FlowT2Populations* >(value)->setNu( n );
    static_cast< RateMatrix_FlowT2Populations* >(value)->setMu( m );
    static_cast< RateMatrix_FlowT2Populations* >(value)->setLambda( l );
    static_cast< RateMatrix_FlowT2Populations* >(value)->setSigma( s );

    value->update();
}



void FlowT2PopulationsRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == nu)
    {
        nu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
    if (oldP == mu)
    {
        mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == lambda)
    {
        lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if (oldP == sigma)
    {
        sigma = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}




