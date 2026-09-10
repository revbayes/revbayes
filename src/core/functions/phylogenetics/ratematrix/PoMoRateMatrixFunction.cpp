#include "PoMoRateMatrixFunction.h"

#include "RateMatrix_PoMo.h"
#include "Cloneable.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PoMoRateMatrixFunction::PoMoRateMatrixFunction(const TypedDagNode< std::int64_t > *ps, const TypedDagNode< RbVector<double> > *mr, const TypedDagNode< RbVector<double>  > *sc) : TypedFunction<RateGenerator>( new RateMatrix_PoMo(4 + 6*(ps->getValue() - 1), ps->getValue(), mr->getValue(), sc->getValue()) ),
    populationSize( ps ),
    mutationRates( mr ),
    selectionCoefficients ( sc )
{

    useMutationMatrix = false;
    // add the lambda parameter as a parent
    addParameter( populationSize );
    addParameter( mutationRates );
    addParameter( selectionCoefficients );

    update();
}


 //MJL 140822: caused compile error
PoMoRateMatrixFunction::PoMoRateMatrixFunction(const TypedDagNode< std::int64_t > *ps, const TypedDagNode< RateGenerator > *mm, const TypedDagNode< RbVector<double>  > *sc) : TypedFunction<RateGenerator>( new RateMatrix_PoMo(4 + 6*(ps->getValue() - 1), ps->getValue(), mm->getValue(), sc->getValue()) ), populationSize( ps ), mutationMatrix( mm ), selectionCoefficients ( sc ) {
    useMutationMatrix = true;
    // add the lambda parameter as a parent
    addParameter( populationSize );
    addParameter( mutationMatrix );
    addParameter( selectionCoefficients );
    
    update();
}


PoMoRateMatrixFunction::~PoMoRateMatrixFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PoMoRateMatrixFunction* PoMoRateMatrixFunction::clone( void ) const
{
    return new PoMoRateMatrixFunction( *this );
}


void PoMoRateMatrixFunction::update( void )
{
    //For the moment, we do not allow that the population size could be updated.
    std::vector<double> r ;
    if (useMutationMatrix)
    {
        // get the information from the arguments for reading the file
        r = setMutationRates( mutationMatrix->getValue() );
    }
    else
    {
        // get the information from the arguments for reading the file
        r = mutationRates->getValue();

    }

    const std::vector<double>& s = selectionCoefficients->getValue();
    
    
    // set the base frequencies
    static_cast< RateMatrix_PoMo* >(value)->setMutationRates( r );
    static_cast< RateMatrix_PoMo* >(value)->setSelectionCoefficients( s );
    static_cast< RateMatrix_PoMo* >(value)->update();
    
}



void PoMoRateMatrixFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == mutationRates)
    {
        mutationRates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == selectionCoefficients)
    {
        selectionCoefficients = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    
}


std::vector<double> PoMoRateMatrixFunction::setMutationRates(const RateGenerator& mm)
{
    
    double age = 0.0;
    double rate = 1.0;
    std::vector<double> r;
    r.push_back( mm.getRate(0,1,age,rate) );
    r.push_back( mm.getRate(0,2,age,rate) );
    r.push_back( mm.getRate(0,3,age,rate) );
    r.push_back( mm.getRate(1,0,age,rate) );
    r.push_back( mm.getRate(1,2,age,rate) );
    r.push_back( mm.getRate(1,3,age,rate) );
    r.push_back( mm.getRate(2,0,age,rate) );
    r.push_back( mm.getRate(2,1,age,rate) );
    r.push_back( mm.getRate(2,3,age,rate) );
    r.push_back( mm.getRate(3,0,age,rate) );
    r.push_back( mm.getRate(3,1,age,rate) );
    r.push_back( mm.getRate(3,2,age,rate) );
    return r;
}

