//
//  PoMoStationaryFrequenciesFunction.cpp
//  RevBayesCore
//
//  Created by Rui Borges 
//  Copyright 2024
//

#include <cstddef>

#include "PoMoStationaryFrequenciesFunction.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbVector.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PoMoStationaryFrequenciesFunction::PoMoStationaryFrequenciesFunction(   long na,
                                                                        long nv, 
                                                                        const TypedDagNode< double > *ne, 
                                                                        const TypedDagNode< Simplex > *bf, 
                                                                        const TypedDagNode< RbVector<double> > *ex) : TypedFunction< Simplex > ( new Simplex() ),
    K( na ),
    V( nv ),
    N( ne ), 
    pi( bf ), 
    rho( ex )
{

    long n_states  = na + na*(na-1)*(nv-1)/2;

    for (size_t i = 0; i < n_states; ++i) {
        static_cast< std::vector<double>* >(value)->push_back(0.0);
    }

    // add the lambda parameter as a parent
    addParameter( N );
    addParameter( pi );
    addParameter( rho );

    update();
}



PoMoStationaryFrequenciesFunction::~PoMoStationaryFrequenciesFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PoMoStationaryFrequenciesFunction* PoMoStationaryFrequenciesFunction::clone( void ) const
{
    return new PoMoStationaryFrequenciesFunction( *this );
}


void PoMoStationaryFrequenciesFunction::update( void )
{

    // get the information from the arguments
    long n_alleles = K;
    long n_virtual = V;
    const double& n_effective = (double)N -> getValue();
    const std::vector<double>& base_frequencies = pi  -> getValue();
    const std::vector<double>& exchangeabilitites = rho -> getValue();

    // harmonic numbers
    double harmonic_number_virtual = RbMath::fastHarmonicNumber(n_virtual-1.0);
    double harmonic_number_effective;
    
    // Talk to Sebastian: Set N=1 if nothing is given.
    double the_correction; 
    if (n_effective == 1) { 
        harmonic_number_effective = 1.0; 
        the_correction    = 1.0;
    } else {
        harmonic_number_effective = RbMath::fastHarmonicNumber(n_effective-1.0); 
        the_correction    = n_effective*harmonic_number_effective/(n_virtual*harmonic_number_virtual);
    }

    // normalization constant
    double normalization_constant = 0.0;
    size_t i_exchangeabilities = 0;
    for (int i =0; i<n_alleles; ++i){
        for (int j =(i+1); j<n_alleles; ++j){
            normalization_constant += base_frequencies[i]*base_frequencies[j]*exchangeabilitites[i_exchangeabilities];
            i_exchangeabilities += 1;
        }
    }
    normalization_constant = 1.0/(1.0 + 2.0*n_virtual*harmonic_number_virtual*normalization_constant*the_correction);

    // stationary frequencies
    // monomorphic states
    for (int i=0; i<n_alleles; ++i) {
        (*value)[i] = base_frequencies[i] * normalization_constant;
    }

    // polymorphic states
    i_exchangeabilities = 0;
    size_t i_states = n_alleles;
    double base;
    for (int i =0; i<n_alleles; ++i){
        for (int j =(i+1); j<n_alleles; ++j){
            base = base_frequencies[i]*base_frequencies[j]*exchangeabilitites[i_exchangeabilities]*n_virtual*n_virtual*normalization_constant*the_correction;
            i_exchangeabilities += 1;
            for (int v=1; v<n_virtual; ++v){
                (*value)[i_states] = base/(v*(V-v)); 
                i_states += 1;
            }

        }
    }

    std::cout << "Seems to be working: test a bit more with mut bias and compare with R code!\n";

}




void PoMoStationaryFrequenciesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == N)
    {
        N = static_cast<const TypedDagNode< double >* >( newP );
    }

    if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}



