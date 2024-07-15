//
//  PoMoReversibleMutationRatesFunction.cpp
//  RevBayesCore
//
//  Created by Rui Borges 
//  Copyright 2024
//

#include <cstddef>
#include <vector>

#include "PoMoReversibleMutationRatesFunction.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbVector.h"
//#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;

PoMoReversibleMutationRatesFunction::PoMoReversibleMutationRatesFunction(   long na, 
                                                                            const TypedDagNode< Simplex > *bf, 
                                                                            const TypedDagNode< RbVector<double> > *ex) : 
    TypedFunction< RbVector<double> > ( new RbVector<double>( na*(na - 1), 0.0) ),
    K( na ),
    pi( bf ), 
    rho( ex )
{

    //long n_mutation_rates  = na*(na-1);

    //for (size_t i = 0; i < n_mutation_rates; ++i) {
    //    static_cast< std::vector<double>* >(value)->push_back(0.0);
    //}

    // add the lambda parameter as a parent
    addParameter( pi );
    addParameter( rho );

    update();
}



PoMoReversibleMutationRatesFunction::~PoMoReversibleMutationRatesFunction( void ) {
    // We don't delete the parameters, because they might be used somewhere else too. The model needs to do that!
}



PoMoReversibleMutationRatesFunction* PoMoReversibleMutationRatesFunction::clone( void ) const
{
    return new PoMoReversibleMutationRatesFunction( *this );
}


void PoMoReversibleMutationRatesFunction::update( void )
{

    // get the information from the arguments
    long n_alleles = K;
    const std::vector<double>& base_frequencies   = pi  -> getValue();
    const std::vector<double>& exchangeabilitites = rho -> getValue();

    // calculating the mutation rates
    size_t i_mutation_rates    = 0;
    size_t i_exchangeabilities = 0;

    for (int i =0; i<n_alleles; ++i){
        for (int j =(i+1); j<n_alleles; ++j){

            // calculating the mutation rates: mu_ij=rho_ij*pi_j
            // following the order expected by PoMoNK
            // mu : Vector of mutation rates: mu=(mu_a0a1,mu_a1a0,mu_a0a2,mu_a2a0,...)
            (*value)[i_mutation_rates  ] = exchangeabilitites[i_exchangeabilities]*base_frequencies[j];
            (*value)[i_mutation_rates+1] = exchangeabilitites[i_exchangeabilities]*base_frequencies[i];

            // updating the indexes
            i_mutation_rates    += 2;
            i_exchangeabilities += 1;
        }
    }
    
    // some more testing is needed
    // std::cout << "Seems to be working: test for differente K 2 and 3 and nonuniform pi and rho!\n";

}




void PoMoReversibleMutationRatesFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    

    if (oldP == pi)
    {
        pi = static_cast<const TypedDagNode< Simplex >* >( newP );
    }
    
    if (oldP == rho)
    {
        rho = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

}



