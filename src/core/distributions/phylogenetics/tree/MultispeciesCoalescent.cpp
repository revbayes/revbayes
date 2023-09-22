#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "MultispeciesCoalescent.h"
#include "AbstractMultispeciesCoalescent.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { class Taxon; }
namespace RevBayesCore { class Tree; }

using namespace RevBayesCore;

MultispeciesCoalescent::MultispeciesCoalescent(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : AbstractMultispeciesCoalescent(sp, t),
    Ne( NULL ),
    Nes( NULL )
{

}


MultispeciesCoalescent::~MultispeciesCoalescent()
{

}





MultispeciesCoalescent* MultispeciesCoalescent::clone( void ) const
{

    return new MultispeciesCoalescent( *this );
}


double MultispeciesCoalescent::computeLnCoalescentProbability(size_t k, const std::vector<double> &times, double begin_age, double end_age, size_t index, bool add_final_interval)
{
    if ( k == 1 ) return 0.0;

    double theta = 2.0 / getNe( index );

    // if (index == 0) {
    //     std::cout << "----------------\nstart age: " << begin_age << std::endl;
    //     std::cout << "end age: " << end_age << std::endl;
    //     std::cout << "N: " << getNe(index) << std::endl;
    // }

    double ln_prob_coal = 0;
    double current_time = begin_age;

    for (size_t i=0; i<times.size(); ++i)
    {
        // now we do the computation
        // a is the time between the previous and the current coalescences
        double a = times[i] - current_time;

        // std::cout << "t: " << a << std::endl;

        current_time = times[i];

        // get the number j of individuals we had before the current coalescence
        size_t j = k - i;

        // std::cout << "j: " << j << std::endl;

        // compute the number of pairs: pairs = C(j choose 2) = j * (j-1.0) / 2.0
        double n_pairs = j * (j-1.0) / 2.0;
        double lambda = n_pairs * theta;

        // add the density for this coalescent event
        //lnProbCoal += log( lambda ) - lambda * a;
        //Corrected version:

        // std::cout << "term: " << log( theta ) - lambda * a << std::endl;

        ln_prob_coal += log( theta ) - lambda * a;
    }

    // compute the probability of no coalescent event in the final part of the branch
    // only do this if the branch is not the root branch
    if ( add_final_interval == true )
    {
        double final_interval = end_age - current_time;
        size_t j = k - times.size();

        // std::cout << "final: " << final_interval << std::endl;
        // std::cout << "j: " << j << std::endl;

        double n_pairs = j * (j-1.0) / 2.0;
        ln_prob_coal -= n_pairs * theta * final_interval;

        // std::cout << "term: " << n_pairs * theta * final_interval << std::endl;
    }

    // std::cout << "ln prob coal: " << ln_prob_coal << std::endl;

    return ln_prob_coal;
}


double MultispeciesCoalescent::drawNe( size_t index )
{

    return getNe( index );
}



double  MultispeciesCoalescent::getNe(size_t index) const
{

    if ( Ne != NULL )
    {
        return Ne->getValue();
    }
    else if ( Nes != NULL )
    {
        return Nes->getValue()[index];
    }
    else
    {
        exit(-1);
    }
}


void MultispeciesCoalescent::setNes(TypedDagNode< RbVector<double> >* input_nes)
{

    removeParameter( Nes );
    removeParameter( Ne );

    Nes = input_nes;
    Ne  = NULL;

    addParameter( Nes );

}


void MultispeciesCoalescent::setNe(TypedDagNode<double>* input_ne)
{

    removeParameter( Ne );
    removeParameter( Nes );

    Ne  = input_ne;
    Nes = NULL;

    addParameter( Ne );
}


/** Swap a parameter of the distribution */
void MultispeciesCoalescent::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{

    if ( oldP == Nes )
    {
        Nes = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }

    if ( oldP == Ne )
    {
        Ne = static_cast<const TypedDagNode< double >* >( newP );
    }

    AbstractMultispeciesCoalescent::swapParameterInternal(oldP, newP);

}
