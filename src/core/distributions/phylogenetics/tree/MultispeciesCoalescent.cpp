#include "Clade.h"
#include "ConstantNode.h"
#include "MultispeciesCoalescent.h"
#include "DistributionExponential.h"
#include "DistributionUniform.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "TopologyNode.h"

#include <algorithm>
#include <cmath>

using namespace RevBayesCore;

MultispeciesCoalescent::MultispeciesCoalescent(const TypedDagNode<Tree> *sp, const std::vector<Taxon> &t) : AbstractMultispeciesCoalescent(sp, t)
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
    
    double theta = 1.0 / getNe( index );
    
    double ln_prob_coal = 0;
    double current_time = begin_age;
    for (size_t i=0; i<times.size(); ++i)
    {
        
        // now we do the computation
        //a is the time between the previous and the current coalescences
        double a = times[i] - current_time;
        
        current_time = times[i];
        
        // get the number j of individuals we had before the current coalescence
        size_t j = k - i;
        
        // compute the number of pairs: pairs = C(j choose 2) = j * (j-1.0) / 2.0
        double n_pairs = j * (j-1.0) / 2.0;
        double lambda = n_pairs * theta;
        
        // add the density for this coalescent event
        //lnProbCoal += log( lambda ) - lambda * a;
        //Corrected version:
        ln_prob_coal += log( theta ) - lambda * a;
        
    }
    
    // compute the probability of no coalescent event in the final part of the branch
    // only do this if the branch is not the root branch
    if ( add_final_interval == true )
    {
        double final_interval = end_age - current_time;
        size_t j = k - times.size();
        ln_prob_coal -= j*(j-1.0)/2.0 * final_interval * theta;
        
    }
    
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
        std::cerr << "Error: Null Pointers for Ne and Nes." << std::endl;
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
