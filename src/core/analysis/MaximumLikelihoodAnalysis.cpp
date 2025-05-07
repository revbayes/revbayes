#include <cstddef>
#include <ostream>

#include "RandomNumberFactory.h"
#include "RbConstants.h"
#include "MaximumLikelihoodAnalysis.h"
#include "RlUserInterface.h"
#include "Cloneable.h"
#include "MaximumLikelihoodEstimation.h"
#include "Parallelizable.h"

namespace RevBayesCore { class Model; }

#ifdef RB_MPI
#include <mpi.h>
#endif


using namespace RevBayesCore;


/**
 * Constructor.
 *
 * \param[in]    m    The monte carlo sampler.
 */
MaximumLikelihoodAnalysis::MaximumLikelihoodAnalysis(MaximumLikelihoodEstimation *m) : Cloneable(), Parallelizable(),
    estimator( m )
{

    estimator->setActivePID( active_PID, num_processes );

}


MaximumLikelihoodAnalysis::MaximumLikelihoodAnalysis(const MaximumLikelihoodAnalysis &a) : Cloneable(), Parallelizable(a),
    estimator(NULL)
{

    if ( a.estimator != NULL )
    {
        estimator = a.estimator->clone();
    }

}




/**
 * Destructor. Nothing to do here
 */
MaximumLikelihoodAnalysis::~MaximumLikelihoodAnalysis(void)
{

    // free the estimator
    delete estimator;

}


/**
 * Overloaded assignment operator.
 * We need to keep track of the MaximumLikelihoodEstimations
 */
MaximumLikelihoodAnalysis& MaximumLikelihoodAnalysis::operator=(const MaximumLikelihoodAnalysis &a)
{
    Parallelizable::operator=(a);

    if ( this != &a )
    {

        // free the runs
        delete estimator;
        estimator = NULL;

        // create replicate Monte Carlo samplers
        if ( a.estimator != NULL )
        {
            estimator = a.estimator->clone();
        }
    }

    return *this;
}



MaximumLikelihoodAnalysis* MaximumLikelihoodAnalysis::clone( void ) const
{

    return new MaximumLikelihoodAnalysis( *this );
}


void MaximumLikelihoodAnalysis::disableScreenMonitors(bool all)
{

    if ( all == true || process_active == false )
    {
        return estimator->disableScreenMonitor();
    }

}


void MaximumLikelihoodAnalysis::finishMonitors( void )
{

    estimator->finishMonitors();

}


//size_t MaximumLikelihoodAnalysis::getCurrentGeneration( void ) const
//{
//
//    return estimator->getCurrentGeneration();
//}


Model& MaximumLikelihoodAnalysis::getModel( void )
{

    return estimator->getModel();
}


const Model& MaximumLikelihoodAnalysis::getModel( void ) const
{

    return estimator->getModel();
}


void MaximumLikelihoodAnalysis::monitor( size_t i ) const
{

    estimator->monitor( i );
}


void MaximumLikelihoodAnalysis::run( double epsilon, bool verbose )
{

#ifdef RB_MPI
    MPI_Comm_split(MPI_COMM_WORLD, active_PID, pid, &analysis_comm);
#endif

    size_t tuning_interval = 100;
//    double min_acceptance_ratio = 0.01;
    double min_improvement = epsilon;

    // get the current generation
    size_t gen = estimator->getCurrentGeneration();

    estimator->disableScreenMonitor();

    // Let user know what we are doing
    if ( verbose == true )
    {
        std::stringstream ss;
        if ( estimator->getCurrentGeneration() == 0 )
        {
            ss << "\n";
            ss << "Running MLE simulation\n";
        }
        else
        {
            ss << "Appending to previous MLE simulation of " << estimator->getCurrentGeneration() << " iterations\n";
        }
        ss << estimator->getStrategyDescription();
        RBOUT( ss.str() );

        RBOUT( "\n" );
        RBOUT( "Step\t -- \tLnProbality" );
    }

    // Monitor
   estimator->startMonitors( 1, false );
   estimator->writeMonitorHeaders();
   estimator->monitor(0);

    // reset the counters for the move schedules
    estimator->reset();

    // Run the chain
    bool converged = false;
    double previous_ln_likelihood = RbConstants::Double::neginf;
    do {
        ++gen;

        estimator->nextCycle();

        // check for autotuning and convergence
        converged = false;
        if ( gen % tuning_interval == 0 )
        {

            double current_ln_likelihood = estimator->getModelLnProbability(false);
            converged = (current_ln_likelihood - previous_ln_likelihood) < min_improvement;

//            converged &= estimator->hasConverged( min_acceptance_ratio );

            previous_ln_likelihood = current_ln_likelihood;

            estimator->tune();

            if ( verbose == true )
            {
                std::stringstream ss;
                ss << gen << "\t -- \t" << current_ln_likelihood;
                RBOUT( ss.str() );
            }

        }


    } while ( converged == false );


   // Monitor
   estimator->monitor(1);
   // estimator->finishMonitors();


#ifdef RB_MPI
    // wait until all replicates complete
    MPI_Barrier( analysis_comm );
    MPI_Comm_free(&analysis_comm);
#endif

}


/**
 * Set the active PID of this specific Monte Carlo analysis.
 */
void MaximumLikelihoodAnalysis::setActivePIDSpecialized(size_t a, size_t n)
{
    estimator->setActivePID( a, n );
}


/**
 * Set the model by delegating the model to the Monte Carlo samplers (replicates).
 */
void MaximumLikelihoodAnalysis::setModel(Model *m)
{

    // reset the counters for the move schedules
    estimator->setModel( m );

}


void MaximumLikelihoodAnalysis::startMonitors( void )
{

    estimator->startMonitors( 1, false );
    estimator->writeMonitorHeaders();

}
