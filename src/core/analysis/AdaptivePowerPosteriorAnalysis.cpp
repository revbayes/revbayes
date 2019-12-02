#include <stddef.h>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <typeinfo>
#include <map>

#include "MonteCarloSampler.h"
#include "MoveSchedule.h"
#include "MpiUtilities.h"
#include "AdaptivePowerPosteriorAnalysis.h"
#include "ProgressBar.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "Cloneable.h"
#include "MonteCarloAnalysisOptions.h"
#include "Parallelizable.h"
#include "StringUtilities.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;

AdaptivePowerPosteriorAnalysis::AdaptivePowerPosteriorAnalysis(MonteCarloSampler *m, const std::string &fn, size_t k) : Cloneable( ), Parallelizable(),
    filename( fn ),
    old_marginal_likelihood( RbConstants::Double::neginf ),
    sampler( m ),
    sampleFreq( 100 ),
    processors_per_likelihood( k )
{
 
    // initialize the powers and the maps
    powers.clear();
    powers.push_back(1.0);
    powers.push_back(0.5);
    powers.push_back(0.0);
    
    average_likelihoods.clear();
    average_likelihoods.insert( std::pair<double, double>(1.0, 0.0) );
    average_likelihoods.insert( std::pair<double, double>(0.5, 0.0) );
    average_likelihoods.insert( std::pair<double, double>(0.0, 0.0) );

    dirty_stones.clear();
    dirty_stones.insert( std::pair<double, bool>(1.0, true) );
    dirty_stones.insert( std::pair<double, bool>(0.5, true) );
    dirty_stones.insert( std::pair<double, bool>(0.0, true) );

    indices.clear();
    indices.insert( std::pair<double, size_t>(1.0, 1) );
    indices.insert( std::pair<double, size_t>(0.5, 2) );
    indices.insert( std::pair<double, size_t>(0.0, 3) );

    initMPI();
    
}


AdaptivePowerPosteriorAnalysis::AdaptivePowerPosteriorAnalysis(const AdaptivePowerPosteriorAnalysis &a) : Cloneable( a ), Parallelizable( a ),
    filename( a.filename ),
    old_marginal_likelihood( a.old_marginal_likelihood ),
    powers( a.powers ),
    average_likelihoods( a.average_likelihoods ),
    dirty_stones( a.dirty_stones ),
    indices( a.indices ),
    sampler( a.sampler->clone() ),
    sampleFreq( a.sampleFreq ),
    tolerance( a.tolerance ),
    processors_per_likelihood( a.processors_per_likelihood )
{
    
}


AdaptivePowerPosteriorAnalysis::~AdaptivePowerPosteriorAnalysis(void)
{
    delete sampler;
}


/**
 * Overloaded assignment operator.
 * We need to keep track of the MonteCarloSamplers
 */
AdaptivePowerPosteriorAnalysis& AdaptivePowerPosteriorAnalysis::operator=(const AdaptivePowerPosteriorAnalysis &a)
{
    Parallelizable::operator=( a );
    
    if ( this != &a )
    {
        
        // free the sampler
        delete sampler;
        
        filename                        = a.filename;
        old_marginal_likelihood         = a.old_marginal_likelihood;
        powers                          = a.powers;
        average_likelihoods             = a.average_likelihoods;
        dirty_stones                    = a.dirty_stones;
        indices                         = a.indices;
        sampler                         = a.sampler->clone();
        sampleFreq                      = a.sampleFreq;
        tolerance                       = a.tolerance;
        processors_per_likelihood       = a.processors_per_likelihood;
        
    }
    
    return *this;
}


/** Run burnin and autotune */
void AdaptivePowerPosteriorAnalysis::burnin(size_t generations, size_t tuningInterval)
{
    
//    initMPI();
    
    // Initialize objects needed by chain
    sampler->initializeSampler();
    
    // reset the counters for the move schedules
    sampler->reset();
    
	// start the progress bar
	ProgressBar progress = ProgressBar(generations, 0);

    if ( process_active == true )
    {
		// Let user know what we are doing
        std::stringstream ss;
        ss << "\n";
        ss << "Running burn-in phase of Adaptive Power Posterior sampler for " << generations << " iterations.\n";
        ss << sampler->getStrategyDescription();
        std::cout << ss.str() << std::endl;
    
        // Print progress bar (68 characters wide)
        progress.start();
    }
    
    // Run the chain
    size_t numStars = 0;
    for (size_t k=1; k<=generations; k++)
    {
		if ( process_active == true )
        {
            progress.update( k );
        }
        
        sampler->nextCycle(false);
            
        // check for autotuning
        if ( k % tuningInterval == 0 && k != generations )
        {
                
            sampler->tune();
            
        }
        
    }
    
    
    if ( process_active == true )
    {
		progress.finish();
    }
    
}



AdaptivePowerPosteriorAnalysis* AdaptivePowerPosteriorAnalysis::clone( void ) const
{
    
    return new AdaptivePowerPosteriorAnalysis( *this );
}


void AdaptivePowerPosteriorAnalysis::initMPI( void )
{
    
    
    size_t active_proc = size_t( floor( pid   / double(processors_per_likelihood)) ) * processors_per_likelihood;
    sampler->setActivePID( active_proc, processors_per_likelihood );
    
    
#ifdef RB_MPI
    // wait until all processes are complete
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
}


void AdaptivePowerPosteriorAnalysis::runAll(size_t gen, size_t bg, size_t ti)
{

//    initMPI();
    
    if( gen < sampleFreq )
    {
        throw(RbException("Trying to run adaptive power posterior analysis for fewer generations than sampleFreq, no samples will be stored"));
    }

    // disable the screen monitor(s) if any
    sampler->disableScreenMonitor(true, 0);
    
    // print some information to the screen but only if we are the active process
    if ( process_active == true )
    {
        std::cout << std::endl;
        std::cout << "Running adaptive power posterior analysis ..." << std::endl;
    }
    
    double new_marginal_likelihood = 0.0;
    bool first = true;
    while ( true ) {
        
        // run the analysis for each dirty stone
        for(size_t i = 0; i < powers.size(); ++i)
        {
            
            // get the power for this stone
            double b = powers.at(i);
            
            // if the stone is dirty, run the sampler on this stone
            if (dirty_stones.find(b)->second == true)
            {
                runStone(i, gen, bg, ti);
            }
            
        }
        
        if ( first == true )
        {
            
            // compute the marginal likelihood with the trapezoid rule
            // for the "first" estimate of the marginal likelihood
            
            double a = powers.at(0);
            double b = powers.at(2);
            double c = powers.at(1);
            
            double fa = average_likelihoods.find(a)->second;
            double fb = average_likelihoods.find(b)->second;
            double fc = average_likelihoods.find(c)->second;
            
            old_marginal_likelihood = trapezoidRule(fa, fc, a, c) + trapezoidRule(fc, fb, c, b);
            
            first = false;
            
        }
        
        // compute the new estimate of the marginal likelihood
        // and determine placement of new stones
        new_marginal_likelihood = adapt();
        
        double diff = fabs(new_marginal_likelihood - old_marginal_likelihood);
        
//        std::cout << old_marginal_likelihood << " -- " << new_marginal_likelihood << " -- " << fabs(new_marginal_likelihood - old_marginal_likelihood) << std::endl;
        
        std::cout << "\t" << diff << std::endl;
        
        // check for convergence of the marginal likelihood
        if ( diff < tolerance ) {
            std::cout << std::endl;
            break;
        }
        
        // add the new stones
//        std::cout << std::endl << "Adapting power posterior analysis ..." << std::endl;

        addStones();
        
        // update the marginal likelihood estimate
        old_marginal_likelihood = new_marginal_likelihood;
        
    }
    
    if ( process_active == true )
    {
        summarizeStones();
    }
    
}

void AdaptivePowerPosteriorAnalysis::runStone(size_t idx, size_t gen, size_t bg, size_t ti)
{
    
    // get the power for this sampler
    double beta  = powers[idx];
    size_t index = indices.find(beta)->second;
    size_t file_index = index - 1;

    // create the directory if necessary
    RbFileManager fm = RbFileManager(filename);
    if (fm.getFileName() == "")
    {
        throw(RbException("Please provide a filename with an extension"));
    }
    std::string stoneFileName = fm.getFileNameWithoutExtension() + "_stone_" + file_index + "." + fm.getFileExtension();

    RbFileManager f = RbFileManager(fm.getFilePath(), stoneFileName);
    f.createDirectoryForFile();
    
    std::fstream outStream;
    outStream.open( f.getFullFileName().c_str(), std::fstream::out);
    outStream << "state\t" << "power\t" << "likelihood" << std::endl;
    
    // reset the sampler
    sampler->reset();

//    size_t burnin = size_t( ceil( 0.25*gen ) );
    
    size_t printInterval = size_t( round( fmax(1,(gen+bg)/40.0) ) );
    size_t digits = size_t( ceil( log10( powers.size() ) ) );
    
    // print output for users
    if ( process_active ==true )
    {
        std::cout << "Step ";
        for (size_t d = size_t( ceil( log10( index+1.1 ) ) ); d < digits; d++ )
        {
            std::cout << " ";
        }
//        std::cout << index << ", power = " << beta;
        std::cout << index;
        std::cout << "\t\t";
    }
    
    // set the power of this sampler
    sampler->setLikelihoodHeat( beta );
    
    // get the average likelihood for this stone
    double& this_average_likelihood = average_likelihoods.find(beta)->second;
    
    std::stringstream ss;
    ss << "_stone_" << file_index;
    sampler->addFileMonitorExtension( ss.str(), false);
    
    // run a burnin
    sampler->tune();
    for (size_t k=1; k<=bg; k++)
    {
        if ( process_active == true )
        {
            if ( k % printInterval == 0 )
            {
                std::cout << "-";
                std::cout.flush();
            }
//            size_t progress = 68 * (double) k / (double) generations;
//            if ( progress > numStars )
//            {
//                for ( ;  numStars < progress; ++numStars )
//                {
//                    std::cout << "*";
//                }
//                
//                std::cout.flush();
//            }
        }
        
        sampler->nextCycle(false);
        
        // check for autotuning
        if ( k % ti == 0 && k != bg )
        {
            sampler->tune();
        }
        
    }

    
    // Monitor
    sampler->startMonitors(gen, false);
    sampler->writeMonitorHeaders( false );
    sampler->monitor(0);
    
    int num_sampled_likelihoods = 0;
    double total_likelihood = 0.0;
    for (size_t k=1; k<=gen; ++k)
    {
        
        if ( process_active == true )
        {
            if ( k % printInterval == 0 )
            {
                std::cout << "*";
                std::cout.flush();
            }
        }
        
        sampler->nextCycle( true );

        // Monitor
        sampler->monitor(k);
        
        // sample the likelihood
//        if ( k > burnin && k % sampleFreq == 0 )
        if ( k % sampleFreq == 0 )
        {
            // compute the joint likelihood
            double likelihood = sampler->getModelLnProbability(true);
            
            num_sampled_likelihoods++;
            total_likelihood += likelihood;
            
            outStream << k << "\t" << beta << "\t" << likelihood << std::endl;
        }
            
    }
    
    this_average_likelihood = total_likelihood / num_sampled_likelihoods;
    
    int num_dirty = 0;
    for(int i = 0; i < powers.size(); ++i)
    {
        if( dirty_stones.find(powers.at(i))->second == true ) {
            num_dirty++;
        }
    }
    bool last_stone = num_dirty == 1;
    
    if ( process_active == true & last_stone == false )
    {
        std::cout << std::endl;
    }
    
    outStream.close();
    
    // Monitor
    sampler->finishMonitors( 1, MonteCarloAnalysisOptions::NONE );
    
    // set this stone as dirty
    dirty_stones.find(beta)->second = false;
    
}

double AdaptivePowerPosteriorAnalysis::adapt( void )
{
    
    double new_lnl = 0.0;
    
    double worst_fit = RbConstants::Double::neginf;

    // compute the integral
    // and figure out where the fit is the worst
    double a, b, c, fa, fb, fc;
    double this_trap, this_simp, fit;
    for(size_t i = 1; i < powers.size(); i += 2 )
    {
        
        a = powers.at(i-1);
        b = powers.at(i+1);
        c = powers.at(i);
        
        fa = average_likelihoods.find(a)->second;
        fb = average_likelihoods.find(b)->second;
        fc = average_likelihoods.find(c)->second;
        
        this_trap = trapezoidRule(fa, fc, a, c) + trapezoidRule(fc, fb, c, b);
        this_simp = simpsonsRule(fa, fb, fc, a, b, c);
        
        fit = fabs(this_trap - this_simp);
        
        if (fit > worst_fit)
        {
            worst_fit = fit;
            worst_stone = i;
        }
        
        new_lnl += this_simp;
        
    }
    
    return new_lnl;
    
    
}

void AdaptivePowerPosteriorAnalysis::addStones(void)
{

    size_t num_stones = powers.size();
    double a = powers.at(worst_stone-1);
    double b = powers.at(worst_stone+1);
    double c = powers.at(worst_stone);

    double new_stone_left  = (a + c) / 2;
    double new_stone_right = (c + b) / 2;

    powers.push_back(new_stone_left);
    average_likelihoods.insert( std::pair<double, double>(new_stone_left, 0.0) );
    dirty_stones.insert( std::pair<double, bool>(new_stone_left, true) );
    indices.insert( std::pair<double, size_t>(new_stone_left, num_stones+1) );

    powers.push_back(new_stone_right);
    average_likelihoods.insert( std::pair<double, double>(new_stone_right, 0.0) );
    dirty_stones.insert( std::pair<double, bool>(new_stone_right, true) );
    indices.insert( std::pair<double, size_t>(new_stone_right, num_stones+2) );

    // sort the powers in decreasing order
    std::sort(powers.begin(), powers.end(), std::greater<double>());

}

double AdaptivePowerPosteriorAnalysis::trapezoidRule(double fa, double fb, double a, double b)
{
    
    double h = fabs(b - a) / 2.0;
    double s = h * (fb + fa);
    return s;
    
}

double AdaptivePowerPosteriorAnalysis::simpsonsRule(double fa, double fb, double fc, double a, double b, double c)
{

//    double h = fabs(b - a) / 6.0;
//    double s = h * (fa + 4.0 * fc + fb);
//    return s;

    // simpson's rule for non-equal intervals
    // otherwise, use the generalized simpson's rule
    double x0 = b;
    double x1 = c;
    double x2 = a;
    
    double y0 = fb;
    double y1 = fc;
    double y2 = fa;
    
//    double x0 = powers[i+1];
//    double x1 = powers[i];
//    double x2 = powers[i-1];
//    
//    double y0 = pathValues[i+1];
//    double y1 = pathValues[i];
//    double y2 = pathValues[i-1];
//    
    double s = (x2 - x0) * (y0 + ( (x2 - x0) / (x1 - x0) ) * (y1 - y0) / 2 );
    s += (2.0 * x2 * x2 - x0 * x2 - x0 * x0 + 3 * x0 * x1 - 3 * x1 * x2 ) * ( ((y2 - y1) / (x2 - x1)) - ((y1 - y0) / (x1 - x0)) ) / 6.0 ;

    return s;
    
}

void AdaptivePowerPosteriorAnalysis::summarizeStones( void )
{
    // create the directory if necessary
    RbFileManager f = RbFileManager(filename);
    f.createDirectoryForFile();
    
    std::ofstream outStream;
    outStream.open( f.getFullFileName().c_str(), std::fstream::out);
    outStream << "state\t" << "power\t" << "likelihood" << std::endl;

    // Append each stone
    for (size_t idx = 0; idx < powers.size(); ++idx)
    {
        
        // get the beta for this stone
        double beta = powers.at(idx);
        
        // get the index for this stone
        size_t stone_index = indices.find(beta)->second;
        size_t file_index  = stone_index - 1;
        
        RbFileManager fm = RbFileManager(filename);
        std::string stoneFileName = fm.getFileNameWithoutExtension() + "_stone_" + file_index + "." + fm.getFileExtension();
        
        RbFileManager f = RbFileManager(fm.getFilePath(), stoneFileName);

        // read the i-th stone
        std::ifstream inStream;
        inStream.open( f.getFullFileName().c_str(), std::fstream::in);
        if (inStream.is_open())
        {
            bool header = true;
            std::string line = "";
            while ( std::getline(inStream,line) )
            {
                // we need to skip the header line
                if ( header == true )
                {
                    header  = false;
                }
                else
                {
                    outStream << line << std::endl;
                }
                
            }
            inStream.close();
        }
        else
        {
            std::cerr << "Problem reading stone " << stone_index << " from file " << stoneFileName << "." << std::endl;
        }

    }
    
    // closing the file stream
    outStream.close();

}

void AdaptivePowerPosteriorAnalysis::setSampleFreq(size_t sf)
{
    sampleFreq = sf;
}

void AdaptivePowerPosteriorAnalysis::setTolerance(double tol)
{
    tolerance = tol;
}
