#include "MarginalLikelihoodEstimator.h"

#include <algorithm>
#include <cstdlib>
#include <ostream>
#include <string>

#include "DistributionGeometric.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "StringUtilities.h"

using namespace RevBayesCore;


MarginalLikelihoodEstimator::MarginalLikelihoodEstimator(const path &fn, const std::string &pn, const std::string &ln, const std::string &del) :
    filename( fn ),
    powers(),
    likelihoodSamples()
{
    
    setActivePID( 0, 1 );
    
    
    if ( process_active == true )
    {
        /********************/
        /* read in the file */
        /********************/
    
        if ( not is_regular_file(filename) )
        {
            throw RbException() << "Could not not file " << filename;
        }
    
    
        bool hasHeaderBeenRead = false;
        
        // Open file
        std::ifstream inFile( filename.string() );
        
        if ( !inFile )
        {
            throw RbException() << "Could not open file " << filename;
        }
    
        // Initialize
        std::string commandLine;
    
        //    RBOUT("Processing file \"" + filename + "\"");
    
        int powerColumnIndex = -1;
        int likelihoodColumnIndex = -1;
        size_t index = 0;
    
        double previousPower = -1.0;
        // loop over file content
        while ( inFile.good() )
        {

            // Read a line
            std::string line;
            safeGetline( inFile, line );
            
            // skip empty lines
            if (line.length() == 0)
            {
                continue;
            }
        
            // removing comments
            if (line[0] == '#')
            {
                continue;
            }
            
            // splitting every line into its columns
            std::vector<std::string> columns;
            StringUtilities::stringSplit(line, del, columns);
            
            // we assume a header at the first line of the file
            if ( hasHeaderBeenRead == false )
            {
            
                for (size_t j=0; j<columns.size(); j++)
                {
                
                    if ( columns[j] == pn )
                    {
                        powerColumnIndex = (int)j;
                    }
                    else if ( columns[j] == ln )
                    {
                        likelihoodColumnIndex = (int)j;
                    }
                
                }
            
                hasHeaderBeenRead = (powerColumnIndex != -1 && likelihoodColumnIndex != -1);
            
                continue;
            }

            // check for broken lines
            if( columns.size() <= powerColumnIndex || columns.size() <= likelihoodColumnIndex ) 
                throw RbException() << "Please check format of file " << filename << ", missing power and likelihood columns in some lines";

            double p, l;
            // get the power entry
            std::string tmp = columns[powerColumnIndex];
            try {
                p = std::stod(tmp);
            } catch (std::invalid_argument&) {
                throw RbException() << "Please check format of file " << filename << ", non-numeric input in power column";
            }
            if ( p != previousPower )
            {
                previousPower = p;
                powers.push_back( p );
                likelihoodSamples.push_back( std::vector<double>() );
                index++;
            }
        
            // get the likelihood entry
            tmp = columns[likelihoodColumnIndex];
            try {
                l = std::stod(tmp);
            } catch (std::invalid_argument&) {
                throw RbException() << "Please check format of file " << filename << ", non-numeric input in likelihood column";
}
            likelihoodSamples[index-1].push_back( l );

        }
    
        inFile.close();
    }
    
}


MarginalLikelihoodEstimator::~MarginalLikelihoodEstimator()
{
    
}


/**
 * This function is adapted from TraceNumeric::update(): the only difference is that it works directly on a vector of numeric
 * values representing recorded samples, as opposed to a Trace object. The calculation is equivalent to the one used in
 * Tracer, and also (except for the choice of the maximum lag value) to convenience:::essTracer.
 *
 * @return The effective sample size
 */
double MarginalLikelihoodEstimator::getESS(const std::vector<double> values) const
{
    size_t MAX_LAG = 1000;
    size_t samples = values.size();
    size_t maxLag = (samples - 1 < MAX_LAG ? samples - 1 : MAX_LAG);
    
    // get the mean
    double mean = 0.0;
    for (size_t i = 0; i < samples; ++i)
    {
        mean += values[i];
    }
    mean /= samples;
    
    double* gammaStat = new double[maxLag];
    for (size_t j = 0; j < maxLag; j++)
    {
        gammaStat[j] = 0;
    }
    double varStat = 0.0;
    
    for (size_t lag = 0; lag < maxLag; lag++) {
        for (size_t j = 0; j < samples - lag; j++) {
            double del1 = values[j] - mean;
            double del2 = values[j + lag] - mean;
            gammaStat[lag] += (del1 * del2);
        }

        gammaStat[lag] /= ((double) (samples - lag));

        if (lag == 0) {
            varStat = gammaStat[0];
        } else if (lag % 2 == 0) {
            if (gammaStat[lag - 1] + gammaStat[lag] > 0) {
                varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
            } else {
                maxLag = lag;
            }
        }
    }
    
    // autocorrelation time
    double act = varStat / gammaStat[0];

    // effective sample size
    double ess = samples / act;
    
    return ess;
}


std::int64_t MarginalLikelihoodEstimator::offsetModulo(std::int64_t i, std::int64_t n)
{
    return 1 + (i - 1) % n;
}


/**
 * Get the indices of those elements of the vector of sampled log likelihoods that will be used for a given block bootstrap replicate.
 * The function is based on mcmc3r:::make.ends() (see https://github.com/dosreislab/mcmc3r/blob/main/R/marginal-lhd.R#L493),
 * whose code was in turn lifted without modifications from the eponymous function in the package "boot", written by Angelo J. Canty
 * and corrected by B. D. Ripley (see https://github.com/cran/boot/blob/master/R/bootfuns.q#L3203). It takes
 * as its argument (1) a pair of integers representing the start index and length of a given block, and (2) the total length of the time
 * series from which the blocks are to be drawn. The series is treated as circular, i.e., if we reach the end before the number of selected
 * indices has reached the length of the block, we keep adding indices from the beginning,
 *
 * @param a Pair of integers denoting the start index of a block and its length, respectively
 * @param n Total length of the time series
 * @return Vector of indices to select from a time series (of sampled log likelihoods)
 */
std::vector<std::int64_t> MarginalLikelihoodEstimator::getIndices(std::pair<std::int64_t, std::int64_t> a, std::int64_t n)
{
    std::vector<std::int64_t> out;
    if (a.second != 0)
    {
        std::vector<std::int64_t> seq( a.second );
        for (size_t i = 0; i < a.second; i++)
        {
            seq[i] = a.first + i;
            out.push_back( offsetModulo(seq[i], n) );
        }
    }
    
    return out;
}


std::vector< std::vector< std::vector<double> > > MarginalLikelihoodEstimator::blockBootstrap(size_t repnum, double prop)
{
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // vector to be returned (repnum replicates for each of the k powers)
    std::vector< std::vector< std::vector<double> > > res( powers.size() );
    
    for (size_t i = 0; i < powers.size(); i++)
    {
        // vector of bootstrap replicates
        std::vector< std::vector<double> > res_inner( repnum );
        
        // we keep n and n_sim distinct for conceptual reasons: technically, the simulated time series could contain a different
        // number of samples than the original time series
        size_t n = likelihoodSamples[i].size();
        std::int64_t n_sim = n;
        
        for (size_t j = 0; j < repnum; j++)
        {
            // vector of resampled likelihood values to be returned for this bootstrap replicate
            std::vector<double> res_innermost( n );
            
            // we always use end correction, so the endpoint is once again equal to n rather than n - (n * prop) + 1
            size_t endpt = n;
            bool cont = true;
            
            std::vector<std::int64_t> len_tot(repnum, 0);
            std::vector< std::vector<std::int64_t> > lens;
            
            while (cont)
            {
                std::vector<std::int64_t> temp0( repnum );
                std::vector<std::int64_t> temp1( repnum );
                std::vector<bool> test( repnum );
                
                for (size_t k = 0; k < repnum; k++)
                {
                    temp0[k] = 1 + RbStatistics::Geometric::rv( 1/(n * prop), *rng );
                    temp1[k] = std::min( temp0[k], n_sim - len_tot[k] );
                }
                
                lens.push_back( temp1 );
                
                for (size_t k = 0; k < repnum; k++)
                {
                    len_tot[k] = len_tot[k] + temp1[k];
                    test[k] = len_tot[k] < n_sim;
                }
                
                cont = std::any_of(test.begin(), test.end(), [](bool x) { return x; });
            }
            
            size_t nn = lens.size();
            std::vector< std::vector<std::int64_t> > st;
            
            for (size_t k = 0; k < nn; k++)
            {
                std::vector<std::int64_t> elem( repnum );
                for(size_t l = 0; l < repnum; l++)
                {
                    // transform uniform01() to integer in range [1, endpt]
                    elem[l] = static_cast<std::int64_t>( std::floor(rng->uniform01() * endpt) ) + 1;
                }
                
                st.push_back(elem);
            }
            
            std::vector< std::pair<std::int64_t, std::int64_t> > ends;
            std::vector< std::vector<std::int64_t> > inds;
            
            for (size_t k = 0; k < nn; k++)
            {
                // paranoid pair constructor syntax compatible with C++14 and earlier
                std::pair<std::int64_t, std::int64_t> tmp0( st[k][0], lens[k][0] );
                ends.push_back( tmp0 );
            }
            
            for (size_t k = 0; k < nn; k++)
            {
                std::vector<std::int64_t> tmp1 = getIndices( ends[k], static_cast<std::int64_t>(n) );
                inds.push_back( tmp1 );
            }
            
            // flatten the nested vector of indices
            size_t total_size = 0;
            for (const auto& inner : inds)
            {
                total_size += inner.size();
            }
            
            std::vector<std::int64_t> inds_flattened;
            inds_flattened.reserve(total_size);
            
            for (const auto& inner : inds)
            {
                inds_flattened.insert(inds_flattened.end(), inner.begin(), inner.end());
            }
            
            // fill in the output vector; if the total length of the flattened index vector exceeds n_sim, truncate
            for (size_t k = 0; k < n_sim; k++)
            {
                size_t idx = inds_flattened[k];
                res_innermost[k] = likelihoodSamples[i][idx];
            }
            
            res_inner[j] = res_innermost;
        }
        
        res[i] = res_inner;
    }
    
    return res;
}
