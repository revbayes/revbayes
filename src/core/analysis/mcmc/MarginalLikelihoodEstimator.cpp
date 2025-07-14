#include "MarginalLikelihoodEstimator.h"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <string>

#include "DistributionGeometric.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlUserInterface.h"
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
            throw RbException() << "Could not find file " << filename;
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
 * Apply the general function to the current sampler object.
 */
double MarginalLikelihoodEstimator::marginalLikelihood( void ) const
{
    double out = marginalLikelihoodGeneral(powers, likelihoodSamples);
    return out;
}


/**
 * This function is adapted from TraceNumeric::update(): the only difference is that it works directly on a vector of numeric
 * values representing recorded samples, as opposed to a Trace object. The calculation is equivalent to the one used in
 * Tracer, and also (except for the choice of the maximum lag value) to convenience:::essTracer().
 *
 * @param values Vector of numeric parameter values
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


std::int64_t MarginalLikelihoodEstimator::offsetModulo(std::int64_t i, std::int64_t n) const
{
    // we do not add 1 to the result since we are using 0-based indexing, not 1-based indexing like in R
    return (i - 1) % n;
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
std::vector<std::int64_t> MarginalLikelihoodEstimator::getIndices(std::pair<std::int64_t, std::int64_t> a, std::int64_t n) const
{
    std::vector<std::int64_t> out;
    if (a.second != 0)
    {
        for (size_t i = 0; i < a.second; i++)
        {
            out.push_back( offsetModulo(a.first + i, n) );
        }
    }
    
    return out;
}


/**
 * Generate block bootstrap replicates of log likelihoods sampled by the power posterior analysis. This is done using the stationary
 * boostrap method of Politis & Romano (1994; J. Am. Stat. Assoc. 89(428): 1303--1313). The code below combines the functions
 * mcmc3r::block.boot() and mcmc3r:::ts.array() from Mario dos Reis' companion R package for MCMCTree (see
 * https://github.com/dosreislab/mcmc3r/blob/main/R/marginal-lhd.R), the latter of which was in turn lifted without
 * modifications from the eponymous function in the package "boot", written by Angelo J. Canty and corrected by B. D. Ripley (see
 * https://github.com/cran/boot/blob/master/R/bootfuns.q#L3160). Combining these two functions means that some
 * of the variables that were user-defined in ts.array() are here hardcoded to fixed values: specifically, n_sim (length of simulated
 * time series) is always equal to n (length of original time series), l (mean block length) is set to n * prop, the actual block length is
 * never fixed but rather always drawn from a geometric distribution with mean l, and end correction is always used.
 *
 * When print_to_file = true, new directories are created in the parent directory of the file storing the original log likelihoods for each
 * of the K powers / temperatures, called stone_1, stone_2, ..., stone_K. In each directory, text files are created for each of the
 * R = repnum bootstrap replicates, with basenames bootrep_1, bootrep_2, ..., bootrep_R, and the same extension as the file
 * storing the original log likelihoods. In addition, a file with basename bootrep_0 stores a copy of the original log-likelihood samples.
 *
 * @param repnum Number of bootstrap replicates
 * @param prop Mean block length, given as a proportion of the MCMC chain length n
 * @param print_to_file Should we create text files recording the resampled log likelihood values?
 * @return Nested vector of n log-likelihood samples (innermost) for each of the R replicates (inner) for each of the K powers (outer)
 */
std::vector< std::vector< std::vector<double> > > MarginalLikelihoodEstimator::blockBootstrap(size_t repnum, double prop, bool print_to_file) const
{
    // Get random number generator
    RandomNumberGenerator* rng = GLOBAL_RNG;
    
    // vector to be returned (repnum replicates for each of the k powers)
    std::vector< std::vector< std::vector<double> > > res( powers.size() );
    
    for (size_t i = 0; i < powers.size(); i++)
    {
        // vector of bootstrap replicates
        std::vector< std::vector<double> > res_inner( repnum + 1 );
        
        // the first element just corresponds to the original (non-resampled) vector of log likelihoods
        res_inner[0] = likelihoodSamples[i];
        
        // we keep n and n_sim distinct for conceptual reasons: technically, the simulated time series could contain a different
        // number of samples than the original time series
        size_t n = likelihoodSamples[i].size();
        std::int64_t n_sim = n;
        
        for (size_t j = 1; j < repnum + 1; j++)
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
                bool cond = false;
                
                for (size_t k = 0; k < repnum; k++)
                {
                    temp0[k] = 1 + RbStatistics::Geometric::rv( 1/(n * prop), *rng );
                    temp1[k] = std::min( temp0[k], n_sim - len_tot[k] );
                    len_tot[k] = len_tot[k] + temp1[k];
                    cond |= len_tot[k] < n_sim;
                }
                
                lens.push_back( temp1 );
                cont = cond;
            }
            
            std::vector<std::int64_t> inds_flat; // flat vector of indices
            
            for (size_t k = 0; k < lens.size(); k++)
            {
                std::vector<std::int64_t> elem( repnum );
                for(size_t l = 0; l < repnum; l++)
                {
                    // transform uniform01() to integer in range [1, endpt]
                    elem[l] = static_cast<std::int64_t>( std::floor(rng->uniform01() * endpt) ) + 1;
                }
                
                std::pair<std::int64_t, std::int64_t> ends( elem[0], lens[k][0] );
                std::vector<std::int64_t> inner = getIndices( ends, static_cast<std::int64_t>(n) );
                
                // add the vector resulting from this iteration (whose size is stochastic, and therefore unknown in advance) to our indices
                inds_flat.insert(inds_flat.end(), inner.begin(), inner.end());
            }
            
            // fill in the output vector; if the total length of the flat index vector exceeds n_sim, truncate
            for (size_t k = 0; k < n_sim; k++)
            {
                size_t idx = inds_flat[k];
                res_innermost[k] = likelihoodSamples[i][idx];
            }
            
            res_inner[j] = res_innermost;
        }
        
        res[i] = res_inner;
    }
    
    if (print_to_file)
    {
        path parent_dir = filename.parent_path();
        
        for (size_t i = 0; i < powers.size(); i++)
        {
            // create stone-specific directory
            std::string stone_dir_name = "stone_" + std::to_string(i + 1); // number the stones from 1 to n, not from 0 to (n - 1)
            path target_dir = parent_dir / stone_dir_name;
            create_directory(target_dir);
            
            for (size_t j = 0; j < repnum + 1; j++)
            {
                // create bootstrap replicate file
                std::string rep_name = "bootrep_" + std::to_string(j) + filename.extension().string();
                path rep_file_path = target_dir / rep_name;
                
                std::ofstream outStream( rep_file_path.string() );
                outStream << "power\t" << "likelihood" << std::endl;
                
                for (size_t k = 0; k < res[i][j].size(); k++)
                {
                    outStream << powers[i] << "\t" << res[i][j][k] << std::endl;
                }
                
                outStream.close();
            }
        }
    }
    
    return res;
}


/**
 * Generate block bootstrap replicates using blockBootstrap() (see above), estimate marginal likelihood on them, and use the
 * resulting values to place a standard error (and optionally, a credible interval) on the original estimate.
 *
 * @param repnum Number of bootstrap replicates
 * @param prop Mean block length, given as a proportion of the MCMC chain length n
 * @param print Should we create text files recording the resampled log likelihood values?
 * @param verbose Should we print the following the screen? (1) Vector of marginal likelihood estimates on the boostrap replicates and (2) 95% credible interval calculated from these replicates
 * @return Standard error of the marginal likelihood estimate calculated from the boostrap replicates
 */
double MarginalLikelihoodEstimator::standardErrorBlockBootstrap(size_t repnum, double prop, bool print, bool verbose) const
{
    std::vector< std::vector< std::vector<double> > > bootreps = blockBootstrap(repnum, prop, print);
    std::vector<double> marg_lnl_estimates( repnum + 1 );
    
    for (size_t i = 0; i < repnum + 1; i++)
    {
        std::vector< std::vector<double> > boot_rep( bootreps.size() );
        
        for (size_t j = 0; j < bootreps.size(); j++)
        {
            boot_rep[j] = bootreps[j][i];
        }
        
        marg_lnl_estimates[i] = marginalLikelihoodGeneral(powers, boot_rep);
    }
    
    // calculate the variance of the bootstraped estimates (not including the estimate derived from original samples)
    double mean_bootreps = 0.0;
    double var_bootreps = 0.0;
    
    for (size_t i = 1; i < repnum + 1; i++)
    {
        mean_bootreps += marg_lnl_estimates[i];
    }
    mean_bootreps /= repnum;
    
    for (size_t i = 1; i < repnum + 1; i++)
    {
        var_bootreps += (marg_lnl_estimates[i] - mean_bootreps)*(marg_lnl_estimates[i] - mean_bootreps);
    }
    var_bootreps /= (int(repnum) - 1);
    
    if (verbose)
    {
        std::stringstream ss;
        ss << "Bootstrapped marginal likelihood estimates:\n\n" << std::fixed << std::setprecision(3);
        
        // print 10 values per line
        size_t iter_num = floor(repnum / 10);
        for (size_t i = 1; i <= iter_num; i++)
        {
            size_t start_idx = 1 + (i - 1) * 10;
            size_t end_idx = i * 10;
            
            for (size_t j = start_idx; j < end_idx; j++)
            {
                ss << marg_lnl_estimates[j] << ", ";
            }
            
            if (end_idx == repnum) {
                ss << marg_lnl_estimates[end_idx] << "\n\n";
            }
            else
            {
                ss << marg_lnl_estimates[end_idx] << ",\n";
            }
        }
        
        if ((double)iter_num != (double)repnum / 10)
        {
            for (size_t i = 1 + iter_num * 10; i < repnum; i++)
            {
                ss << marg_lnl_estimates[i] << ", ";
            }
            
            ss << marg_lnl_estimates[repnum] << "\n\n";
        }
        
        ss << "Original estimate: " << marg_lnl_estimates[0];
        
        // get the 95% credible interval using a crude quantile function
        std::sort(marg_lnl_estimates.begin(), marg_lnl_estimates.end());
        size_t idx_lower = round( (marg_lnl_estimates.size() - 1) * 0.025 );
        size_t idx_upper = round( (marg_lnl_estimates.size() - 1) * 0.975 );
        
        ss << "; 95% CI: [" << marg_lnl_estimates[idx_lower] << ", " << marg_lnl_estimates[idx_upper] << "]";
        
        RBOUT( ss.str() );
    }
    
    // standard error = standard deviation of bootstrap replicates = square root of their variance
    return sqrt(var_bootreps);
}
