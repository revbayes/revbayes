#include "MarginalLikelihoodEstimator.h"

#include <cstdlib>
#include <ostream>
#include <string>

#include "RbException.h"
#include "RbFileManager.h"
#include "StringUtilities.h"

using namespace RevBayesCore;


MarginalLikelihoodEstimator::MarginalLikelihoodEstimator(const path &fn, const std::string &pn, const std::string &ln, const std::string &del) :
    powers(),
    likelihoodSamples()
{
    
    setActivePID( 0, 1 );
    
    
    if ( process_active == true )
    {
        /********************/
        /* read in the file */
        /********************/
    
        if ( not is_regular_file(fn) )
        {
            throw RbException()<< "Could not not file " << fn;
        }
    
    
        bool hasHeaderBeenRead = false;
        
        // Open file
        std::ifstream inFile( fn.string() );
        
        if ( !inFile )
        {
            throw RbException()<<"Could not open file "<<fn;
        }
    
        // Initialize
        std::string commandLine;
    
        //    RBOUT("Processing file \"" + fn + "\"");
    
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
                throw RbException() << "Please check format of file " << fn << ", missing power and likelihood columns in some lines";

            double p, l;
            // get the power entry
            std::string tmp = columns[powerColumnIndex];
            try {
                p = std::stod(tmp);
            } catch (std::invalid_argument&) {
                throw RbException() << "Please check format of file " << fn << ", non-numeric input in power column";
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
                throw RbException() << "Please check format of file " << fn << ", non-numeric input in likelihood column";
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


size_t MarginalLikelihoodEstimator::offsetModulo(size_t i, size_t n)
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
std::vector<size_t> MarginalLikelihoodEstimator::getIndices(std::pair<size_t, size_t> a, size_t n)
{
    std::vector<size_t> out;
    if (a.second != 0)
    {
        std::vector<size_t> seq( a.second );
        for (size_t i = 0; i < a.second; i++)
        {
            seq[i] = a.first + i;
            out.push_back( offsetModulo(seq[i], n) );
        }
    }
    
    return out;
}
