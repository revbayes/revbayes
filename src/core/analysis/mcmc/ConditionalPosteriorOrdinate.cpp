#include "ConditionalPosteriorOrdinate.h"

#include <stddef.h>
#include <cmath>
#include <vector>
#include <map>

#include "Cloneable.h"
#include "RbException.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "RbConstants.h"
#include "RbMathLogic.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;


ConditionalPosteriorOrdinate::ConditionalPosteriorOrdinate(const path& fn, const std::string& del, const std::vector<std::string>& skip_col_names)
{
    
    if ( process_active == true )
    {
        /********************/
        /* read in the file */
        /********************/
    
        if ( not is_regular_file(fn) )
        {
            throw RbException()<< "Could not open file " << fn;
        }
            
            
        long thinning = 1;
        
        bool has_header_been_read = false;
        
        // Open file
        std::ifstream input_file( fn.string() );
        
        if ( !input_file )
        {
            throw RbException()<<"Could not open file "<<fn;
        }
        
        /* Initialize */
        std::string command_line;
        RBOUT("Processing file \"" + fn.string() + "\"");
        size_t n_samples = 0;
        
        // flags which columns to skip
        // we skip columns like the iteration number, posterior, likelihood and prior, or any other user defined column
        std::vector<bool> skip_column;
                
        /* Command-processing loop */
        while ( input_file.good() )
        {
                    
            // Read a line
            std::string line;
            safeGetline( input_file, line );
                    
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
            if ( has_header_been_read == false )
            {
                
                // re-initialize the flags vector which columns to skip
                skip_column = std::vector<bool>( columns.size(), false );
                        
                for (size_t j=0; j<columns.size(); j++)
                {
                    // check if we should skip this column
                    for (size_t n=0; n<skip_col_names.size(); ++n)
                    {
                        if ( skip_col_names[n] == columns[j] )
                        {
                            skip_column[j] = true;
                            break;
                        }
                    }
                    
                    // only add if we don't skip this column
                    if ( skip_column[j] == false )
                    {
                        std::vector<double> t;
                        samples.push_back( t );
                    }
                    
                }
                        
                has_header_been_read = true;
                        
                continue;
            }
                
                
            // increase our sample counter
            ++n_samples;
                
            // we need to check if we skip this sample in case of thinning.
            if ( (n_samples-1) % thinning > 0 )
            {
                continue;
            }
                
            // adding values to the Tracess
            size_t column_index = 0;
            for (size_t j=0; j<columns.size(); j++)
            {
                // only if we don't skip this column
                if ( skip_column[j] == false )
                {
                    std::vector<double>& t = samples[column_index];
                    std::string tmp = columns[j];
                    double d = atof( tmp.c_str() );
                    t.push_back(d);
                    
                    ++column_index;
                }
            }

        } // end-while over the input file content

    } // end-if this process is active
    
}



ConditionalPosteriorOrdinate::~ConditionalPosteriorOrdinate()
{
    
}



ConditionalPosteriorOrdinate* ConditionalPosteriorOrdinate::clone( void ) const
{
    return new ConditionalPosteriorOrdinate(*this);
}


double ConditionalPosteriorOrdinate::predictiveProbability( const std::vector<double>& counts, bool site_probs_as_log, bool folded ) const
{
    
    double pred_prob = 0.0;
    
    if ( process_active == true )
    {
        
        std::vector<double> cpo     = std::vector<double>( samples.size(), 0.0 );
        std::vector<bool>   zero    = std::vector<bool>( samples.size(), true );

        for (size_t i = 0; i < samples.size(); ++i)
        {
        
            std::vector<double> lnl_per_site = samples[i];
            size_t num_mcmc_samples = lnl_per_site.size();
            
            // check if all probabilities are zero
            for (size_t j = 0; j < num_mcmc_samples; ++j)
            {
                if ( site_probs_as_log )
                {
                    zero[i] = zero[i] & RbMath::isFinite( lnl_per_site[j] );
                }
                else
                {
                    zero[i] = zero[i] & (lnl_per_site[j] == 0);
                }
            }
            
            if ( site_probs_as_log == false )
            {
                for (size_t j = 0; j < num_mcmc_samples; ++j)
                {
                    lnl_per_site[j] = log(lnl_per_site[j]);
                }
            }
            
            // first, we need to find the minimum ln(p) to make sure we don't get underflows
            double min = lnl_per_site[0];
            for (size_t j = 1; j < num_mcmc_samples; ++j)
            {
                if (min > lnl_per_site[j])
                {
                    min = lnl_per_site[j];
                }
            }
        
            // now we can sum the "normalized" ln(p)
            double sum_per_site_probs = 0.0;
            for (size_t j = 0; j < num_mcmc_samples; ++j)
            {
                sum_per_site_probs += exp( (min - lnl_per_site[j]) );
            }
        
            // add the per site predictive log probability
            cpo[i] = ( min + log(num_mcmc_samples) - log(sum_per_site_probs) );
        
        }
        
        size_t max_state = samples.size();
        if ( folded == true )
        {
            max_state = floor( samples.size()/2.0 + 1);
        }
        
        
        for (size_t i = 0; i < max_state; ++i)
        {
            
            
            // get the number of observations for this site (if there are multiple times we observed this site)
            double count = 1.0;
            if (  counts.size() > 0 )
            {
                count = counts[i];
            }
            
            double this_prob = 0.0;
            
            if ( folded == true )
            {
                if ( (2.0*(i+1)) < samples.size() )
                {
                    
                    if ( zero[i] == true && zero[samples.size()-i-1]  == true )
                    {
                        this_prob = RbConstants::Double::neginf;
                    }
                    else if ( zero[i] == true  )
                    {
                        this_prob =cpo[samples.size()-i-1];
                    }
                    else if ( zero[samples.size()-i-1]  == true )
                    {
                        this_prob = cpo[i];
                    }
                    else
                    {
                        this_prob = log( 1 + exp(cpo[i] - cpo[samples.size()-i-1]) ) + cpo[samples.size()-i-1];
                    }
                }
                else
                {
                    if ( zero[i] == true )
                    {
                        this_prob = 0.0;
                    }
                    else
                    {
                        this_prob = cpo[i];
                    }
                }
                
            }
            else
            {
                
                if ( zero[i] == true )
                {
                    this_prob = 0.0;
                }
                else
                {
                    this_prob = cpo[i];
                }
            }
        
            // add the per site predictive log probability
            if ( count > 0 )
            {
                pred_prob += this_prob * count;
            }
        }

    }
    
#ifdef RB_MPI
    MPI_Bcast(&pred_prob, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return pred_prob;
}

