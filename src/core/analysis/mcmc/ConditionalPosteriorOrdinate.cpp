#include "ConditionalPosteriorOrdinate.h"

#include <stddef.h>
#include <cmath>
#include <vector>
#include <map>

#include "Cloneable.h"
#include "RbException.h"
#include "RbFileManager.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"

#ifdef RB_MPI
#include <mpi.h>
#endif

using namespace RevBayesCore;


ConditionalPosteriorOrdinate::ConditionalPosteriorOrdinate(const std::string &fn, const std::string &del, const std::vector<std::string>& skip_col_names)
{
    
    if ( process_active == true )
    {
        
        // check that the file/path name has been correctly specified
        RbFileManager my_file_manager( fn );
        if ( !my_file_manager.testFile() || !my_file_manager.testDirectory() )
        {
            std::string error_str = "";
            my_file_manager.formatError( error_str );
            throw RbException(error_str);
        }
            
            
        long thinning = 1;
        
        bool has_header_been_read = false;
            
        /* Open file */
        std::ifstream input_file( my_file_manager.getFullFileName().c_str() );
            
        if ( !input_file )
        {
            throw RbException( "Could not open file \"" + fn + "\"" );
        }
        
        /* Initialize */
        std::string command_line;
        RBOUT("Processing file \"" + fn + "\"");
        size_t n_samples = 0;
        
        // flags which columns to skip
        // we skip columns like the iteration number, posterior, likelihood and prior, or any other user defined column
        std::vector<bool> skip_column;
                
        /* Command-processing loop */
        while ( input_file.good() )
        {
                    
            // Read a line
            std::string line;
            my_file_manager.safeGetline(input_file, line);
                    
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


double ConditionalPosteriorOrdinate::predictiveProbability( const std::vector<double>& counts, bool site_probs_as_log ) const
{
    
    double pred_prob = 0.0;
    
    if ( process_active == true )
    {
    
        for (size_t i = 0; i < samples.size(); ++i)
        {
        
            std::vector<double> lnl_per_site = samples[i];
            size_t num_mcmc_samples = lnl_per_site.size();
            
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
            
            // get the number of observations for this site (if there are multiple times we observed this site)
            double count = 1.0;
            if (  counts.size() > 0 )
            {
                count = counts[i];
            }
        
            // add the per site predictive log probability
            pred_prob += ( min + log(num_mcmc_samples) - log(sum_per_site_probs) ) * count;
        
        }

    }
    
#ifdef RB_MPI
    MPI_Bcast(&pred_prob, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    return pred_prob;
}

