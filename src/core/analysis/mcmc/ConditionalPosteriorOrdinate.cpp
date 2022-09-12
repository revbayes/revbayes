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


ConditionalPosteriorOrdinate::ConditionalPosteriorOrdinate(const std::string &fn, const std::string &pn, const std::string &ln,  const std::string &del) 
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
            throw RbException( "Could not open file \"" + fn + "\"" );
                
        /* Initialize */
        std::string command_line;
        RBOUT("Processing file \"" + fn + "\"");
        size_t n_samples = 0;
                
        /* Command-processing loop */
        while ( input_file.good() )
        {
                    
            // Read a line
            std::string line;
            my_file_manager.safeGetline(input_file, line);
                    
            // skip empty lines
            //line = stringByTrimmingCharactersInSet:[NSCharacterSet whitespaceAndNewlineCharacterSet]];
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
                        
                for (size_t j=0; j<columns.size(); j++)
                {
                    std::vector<double> t;
                    samples.push_back( t );
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
            for (size_t j=0; j<columns.size(); j++)
            {
                std::vector<double>& t = samples[j];
                std::string tmp = columns[j];
                double d = atof( tmp.c_str() );
                t.push_back(d);
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


double ConditionalPosteriorOrdinate::marginalLikelihood( void ) const
{
    
    double marginal = 0.0;
    
    if ( process_active == true )
    {
    
//        std::vector<double> pathValues;
//        for (size_t i = 1; i < powers.size(); ++i)
//        {
//        
//            size_t samplesPerPath = likelihoodSamples[i].size();
//            double max = likelihoodSamples[i][0];
//            for (size_t j = 1; j < samplesPerPath; ++j)
//            {
//                if (max < likelihoodSamples[i][j])
//                {
//                    max = likelihoodSamples[i][j];
//                }
//            }
//        
//            // mean( exp(samples-max)^(beta[k-1]-beta[k]) )     or
//            // mean( exp( (samples-max)*(beta[k-1]-beta[k]) ) )
//            double mean = 0.0;
//            for (size_t j = 0; j < samplesPerPath; ++j)
//            {
//                mean += exp( (likelihoodSamples[i][j]-max)*(powers[i-1]-powers[i]) ) / samplesPerPath;
//            }
//        
//            marginal += log(mean) + (powers[i-1]-powers[i])*max;
//        
//        }

    }
    
#ifdef RB_MPI
    MPI_Bcast(&marginal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    
    
    return marginal;
}

