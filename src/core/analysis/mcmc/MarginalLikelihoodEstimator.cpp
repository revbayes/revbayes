#include "MarginalLikelihoodEstimator.h"

#include <cstdlib>
#include <ostream>
#include <string>

#include "RbException.h"
#include "RbFileManager.h"
#include "RlUserInterface.h"
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
            throw RbException()<< "Could not open file " << fn;
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
    
        RBOUT("Processing file \"" + fn.string() + "\"");
    
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

