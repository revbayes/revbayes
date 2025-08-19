#include <cstdlib>
#include <iosfwd>
#include <vector>

#include "RbFileManager.h"
#include "TraceContinuousReader.h"
#include "DelimitedDataReader.h"
#include "Trace.h"
#include "TraceNumeric.h"


using namespace RevBayesCore;


/**
 * Constructor. Here we read in immediately the file and then we parse through each line
 * and extract the TraceContinuous information.
 *
 * \param[in]     fn       The name of the file where the data is stored.
 * \param[in]     delim    The delimiter between the columns.
 */
TraceContinuousReader::TraceContinuousReader(const path &fn, std::string delim) : DelimitedDataReader( fn, delim )
{
    
    const std::vector<std::string>& headers = chars[0];
    
    size_t numSkippedCols = 1;
    
    for (size_t j=numSkippedCols; j<headers.size(); j++)
    {
        RevBayesCore::TraceNumeric t;
        
        std::string parmName = headers[j];
        t.setParameterName(parmName);
        t.setFileName(fn);
        
        data.push_back( t );
    }
    
    for (size_t i = 1; i < chars.size(); ++i) //going through all the lines
    {
        const std::vector<std::string>& columns = chars[i];
        
        // adding values to the Tracess
        for (size_t j=numSkippedCols; j<columns.size(); j++)
        {
            TraceNumeric& t = data[j-numSkippedCols];
            std::string tmp = columns[j];
            double d = atof( tmp.c_str() );
            t.addObject(d);
        }
        
    }

    
}


/**
 * Get the TraceContinuous information read from the file.
 *
 * \return The vector of taxa.
 */
std::vector<TraceNumeric>& TraceContinuousReader::getTraces( void )
{
    
    return data;
}


/**
 * Get the TraceContinuous information read from the file.
 *
 * \return The vector of taxa.
 */
const std::vector<TraceNumeric>& TraceContinuousReader::getTraces( void ) const
{
    
    return data;
}
