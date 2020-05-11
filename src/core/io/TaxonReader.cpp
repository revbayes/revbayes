#include <stddef.h>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <vector>

#include "RbException.h"
#include "StringUtilities.h"
#include "TaxonReader.h"
#include "DelimitedDataReader.h"
#include "Taxon.h"
#include "TimeInterval.h"

using namespace RevBayesCore;


/**
 * Constructor. Here we read in immediately the file and then we parse through each line 
 * and extract the taxon information.
 *
 * \param[in]     fn       The name of the file where the data is stored.
 * \param[in]     delim    The delimiter between the columns.
 */
TaxonReader::TaxonReader(const std::string &fn, char delim) : DelimitedDataReader( fn, delim )
{
    
    //Reading the header
    std::vector<std::string>& line = chars[0];
    std::map<std::string, int> column_map;

    std::string arr[] = {"taxon","age","species","min","max"};
    std::vector<std::string> fields (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    
    for (size_t i = 0 ; i < line.size() ; ++i)
    {
        std::string tmp = line[i];
        StringUtilities::toLower( tmp );
        if (std::find(fields.begin(), fields.end(), tmp) != fields.end())
        {
            column_map[tmp] = int(i);
        }
        else
        {
            std::stringstream field_stream;
            for (size_t j = 0; j < fields.size(); j++)
            {
                field_stream << "\"" << fields[j] << "\"";
                if (j < fields.size() - 1)
                {
                    field_stream << ", ";
                }
            }
            throw RbException("Unrecognized field: "+tmp+" in the taxon definition file. Allowed fields: "+field_stream.str());
        }
    }
    
    if (column_map.find("taxon") == column_map.end())
    {
        throw RbException("Missing header in the taxon definition file. It has to contain \"taxon\" field.");
    }
    
    std::map<std::string,int>::iterator minit = column_map.find("min");
    std::map<std::string,int>::iterator maxit = column_map.find("max");

    if ( (minit == column_map.end() || maxit == column_map.end()) && minit != maxit)
    {
        throw RbException("Taxon definition file header must contain both \"min\" and \"max\" age fields");
    }
    if ( (minit != column_map.end() || maxit != column_map.end()) && column_map.find("age") != column_map.end())
    {
        throw RbException("Taxon definition file header cannot contain both \"age\" and (\"min\" or \"max\") fields");
    }

    for (size_t i = 1; i < chars.size(); ++i) //going through all the lines
    {
        const std::vector<std::string>& line = chars[i];
        if (line.size() != column_map.size())
        {
            std::stringstream err;
            err << "Line " << i+1 << " in taxon definition file does not contain "<<column_map.size()<<" elements";
            throw RbException(err.str());
        }
        Taxon t = Taxon( line[ column_map["taxon"] ] );
        
        if ( column_map.find("min") != column_map.end() )
        {
            double min = 0,max = 0;
            TimeInterval interval;
            std::stringstream ss;

            ss.str( line[ column_map["min"] ] );
            ss >> min;
            interval.setMin(min);

            ss.clear();
            ss.str( line[ column_map["max"] ] );
            ss >> max;
            interval.setMax(max);

            t.setAgeRange(interval);
        }
        else if ( column_map.find("age") != column_map.end() )
        {
            double age = 0.0;
            std::stringstream ss;
            ss.str( line[ column_map["age"] ] );
            ss >> age;
            t.setAge( age );
        }
        else
        {
            t.setAge( 0.0 );
        }

        
        if ( column_map.find("species") != column_map.end() )
        {
            t.setSpeciesName( line[ column_map["species"] ] );
        }
        
        taxa.push_back( t );
    }

    
    std::set<std::string> found;
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if (found.find(taxa[i].getName()) == found.end())
        {
            found.insert(taxa[i].getName());
        }
        else
        {
            std::stringstream ss;
            ss << "Duplicate taxon name '" << taxa[i].getName() << "' encountered when reading taxon definition file";
            throw(RbException(ss.str()));
        }
    }

}


/**
 * Get the taxon information read from the file.
 *
 * \return The vector of taxa.
 */
const std::vector<Taxon>& TaxonReader::getTaxa( void ) const
{
    
    return taxa;
}
