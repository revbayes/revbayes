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

    std::string arr[] = {"taxon","age","species","min","max", "minmin", "minmax", "maxmin", "maxmax"};
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
        throw RbException("Taxon definition file must contain \"taxon\" field.");
    }
    
    std::map<std::string,int>::iterator speciesit = column_map.find("species");
    std::map<std::string,int>::iterator ageit = column_map.find("age");
    std::map<std::string,int>::iterator minit = column_map.find("min");
    std::map<std::string,int>::iterator maxit = column_map.find("max");
    std::map<std::string,int>::iterator minminit = column_map.find("minmin");
    std::map<std::string,int>::iterator minmaxit = column_map.find("minmax");
    std::map<std::string,int>::iterator maxminit = column_map.find("maxmin");
    std::map<std::string,int>::iterator maxmaxit = column_map.find("maxmax");

    if ( (minit != column_map.end() || minminit != column_map.end() || minmaxit != column_map.end())
            && maxit == column_map.end() && maxminit == column_map.end() && maxmaxit == column_map.end())
    {
        throw RbException("Taxon definition file header containing a \"min*\" age field must also contain a \"max*\" age field");
    }
    if ( (maxit != column_map.end() || maxminit != column_map.end() || maxmaxit != column_map.end())
            && minit == column_map.end() && minminit == column_map.end() && minmaxit == column_map.end())
    {
        throw RbException("Taxon definition file header containing a \"max*\" age field must also contain a \"min*\" age field");
    }
    if ( (minminit == column_map.end() || minmaxit == column_map.end()) && minminit != minmaxit)
    {
        throw RbException("Taxon definition file header must contain both \"minmin\" and \"minmax\" age fields");
    }
    if ( (maxminit == column_map.end() || maxmaxit == column_map.end()) && maxminit != maxmaxit)
    {
        throw RbException("Taxon definition file header must contain both \"maxmin\" and \"maxmax\" age fields");
    }
    if ( (minit != column_map.end() || maxit != column_map.end()
            || minminit != column_map.end() || minmaxit != column_map.end()
            || maxminit != column_map.end() || maxmaxit != column_map.end()) && ageit != column_map.end())
    {
        throw RbException("Taxon definition file header cannot contain both \"age\" and (\"min*\" or \"max*\") fields");
    }
    if ( (minit != column_map.end() || maxit != column_map.end()) &&
            (minminit != column_map.end() || minmaxit != column_map.end()
            || maxminit != column_map.end() || maxmaxit != column_map.end()))
    {
        throw RbException("Taxon definition file header cannot contain both \"min/max\" and (\"min/maxmin\" or \"min/maxmax\") fields");
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
        t.setAge( 0.0 );
        
        if ( ageit != column_map.end() )
        {
            double age = 0.0;
            std::stringstream ss;
            ss.str( line[ column_map["age"] ] );
            ss >> age;
            t.setAge( age );
        }

        if ( minit != column_map.end() || minminit != column_map.end() )
        {
            double age = 0;
            TimeInterval interval;
            std::stringstream ss;

            if ( minminit != column_map.end() )
            {
                ss.str( line[ column_map["minmin"] ] );
                ss >> age;
                interval.setMin(age);

                ss.clear();
                ss.str( line[ column_map["minmax"] ] );
                ss >> age;
                interval.setMax(age);
            }
            else if ( minit != column_map.end() )
            {
                ss.str( line[ column_map["min"] ] );
                ss >> age;
                interval.setMin(age);
                interval.setMax(age);
            }

            t.setMinAgeRange(interval);
        }

        if ( maxit != column_map.end() || maxminit != column_map.end() )
        {
            double age = 0;
            TimeInterval interval;
            std::stringstream ss;

            if ( maxminit != column_map.end() )
            {
                ss.str( line[ column_map["maxmin"] ] );
                ss >> age;
                interval.setMin(age);

                ss.clear();
                ss.str( line[ column_map["maxmax"] ] );
                ss >> age;
                interval.setMax(age);
            }
            else if ( maxit != column_map.end() )
            {
                ss.str( line[ column_map["max"] ] );
                ss >> age;
                interval.setMin(age);
                interval.setMax(age);
            }

            t.setMaxAgeRange(interval);
        }

        
        if ( speciesit != column_map.end() )
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
