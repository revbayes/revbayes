#include <cstddef>
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
TaxonReader::TaxonReader(const std::string &fn, std::string delim) : DelimitedDataReader( fn, delim )
{
    
    //Reading the header
    std::vector<std::string>& line = chars[0];
    std::map<std::string, int> column_map;

    std::string arr[] = {"taxon", "species", "age", "min_age", "max_age", "status", "count"};
    std::vector<std::string> fields (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    
    for (size_t i = 0 ; i < line.size() ; ++i)
    {
        std::string tmp = line[i];
        StringUtilities::toLower( tmp );
        if ( std::find(fields.begin(), fields.end(), tmp) != fields.end())
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
            throw RbException() << "Unrecognized field: \'" << tmp << "\' in the taxon definition file. Allowed fields: " << field_stream.str(); 
        }
    }
    
    std::map<std::string,int>::iterator speciesit = column_map.find("species");
    std::map<std::string,int>::iterator ageit = column_map.find("age");
    std::map<std::string,int>::iterator minit = column_map.find("min_age");
    std::map<std::string,int>::iterator maxit = column_map.find("max_age");
    std::map<std::string,int>::iterator statusit = column_map.find("status");
    std::map<std::string,int>::iterator countit = column_map.find("count");

    if ( column_map.find("taxon") == column_map.end())
    {
        throw RbException("Taxon definition file must contain \"taxon\" field.");
    }
    if ( ageit == column_map.end() && minit == column_map.end())
    {
        throw RbException("Taxon definition file header must contain either \"age\" or (\"min_age\" and \"max_age\") fields");
    }
    if ( minit != column_map.end() && maxit == column_map.end() )
    {
        throw RbException("Taxon definition file header containing a \"min_age\" age field must also contain a \"max_age\" age field");
    }
    if ( maxit != column_map.end() && minit == column_map.end() )
    {
        throw RbException("Taxon definition file header containing a \"max_age\" age field must also contain a \"min_age\" age field");
    }
    if ( (minit != column_map.end() || maxit != column_map.end()) && ageit != column_map.end())
    {
        throw RbException("Taxon definition file header cannot contain both \"age\" and (\"min_age\" or \"max_age\") fields");
    }

    std::map<std::string, Taxon > taxon_map;

    for (size_t i = 1; i < chars.size(); ++i) //going through all the lines
    {
        const std::vector<std::string>& line = chars[i];
        if (line.size() != column_map.size())
        {
            std::stringstream err;
            err << "Line " << i+1 << " in taxon definition file does not contain "<<column_map.size()<<" elements";
            throw RbException(err.str());
        }
        std::string taxon_name = line[ column_map["taxon"] ];
        std::string species_name = taxon_name;

        bool found = ( taxon_map.find(taxon_name) != taxon_map.end() );

        Taxon& taxon = taxon_map[taxon_name];

        if ( speciesit != column_map.end() )
        {
            species_name = line[ column_map["species"] ];

            if ( found == true && species_name != taxon.getSpeciesName() )
            {
                std::stringstream ss;
                ss << "Inconsistent species name for taxon \"" << taxon_name << "\"";
                throw(RbException(ss.str()));
            }
        }

        if ( found == false )
        {
            taxon.setName(taxon_name);
            taxon.setSpeciesName(species_name);
        }
        
        if ( ageit != column_map.end() )
        {
            double age = 0.0;
            std::stringstream ss;
            ss.str( line[ column_map["age"] ] );
            ss >> age;

            TimeInterval interval(age,age);

            if ( found == false )
            {
                taxon.setAgeRange(interval);
            }

            taxon.addOccurrence(interval);
        }

        if ( minit != column_map.end() )
        {
            double min_age, max_age;
            TimeInterval interval;
            std::stringstream ss;

            ss.str( line[ column_map["min_age"] ] );
            ss >> min_age;
            ss.clear();

            ss.str( line[ column_map["max_age"] ] );
            ss >> max_age;
            ss.clear();

            interval.setMin(min_age);
            interval.setMax(max_age);

            if ( found == false )
            {
                taxon.setAgeRange(interval);
            }

            taxon.addOccurrence(interval);

            if ( countit != column_map.end() )
            {
                size_t k = 0;
                std::stringstream ss;

                ss.str( line[ column_map["count"] ] );
                ss >> k;

                for(size_t i = 1; i < k; i++)
                {
                    taxon.addOccurrence(interval);
                }
            }
        }

        if ( statusit != column_map.end() )
        {
            bool extinct = (line[ column_map["status"] ] == "extinct");

            if ( found == true && extinct != taxon.isExtinct() )
            {
                std::stringstream ss;
                ss << "Inconsistent extinction status for taxon \"" << taxon_name << "\"";
                throw(RbException(ss.str()));
            }

            taxon.setExtinct( extinct );
        }
        else
        {
            taxon.setExtinct( taxon.getMinAge() > 0.0 );
        }
    }

    for (std::map<std::string, Taxon>::iterator it = taxon_map.begin(); it != taxon_map.end(); it++ )
    {
        taxa.push_back(it->second);
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
