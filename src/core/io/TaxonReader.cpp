#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <set>
#include <map>
#include <algorithm>
#include <string>
#include <vector>

#include "RbException.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h" // for RBOUT
#include "StringUtilities.h"
#include "TaxonReader.h"
#include "DelimitedDataReader.h"
#include "Taxon.h"
#include "TimeInterval.h"

using namespace RevBayesCore;


/**
 * Parse a numeric field, requiring the whole field to be consumed.
 *
 * Deliberately not `stringstream >> double`, which yields 0 for a field that does not parse
 * and reports it only through failbit, and which does not recognize "Inf" -- a taxon file may
 * use that for an occurrence whose age has no upper bound.
 *
 * \param[in]    s        The raw field.
 * \param[in]    field    Field name, for the error message.
 * \param[in]    line     1-based line number in the file, for the error message.
 */
double TaxonReader::parseNumericField(const std::string &s, const std::string &field, size_t line)
{
    const char *begin = s.c_str();
    char *end = NULL;

    double v = std::strtod(begin, &end);

    // strtod skips leading whitespace; skip trailing before checking for leftovers
    while ( *end == ' ' || *end == '\t' || *end == '\r' || *end == '\n' ) ++end;

    if ( end == begin || *end != '\0' )
    {
        throw RbException() << "Could not parse \'" << s << "\' as a number in the \"" << field
                            << "\" field on line " << line << " of the taxon definition file.";
    }

    return v;
}


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

    // rows declaring an occurrence but counting zero of it
    std::vector<std::string> zero_count_rows;

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
            double age = parseNumericField( line[ column_map["age"] ], "age", i+1 );

            TimeInterval interval(age,age);

            if ( found == false )
            {
                taxon.setAgeRange(interval);
            }

            taxon.addOccurrence(interval);
        }

        if ( minit != column_map.end() )
        {
            TimeInterval interval;

            double min_age = parseNumericField( line[ column_map["min_age"] ], "min_age", i+1 );
            double max_age = parseNumericField( line[ column_map["max_age"] ], "max_age", i+1 );

            // ordering (and its intentional 1e-6 tolerance) is TimeInterval's to enforce
            interval.setMin(min_age);
            interval.setMax(max_age);

            if ( found == false )
            {
                taxon.setAgeRange(interval);
            }

            taxon.addOccurrence(interval);

            if ( countit != column_map.end() )
            {
                double c = parseNumericField( line[ column_map["count"] ], "count", i+1 );

                // 0 is meaningful: an extant species may have no fossil samples
                if ( c < 0.0 || RbMath::isFinite(c) == false || c != std::floor(c) )
                {
                    throw RbException() << "count (" << line[ column_map["count"] ] << ") must be a non-negative whole number on line "
                                        << i+1 << " of the taxon definition file.";
                }

                size_t k = size_t(c);

                // count 0 fits an extant species; on a row placing an occurrence in the
                // past it contradicts itself
                if ( k == 0 && max_age > 0.0 )
                {
                    std::stringstream ss;
                    ss << "\"" << taxon_name << "\" on line " << i+1;
                    zero_count_rows.push_back( ss.str() );
                }

                for(size_t j = 1; j < k; j++)
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
            // only an occurrence at the present proves survival; min_age == 0 can come
            // from a bin edge instead
            bool sampled_at_present = false;
            std::map<TimeInterval, size_t> occs = taxon.getOccurrences();
            for ( std::map<TimeInterval, size_t>::const_iterator it = occs.begin(); it != occs.end(); it++ )
            {
                if ( it->first.getMin() == 0.0 && it->first.getMax() == 0.0 )
                {
                    sampled_at_present = true;
                }
            }
            taxon.setExtinct( sampled_at_present == false );
        }
    }

    if ( zero_count_rows.empty() == false )
    {
        // list a few and count the rest, so the message stays on one line
        size_t shown = std::min( zero_count_rows.size(), size_t(2) );

        std::stringstream ss;
        ss << "Warning: count = 0 but max_age > 0 in the taxon file, row ";
        for (size_t j = 0; j < shown; j++)
        {
            ss << ( j > 0 ? ", " : "" ) << zero_count_rows[j];
        }
        if ( zero_count_rows.size() > shown )
        {
            ss << " (and " << zero_count_rows.size() - shown << " more)";
        }
        ss << ".";

        RBOUT( ss.str() );
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
