#include <cstddef>
#include <string>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "DiscreteTaxonData.h"
#include "PoMoCountFileReader.h"
#include "PoMoState.h"
#include "Cloneable.h"
#include "DelimitedDataReader.h"
#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Taxon.h"


using namespace RevBayesCore;


PoMoCountFileReader::PoMoCountFileReader(   const path &fn, 
                                            const size_t vps, 
                                            FORMAT f,
                                            const string &wm,
                                            const long eps ) : DelimitedDataReader(fn, ""),
    virtual_population_size ( vps ),
    data_format( f ),
    weighting_method( wm ),
    effective_population_size( eps )
{
    if ( data_format == FORMAT::PoMo )
    {
        matrix = new HomologousDiscreteCharacterData<PoMoState>();
    }
    else
    {
        matrix = new HomologousDiscreteCharacterData<NaturalNumbersState>();
    }

    RevBayesCore::PoMoState::WEIGHTING weighting = RevBayesCore::PoMoState::FIXED;
    if ( weighting_method == "Binomial" )
    {
        weighting  = RevBayesCore::PoMoState::BINOMIAL;
    } else if ( weighting_method == "Sampled" )
    {
        weighting  = RevBayesCore::PoMoState::SAMPLED;
    } else if ( weighting_method == "Hypergeometric" )
    {
        weighting  = RevBayesCore::PoMoState::HYPERGEOMETRIC;
    } else if ( weighting_method == "None" )
    {
        weighting  = RevBayesCore::PoMoState::NONE;
    }
    
	// chars is a matrix containing all the lines of the file fn.

	int start = -1;
	// ignoring comments
	do {
			start = start + 1;
	}
	while (chars[start][0] == "#");

    // checking whether the 1st header line of the count file is properly formatted
    // first header should be like: COUNTSFILE  NPOP 5   NSITES N
	if (chars[start][0] != "COUNTSFILE" || chars[start].size() != 5)
    {
        throw RbException()<<"File "<<fn<<" is not a proper PoMo counts file. First line is not properly formatted: it should be similar to \nCOUNTSFILE NPOP X NSITES Y\n.";
	}
	else
    {
        n_taxa  = StringUtilities::asIntegerNumber( chars[start][2] );
        n_sites = StringUtilities::asIntegerNumber( chars[start][4] );
	}
	size_t numberOfFields = 2 + n_taxa;

    // checking whether the 1st header line of the count file is properly formated
	// The second line should look like this:
	// CHROM  POS  Sheep    BlackSheep  RedSheep  Wolf     RedWolf
	if (chars[start+1][0] != "CHROM" || chars[start+1][1] != "POS" || chars[start+1].size() != numberOfFields)
    {
        throw RbException()<<"File "<<fn<<" is not a proper PoMo counts file. Second line is not properly formatted: it should be similar to \nCHROM POS Species 1 Species 2 ... SpeciesX\n.";
	}
	else
    {
		for (size_t i = 2; i < 2 + n_taxa; ++i )
        {
			names.push_back(chars[start+1][i]);
		}
	}

	// setting the taxon names in the data matrix
	std::map<std::string, DiscreteTaxonData<PoMoState> > name_to_taxon_data;
	for (size_t i = 0; i < names.size(); ++i )
    {
		DiscreteTaxonData<PoMoState> tax (names[i]);
		name_to_taxon_data.insert(std::pair< std::string, DiscreteTaxonData<PoMoState> >(names[i], tax) );
	}


    // estimate the number of states
    size_t row=start+2;
    size_t col=start+2;

    // move ahead until we found a non-missing state
    while ( chars[row][col] == "?" || chars[row][col] == "-" )
    {
        ++col;
        if ( col == chars[row].size() )
        {
            // reset at end of row
            ++row;
            col = 2;
        }
        if ( row == chars.size() )
        {
            throw RbException("Could not compute the number of alleles in PoMo counts file because all states are missing.");
        }
    }
    
    // calculating the number of alleles
    std::vector<std::string> tmp_res;
    StringUtilities::stringSplit(chars[row][col], ",", tmp_res);
    size_t n_alleles = tmp_res.size();

    // we have officially reached the counts part of the count file
	for (size_t i = start + 2; i < chars.size(); ++i)
	{
            
        if (chars[i].size() != numberOfFields)
        {
            throw RbException()<<"File "<<fn<<" is not a proper PoMo Counts file. Line "<<i<<" does not have "<<numberOfFields<<" space-separated fields.";
        }

        // getting chromosome and position
		const std::string& chromosome = chars[i][0];
		size_t position               = StringUtilities::asIntegerNumber( chars[i][1] );

        // reading the counts of a genomic position for each species 
		for (size_t j = 2; j < 2 + n_taxa; ++j)
		{
            // creating PoMo state
            PoMoState pState (n_alleles, virtual_population_size, chars[i][j], chromosome, position, weighting, effective_population_size );
            name_to_taxon_data.at(names[j-2]).addCharacter( pState);

            std::vector<double> vec = pState.getWeights();
            
            // // checks the weight
            // std::cout << " Count: " << pState.getStringValue() << "\n";
            // for (int i=0; i < vec.size(); i++)
            // {
            //     std::cout << vec[i] << ' ';
            // }
            // std::cout << "\n\n";
		}
	}

	// We have finished all lines, we fill up the data matrix
    for (std::map<std::string, DiscreteTaxonData<PoMoState> >::iterator tax = name_to_taxon_data.begin(); tax != name_to_taxon_data.end(); ++tax )
    {
        matrix->addTaxonData(tax->second);
    }


}

PoMoCountFileReader::PoMoCountFileReader(const PoMoCountFileReader& r) : DelimitedDataReader(r),
n_taxa(r.n_taxa),
n_sites(r.n_sites),
virtual_population_size(r.virtual_population_size),
names(r.names),
data_format(r.data_format)
{
    matrix = r.matrix->clone();
}


PoMoCountFileReader::~PoMoCountFileReader()
{
    delete matrix;
}

PoMoCountFileReader& PoMoCountFileReader::PoMoCountFileReader::operator=(const PoMoCountFileReader& r)
{
    if (this != &r)
    {
        n_taxa                  = r.n_taxa;
        n_sites                 = r.n_sites;
        virtual_population_size = r.virtual_population_size;
        names                   = r.names;
        data_format             = r.data_format;
        
        matrix                  = r.matrix->clone();
    }
    
    return *this;
}

const size_t PoMoCountFileReader::getNumberOfPopulations( void )
{
	return n_taxa;
}

const size_t PoMoCountFileReader::getNumberOfSites( void )
{
	return n_sites;
}

const AbstractHomologousDiscreteCharacterData& PoMoCountFileReader::getMatrix( void )
{

	return *matrix;

}

const size_t PoMoCountFileReader::getVirtualPopulationSize( void )
{
	return virtual_population_size;
}
