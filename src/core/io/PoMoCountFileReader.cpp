#include <cstddef>
#include <string>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

#include "DiscreteTaxonData.h"
#include "PoMoCountFileReader.h"
#include "PoMoState4.h"
#include "PoMoState.h"
#include "Cloneable.h"
#include "DelimitedDataReader.h"
#include "HomologousDiscreteCharacterData.h"
#include "NaturalNumbersState.h"
#include "RbException.h"
#include "StringUtilities.h"
#include "Taxon.h"


using namespace RevBayesCore;


PoMoCountFileReader::PoMoCountFileReader(const path &fn, const size_t virtualPopulationSize, std::string d, size_t ns) : DelimitedDataReader(fn, d, ns), virtualPopulationSize_ ( virtualPopulationSize )
{
	matrix_ = new HomologousDiscreteCharacterData<PoMoState> ();

	// chars is a matrix containing all the lines of the file fn.
	// First line, with the names of the columns:
	// First line should be like: COUNTSFILE  NPOP 5   NSITES N

	int start = -1;
	// Skip comments.
	do {
			start = start + 1;
	}
	while (chars[start][0] == "#");

	if (chars[start][0] != "COUNTSFILE" || chars[0].size() != 5) {
            throw RbException()<<"File "<<fn<<" is not a proper PoMo Counts file: first line is not correct, it should be similar to \nCOUNTSFILE NPOP 5 NSITES N\n.";
	}
	else {
        numberOfPopulations_ = StringUtilities::asIntegerNumber( chars[0][2] );
        numberOfSites_ = StringUtilities::asIntegerNumber( chars[0][4] );
	}
	size_t numberOfFields = 2 + numberOfPopulations_;

	// The second line should look like this:
	//CHROM  POS  Sheep    BlackSheep  RedSheep  Wolf     RedWolf
	if (chars[start+1][0] != "CHROM" || chars[1][1] != "POS" || chars[1].size() != numberOfFields)
        {
            throw RbException()<<"File "<<fn<<" is not a proper PoMo Counts file: second line is not correct, it should be similar to \nCHROM POS Sheep BlackSheep RedSheep Wolf RedWolf\n.";
	}
	else
        {
            for (size_t i = start+2; i < 2 + numberOfPopulations_; ++i )
            {
                names_.push_back(chars[1][i]);
            }
	}

	// Setting the taxon names in the data matrix
	std::map<std::string, DiscreteTaxonData<PoMoState> > nameToTaxonData;
	for (size_t i = 0; i < names_.size(); ++i )
    {
		DiscreteTaxonData<PoMoState> tax (names_[i]);
		nameToTaxonData.insert(std::pair< std::string, DiscreteTaxonData<PoMoState> >(names_[i], tax) );
	}


	// Reading the other lines to fill up the taxon data
	//std::vector<std::string> chromosomes (numberOfSites_, "") ;
	//std::vector<size_t> positions (numberOfSites_, 0) ;
	//std::vector<std::string> sites (numberOfSites_, "") ;
	//std::map<std::string, std::vector<std::string> > speciesToSites;
	// for (size_t j = 0; j < numberOfPopulations_; ++j)
	// {
	// 	speciesToSites[names_[j]] = sites;
	// }

    std::vector<std::string> tmp_res;
    StringUtilities::stringSplit(chars[2][2], ",", tmp_res);
    size_t num_states = tmp_res.size();
    
	for (size_t i = 2; i < chars.size(); ++i)
	{
            if (chars[i].size() != numberOfFields)
            {
                throw RbException()<<"File "<<fn<<" is not a proper PoMo Counts file: line "<<i<<" is not correct, it does not have "<<numberOfFields<<" space-separated fields.";
            }

		//chromosomes.push_back(chars[i][0]);
    //positions.push_back( StringUtilities::asIntegerNumber( chars[i][1] ) );

		std::string chromosome = chars[i][0];
		size_t position = StringUtilities::asIntegerNumber( chars[i][1] );
		std::cout << "POSITION: "<< chars[i][1] << std::endl;

		for (size_t j = 2; j < 2 + numberOfPopulations_; ++j)
		{

			std::cout << "State: " << chars[i][j] <<std::endl;
            if ( num_states == 4 )
            {
                PoMoState4 pState (chars[i][j], chromosome, position, virtualPopulationSize_ );
                nameToTaxonData.at(names_[j-2]).addCharacter( pState);
            }
            else
            {
                PoMoState pState (num_states, virtualPopulationSize_, chars[i][j], chromosome, position );
                nameToTaxonData.at(names_[j-2]).addCharacter( pState);
            }
		}
	}

	// We have finished all lines, we fill up the data matrix
	for (std::map<std::string, DiscreteTaxonData<PoMoState> >::iterator tax = nameToTaxonData.begin(); tax != nameToTaxonData.end(); ++tax )
    {
	 	matrix_->addTaxonData(tax->second);
	}

	return ;
}

const size_t PoMoCountFileReader::getNumberOfPopulations( void )
{
	return numberOfPopulations_;
}

const size_t PoMoCountFileReader::getNumberOfSites( void )
{
	return numberOfSites_;
}

HomologousDiscreteCharacterData<PoMoState>* PoMoCountFileReader::getMatrix( void )
{

	return matrix_;

}

const size_t PoMoCountFileReader::getVirtualPopulationSize( void )
{
	return virtualPopulationSize_;
}
