#include "BranchLengthDistributionReader.h"
#include "RbFileManager.h"

#include "RbConstants.h"


using namespace RevBayesCore;


BranchLengthDistributionReader::BranchLengthDistributionReader(const std::string &fn, char d, size_t numSkipped) : DelimitedDataReader(fn, d, numSkipped)
{
	
	filename = fn;
	
	for (size_t i = 0; i < chars.size(); ++i)
	{
		//Get the two tip names:
		std::string name1 = chars[i][0];
		std::string name2 = chars[i][1];
		TipPairs.push_back(std::make_pair( name1, name2 ) );
		pairToIndex[std::make_pair( name1, name2 )] = i;
		//Now we get the branch lengths
		//They have the format familyName:length
		std::string delimiter = ":";

		for (size_t j = 2; j < chars[i].size(); ++j)
		{
            std::string s = chars[i][j];
			size_t delimiterPosition = s.find(delimiter) ;
			std::string familyName = s.substr(0, delimiterPosition);
//			double bl = atof ( s.substr(delimiterPosition + 1, s.size() ).c_str() );
			
		}
		
	}
}


const std::vector<std::pair<std::string, std::string > >& BranchLengthDistributionReader::getTipPairs(void)
{
	return TipPairs;
}


const std::map < std::string, std::vector< double > >&   BranchLengthDistributionReader::getBranchLengths( std::pair<std::string, std::string >& tipPair)
{
	size_t id = pairToIndex[ tipPair ];
	return matrix[id];
}


const size_t BranchLengthDistributionReader::getNumberOfGeneFamilies(void) {
	
	return geneFamilies.size();
}


