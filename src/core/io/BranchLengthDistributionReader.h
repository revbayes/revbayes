#ifndef BranchLengthDistributionReader_H
#define BranchLengthDistributionReader_H

#include <cstddef>
#include <map>
#include <vector>
#include <iosfwd>
#include <utility>
#include <string>

#include "DelimitedDataReader.h"

namespace RevBayesCore {
	
	
	/**
	 * Reader for branch length distributions.
	 * Each line corresponds to a particular species tree branch length, defined using a pair of tip names.
	 * After the tip names, there is a list of branch lengths annotated with gene family name.
	 *
	 *
	 *
	 * @copyright Copyright 2009-
	 * @author The RevBayes Development Core Team (Bastien Boussau)
	 * @since 2015-03-03, version 1.0
	 *
	 */
	class BranchLengthDistributionReader : public DelimitedDataReader {
		
		public:
		
		BranchLengthDistributionReader(const std::string &fn, std::string d="", size_t ns=0);
		
		const std::vector< std::pair<std::string, std::string > >&     		 getTipPairs(void);
		const std::map < std::string, std::vector< double > >&   getBranchLengths(std::pair<std::string, std::string >& tipPair );
		const size_t 											 getNumberOfGeneFamilies(void);
		
		
		protected:
		
		std::vector< std::pair < std::string, std::string > >            TipPairs;
		std::vector < std::map < std::string, std::vector< double > > >  matrix;
		std::vector< std::string >                 		 				 geneFamilies;
		std::map< std::pair<std::string, std::string >, size_t >  		 pairToIndex;
		
	};
	
}

#endif
