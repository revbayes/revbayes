//
//  PhylowoodConverter.h
//  revbayes-proj
//
//  Created by Michael Landis on 3/2/15.
//  Copyright (c) 2015 Michael Landis. All rights reserved.
//

#ifndef __revbayes_proj__PhylowoodConverter__
#define __revbayes_proj__PhylowoodConverter__

#include <string>
#include "TopologyNode.h"
#include "RbFileManager.h"

namespace RevBayesCore {
    
    class DelimitedDataReader;
    class TimeAtlas;
    class PhylowoodConverter {
        
    public:
        PhylowoodConverter(const path &sfn, const path &tfn, const path &gfn, const path &pfn, double b, const std::string& ct, const std::string& bt);
        ~PhylowoodConverter(void);
        void                                        convert(void);
        
    private: 
        std::string                                 buildCharacterHistoryString(const TopologyNode* n, unsigned end);
        std::string                                 buildExtendedNewick(const TopologyNode* n);
        std::string                                 buildPhylowoodString(void);
        void                                        makeMarginalAreaProbs(void);
        void                                        makeBits(void);
        void                                        test(void);

        Tree* tree;
        TimeAtlas* atlas;
        DelimitedDataReader* dat;
        
        std::vector<std::vector<unsigned> > bits;
        std::vector<std::vector<double> > marginalStartProbs;
        std::vector<std::vector<double> > marginalEndProbs;
        
        path stateFilename;
        path treeFilename;
        path geoFilename;
        path phwFilename;
        std::string chartype;
        std::string bgtype;
    
        double burn;
        size_t num_nodes;
        size_t num_states;
        size_t numAreas;
        size_t numEpochs;
    };
};



#endif /* defined(__revbayes_proj__PhylowoodConverter__) */
