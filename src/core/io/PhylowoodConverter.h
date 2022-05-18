#ifndef PhylowoodConverter_H
#define PhylowoodConverter_H

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
        std::string                                 buildCharacterHistoryString(TopologyNode* n, unsigned end);
        std::string                                 buildExtendedNewick(TopologyNode* n);
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



#endif 
