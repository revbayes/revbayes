#ifndef TreeCharacterHistoryNhxMonitor_H
#define TreeCharacterHistoryNhxMonitor_H


#include "Monitor.h"
#include "BranchHistoryDiscrete.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"
#include "TimeAtlas.h"
#include "Tree.h"
#include "GeographicArea.h"
#include "AbstractHomologousDiscreteCharacterData.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace RevBayesCore {
    
    template<class charType>
    class TreeCharacterHistoryNhxMonitor : public Monitor {
        
    public:
        // Constructors and Destructors
        TreeCharacterHistoryNhxMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* s, TypedDagNode<Tree> *t, const TimeAtlas* ta, std::uint64_t g, std::uint64_t mg, int burn, const path &fname, const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false, bool sm=true, bool sr=true);
        
        TreeCharacterHistoryNhxMonitor(const TreeCharacterHistoryNhxMonitor& f);
        
        // basic methods
        TreeCharacterHistoryNhxMonitor*   clone(void) const;                                                  //!< Clone the object
        
        // Monitor functions
        void                                    monitor(std::uint64_t gen);                                                  //!< Monitor at generation gen
        void                                    swapNode(DagNode *oldN, DagNode *newN);
        
        // FileMonitor functions
        void                                    closeStream(void);                                                  //!< Close stream after finish writing
        void                                    openStream(bool reopen);                                            //!< Open the stream for writing
        void                                    printHeader(void);                                                  //!< Print header
        std::vector<unsigned int>               getChildCharacterCounts(int idx);
        std::vector<unsigned int>               getParentCharacterCounts(int idx);
        long                                    getNumSamples(void);
        
    private:
        std::string                             buildExtendedNewick();
        std::string                             buildExtendedNewick(TopologyNode* n);
        std::string                             buildCharacterHistoryString(TopologyNode* n, std::string brEnd="child");
        void                                    updateCharacterCounts(TopologyNode* n, std::string brEnd="child");
        std::string                             buildNhxString(void);
        
        // the stream to print
        std::fstream                            outStream;
        
        // parameters
        StochasticNode<AbstractHomologousDiscreteCharacterData>*  variable;
        TypedDagNode<Tree>*                     tree;
        const TimeAtlas*                        timeAtlas;
        std::set<DagNode *>                     nodeVariables;
        
        std::vector<std::vector<unsigned int> > parentCharacterCounts;
        std::vector<std::vector<unsigned int> > childCharacterCounts;
        
        std::vector<GeographicArea*>        areas;
        
        size_t numHistories;
        size_t numCharacters;
        
        path                                filename;
        std::string                         separator;
        bool                                posterior;
        bool                                prior;
        bool                                likelihood;
        bool                                append;
        bool                                showMetadata;
        bool                                showRates;
        long                                numSamples;
        std::uint64_t                       maxGen;
        long                                burn;
        
    };
    
}

#include "TreeHistoryCtmc.h"
#include "Model.h"
#include "RbException.h"

#include <sstream>

/* Constructor */
template<class charType>
RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::TreeCharacterHistoryNhxMonitor(StochasticNode<AbstractHomologousDiscreteCharacterData>* s, TypedDagNode<Tree>* t, const TimeAtlas* ta, std::uint64_t g, std::uint64_t mg, int b, const RevBayesCore::path &fname, const std::string &del, bool pp, bool l, bool pr, bool ap, bool sm, bool sr) :
Monitor(g,t),
outStream(),
variable( s ),
tree( t ),
timeAtlas( ta ),
filename( fname ),
separator( del ),
posterior( pp ),
prior( pr ),
likelihood( l ),
append(ap),
showMetadata(sm),
showRates(sr),
numSamples(0),
maxGen(mg),
burn(b) {
    
    nodes.push_back(s);
    s->incrementReferenceCount();
//    nodes.push_back(t);
    
    numHistories = tree->getValue().getNumberOfNodes();
    numCharacters = variable->getValue().getNumberOfCharacters();
    
    parentCharacterCounts.resize(numHistories);
    childCharacterCounts.resize(numHistories);
    for (size_t i = 0; i < numHistories; i++)
    {
        parentCharacterCounts[i].resize(numCharacters,0);
        childCharacterCounts[i].resize(numCharacters,0);
    }
    areas = timeAtlas->getAreas().back();
}

template<class charType>
RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::TreeCharacterHistoryNhxMonitor(const TreeCharacterHistoryNhxMonitor &m) :
Monitor( m ),
outStream( ),
variable( m.variable ),
tree( m.tree ),
timeAtlas( m.timeAtlas ),
nodeVariables( m.nodeVariables ),
parentCharacterCounts(m.parentCharacterCounts),
childCharacterCounts(m.childCharacterCounts),
areas(m.areas),
numHistories(m.numHistories),
numCharacters(m.numCharacters),
showMetadata(m.showMetadata),
showRates(m.showRates),
numSamples(m.numSamples),
maxGen(m.maxGen),
burn(m.burn) {
    
    filename    = m.filename;
    separator   = m.separator;
    prior       = m.prior;
    posterior   = m.posterior;
    likelihood  = m.likelihood;
    append      = m.append;
}


/* Clone the object */
template<class charType>
RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>* RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::clone(void) const {
    
    return new TreeCharacterHistoryNhxMonitor<charType>(*this);
}

template<class charType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::closeStream() {
    outStream.close();
}

template<class charType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::buildExtendedNewick( void ) {
    //tree->getValue().getRoot().setNewickNeedsRefreshing(true);
    numSamples++;
    std::string newick = buildExtendedNewick( &tree->getValue().getRoot() );
    return newick;
}

template<class charType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::updateCharacterCounts(TopologyNode* n, std::string brEnd)
{
    
    TreeHistoryCtmc<charType>* p = static_cast< TreeHistoryCtmc<charType>* >(&variable->getDistribution());
    const BranchHistory& bh = p->getHistory(*n);

    std::vector<CharacterEvent*> characters;
    if (brEnd=="child")
    {
        characters = bh.getChildCharacters();
        for (size_t i = 0; i < characters.size(); i++)
            childCharacterCounts[n->getIndex()][i] += static_cast<CharacterEventDiscrete*>(characters[i])->getState();
    }
    else if (brEnd=="parent")
    {
        characters = bh.getParentCharacters();
        for (size_t i = 0; i < characters.size(); i++)
            parentCharacterCounts[n->getIndex()][i] += static_cast<CharacterEventDiscrete*>(characters[i])->getState();
    }
    
}

template<class charType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::buildCharacterHistoryString(TopologyNode* n, std::string brEnd)
{
    size_t nd_idx = n->getIndex();
    TreeHistoryCtmc<charType>* p = static_cast< TreeHistoryCtmc<charType>* >(&variable->getDistribution());
    const BranchHistory& bh = p->getHistory(*n);

    std::vector<CharacterEvent*> characters;
    std::stringstream ss;
    
    if (brEnd=="child")
    {
        characters = bh.getChildCharacters();
        for (size_t i = 0; i < characters.size(); i++)
        {
            if (i != 0)
                ss << ",";
            ss << (double)childCharacterCounts[nd_idx][i] / numSamples;
        }
    }
    else if (brEnd=="parent")
    {
        characters = bh.getParentCharacters();
        for (size_t i = 0; i < characters.size(); i++)
        {
            if (i != 0)
                ss << ",";
            ss << (double)parentCharacterCounts[nd_idx][i] / numSamples;
        }
    }
    
    return ss.str();
}

/* Build newick string */
template<class charType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::buildExtendedNewick( TopologyNode* n ) {
    // create the newick string
    std::stringstream o;
    
    // extended data is only found on admixture nodes
    std::string additionalInfo = "";
    
    // loop over admixture nodes per branch
    std::stringstream characterStream;
    
    double br = 1.0;
    
    if (showMetadata)
    {
        updateCharacterCounts(n,"child");
        updateCharacterCounts(n,"parent");
        
        characterStream << "[";
        
        // character history
        characterStream << "&ch={" << buildCharacterHistoryString(n,"child") << "}";
        characterStream << ",&pa={" << buildCharacterHistoryString(n,"parent") << "}";
        
        // ... whatever else
        characterStream << "]";
        
        additionalInfo = characterStream.str();
    }
    
    // test whether this is a internal or external node
    if (n->isTip()) {
        // this is a tip so we just return the name of the node
        o << n->getIndex() << additionalInfo << ":" << n->getBranchLength();
    }
    
    else {
        o << "(";
        for (size_t i=0; i<(n->getNumberOfChildren()-1); i++) {
            o << buildExtendedNewick( &n->getChild(i) ) << ",";
        }
        o << buildExtendedNewick( &n->getChild(n->getNumberOfChildren()-1) ) << ")" << additionalInfo << ":" << n->getBranchLength() * br;
    }
    
    return o.str();
}


/** Monitor value at generation gen */
template<class charType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::monitor(std::uint64_t gen) {
    
    // get the printing frequency
    std::uint64_t samplingFrequency = printgen;
    
    if (gen % samplingFrequency == 0 && gen != maxGen) { // && gen >= burn) {
        
        std::string en = buildExtendedNewick();
        //std::cout << "\n" << en << "\n\n";
        
    }
    
    if (gen == maxGen)
    {
        std::string nhxStr = buildNhxString();
        outStream << nhxStr << std::endl;
    }
}

template<class charType>
std::string RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::buildNhxString(void)
{
    std::stringstream nhxStrm;
    
    // begin nexus file
    nhxStrm << "#NEXUS" << "\n\n";
    
    // Nhx settings block
    nhxStrm << "Begin phylowood;\n";
    nhxStrm << "\tdrawtype pie\n";
    nhxStrm << "\tmodeltype biogeography\n";
    nhxStrm << "\tareatype discrete\n";
    nhxStrm << "\tmaptype clean\n";
    nhxStrm << "\tpieslicestyle even\n";
    nhxStrm << "\tpiefillstyle outwards\n";
    nhxStrm << "\ttimestart -" << tree->getValue().getRoot().getAge() << "\n";
    nhxStrm << "\tmarkerradius " << 200 << "\n";
    nhxStrm << "\tminareaval " << 0.1 << "\n";
    nhxStrm << "End;\n\n";
    
    // bayarea-fig block
    nhxStrm << "Begin bayarea-fig;\n";
    nhxStrm << "\tmapheight\t100\n";
    nhxStrm << "\tmapwidth\t150\n";
    nhxStrm << "\tcanvasheight\t2000\n";
    nhxStrm << "\tcanvaswidth\t1000\n";
    nhxStrm << "\tminareaval\t0.15\n";
    nhxStrm << "\tareacolors black\n";
    nhxStrm << "\tareatypes";
    for (unsigned i = 0; i < numCharacters; i++)
        nhxStrm << " 1";
    nhxStrm << "\n";
    nhxStrm << "\tareanames Default\n";
    nhxStrm << "End;\n\n";
    
    // taxa block
    nhxStrm << "Begin taxa;\n";
    nhxStrm << "\tDimensions ntax=" << tree->getValue().getNumberOfTips() << ";\n";
    nhxStrm << "\tTaxlabels\n";
    for (unsigned i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
        TopologyNode* p = &tree->getValue().getNode(i);
        if (p->isTip())
        {
            nhxStrm << "\t\t" << p->getName() << "\n";
        }
    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";
    
    // geo block
    nhxStrm << "Begin geo;\n";
    nhxStrm << "\tDimensions ngeo=" << numCharacters << ";\n";
    nhxStrm << "\tCoords\n";
    
    for (unsigned i = 0; i < numCharacters; i++)
    {
//        nhxStrm << "\t\t" << i << "\t" << geographicCoordinates[i][0] << "\t" << geographicCoordinates[i][1];
        nhxStrm << "\t\t" << i << "\t" << areas[i]->getLatitude() << "\t" << areas[i]->getLongitude();
        if (i < (numCharacters - 1))
            nhxStrm << ",";
        nhxStrm << "\n";
    }
    nhxStrm << "\t\t;\n";
    nhxStrm << "End;\n\n";
    
    // tree block
    nhxStrm << "Begin trees;\n";
    nhxStrm << "\tTranslate\n";
    for (unsigned i = 0; i < tree->getValue().getNumberOfNodes(); i++)
    {
        TopologyNode* p = &tree->getValue().getNode(i);
        if (p->isTip())
        {
            nhxStrm << "\t\t" << p->getIndex() << "\t" << p->getName();
            if (i < (tree->getValue().getNumberOfNodes() - 1))
                nhxStrm << ",";
            nhxStrm << "\n";
        }
    }
    nhxStrm << "\t\t;\n";
    
    // write tree string
    std::string treeStr = "";
    treeStr = buildExtendedNewick();
    nhxStrm << "tree TREE1 = " << treeStr << ";\n";
    nhxStrm << "End;\n";
    
    return nhxStrm.str();
}


/** open the file stream for printing */
template<class charType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::openStream(bool reopen)
{
    createDirectoryForFile( filename );
    
    // open the stream to the file
    if (append)
        outStream.open( filename.string(), std::fstream::out | std::fstream::app);
    else
        outStream.open( filename.string(), std::fstream::out);
}

/** Print header for monitored values */
template<class charType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::printHeader()
{
    
    // do nothing
    ;
}

template<class charType>
std::vector<unsigned int> RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::getChildCharacterCounts(int idx)
{
    return childCharacterCounts[idx];
}

template<class charType>
std::vector<unsigned int> RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::getParentCharacterCounts(int idx)
{
    return parentCharacterCounts[idx];
}

template<class charType>
long RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::getNumSamples(void)
{
    return numSamples;
}

template<class charType>
void RevBayesCore::TreeCharacterHistoryNhxMonitor<charType>::swapNode(DagNode *oldN, DagNode *newN) {
    
    bool found = false;
    if ( oldN == tree )
    {
        tree = static_cast< TypedDagNode<Tree> *>(newN);
        found = true;
    }
    else if ( oldN == variable )
    {
        variable = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN);
        found = true;
    }
    
//    for (size_t i = 0; i < branchHistories.size(); i++)
//    {
//        if (oldN == branchHistories[i])
//        {
//            branchHistories[i] = static_cast<StochasticNode<BranchHistory>* >(newN);
//            found = true;
//        }
//    }
//    
//    /*
//     if (found == false)
//     {
//     // error catching
//     if ( nodeVariables.find(oldN) == nodeVariables.end() ) {
//     throw RbException("Cannot replace DAG node in this monitor because the monitor doesn't hold this DAG node.");
//     }
//     
//     nodeVariables.erase( oldN );
//     nodeVariables.insert( newN );
//     }
//     */
    
    // delegate to base class
    Monitor::swapNode(oldN, newN);
}


#endif
