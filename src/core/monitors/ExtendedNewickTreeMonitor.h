/**
 * @file
 * This file contains the declaration of a ExtendedNewickTreeMonitor, used to save information
 * to a file about the tree and additional variable on nodes (or branches) using
 * the extended Newick formay.
 *
 * @brief Declaration of ExtendedNewickTreeMonitor
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-11-26 05:34:51 -0800 (Mon, 26 Nov 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-12-11, version 1.0
 *
 * $Id: ExtendedNewickTreeMonitor.h 1867 2012-11-26 13:34:51Z hoehna $
 */

#ifndef ExtendedNewickTreeMonitor_H
#define ExtendedNewickTreeMonitor_H

#include <fstream>
#include <cstdint>
#include <vector>

#include "VariableMonitor.h"

namespace RevBayesCore {
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;
    
    class ExtendedNewickTreeMonitor : public VariableMonitor {
        
    public:
        // Constructors and Destructors
        ExtendedNewickTreeMonitor(TypedDagNode<Tree> *t, const std::vector<DagNode*> &n, bool np, std::uint64_t g, const std::string &fname,
                                  const std::string &del, bool pp=true, bool l=true, bool pr=true, bool ap=false);                                              //!< Constructor with set of DAG node
        
        // basic methods
        ExtendedNewickTreeMonitor*          clone(void) const;                                                      //!< Clone the object
        
        // Monitor functions
        void                                monitorVariables(std::uint64_t gen);                                    //!< Monitor at generation gen
        void                                swapNode(DagNode *oldN, DagNode *newN);

        // FileMonitor functions
        void                                printFileHeader(void);                                                  //!< Print header
        
    private:        
        // parameters
        bool                                isNodeParameter;
        TypedDagNode<Tree>*                 tree;
        std::vector<DagNode*>               nodeVariables;        
    };
    
}

#endif

