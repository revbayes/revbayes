#ifndef Model_H
#define Model_H

#include <map>
#include <set>
#include <vector>

#include "Cloneable.h"
#include "DagNodeMap.h"
#include "Parallelizable.h"

namespace RevBayesCore {

    class DagNode;
    
    /**
     * \brief The Model class holds its independent copy of the model graph (DAG nodes) contained in the model
     * and provides methods for convenient access.
     *
     * A Model object holds the model graph which is an independent copy
     * of the DAG nodes contained in the model. The model graph may or may not be connected.
     * The model graph is obtained by providing a set of DAG nodes to the model and the model
     * object extracts all DAG nodes that are connected to these provided nodes.
     * An independent copy of these nodes is then obtained. 
     * The model object simply provides access to these DAG nodes.
     *
     * \copyright (c) Copyright 2009-2013 (GPL version 3)
     * \author The RevBayes Development Core Team (Sebastian Hoehna)
     * \since Version 1.0, 2012-06-21
     *
     */
    class Model : public Cloneable, public Parallelizable {
        
        public:
                                                                    Model(const DagNode* source);                                   //!< Constructor from a single DAG node from which the model graph is extracted.
                                                                    Model(const std::set<const DagNode*> &sources);                 //!< Constructor from a set of DAG nodes from each of which the model graph is extracted.
                                                                    Model(const Model &m);                                          //!< Copy constructor.
        virtual                                                    ~Model(void);                                                    //!< Destructor.
    
        
        // overloaded operators
        Model&                                                      operator=(const Model& x);                                      //!< Assignment operator.
        
        
        // convenience methods
        virtual Model*                                              clone(void) const;                                              //!< Clone function. This is similar to the copy constructor but useful in inheritance.
    
        // getters and setters
        std::vector<DagNode*>&                                      getDagNodes(void);                                              //!< Non-constant getter function of the set of DAG nodes contained in the model graph.
        const std::vector<DagNode*>&                                getDagNodes(void) const;                                        //!< Constant getter function of the set of DAG nodes contained in the model graph.
        const DagNodeMap&                                           getNodesMap(void) const;                                        //!< Constant getter function of the map between the pointer of the original DAG nodes to the pointers of the copied DAG nodes.
        std::vector<DagNode*>                                       getOrderedStochasticNodes(void);
        
    protected:
        void                                                        setActivePIDSpecialized(size_t i, size_t n);                    //!< Set the number of processes for this class.

        
    private:
        
        // private methods
        void                                                        addSourceNode(const DagNode *sourceNode);                       //!< Add a source node, extract the model graph and create and indepedent copy of it.
        void                                                        getOrderedStochasticNodes(  const DagNode*    dagNode,
                                                                                      std::vector<DagNode*>&      orderedStochasticNodes,
                                                                                      std::set<const DagNode*>&   visitedNodes);
        
        // members
        std::vector<DagNode*>                                       nodes;                                                          //!< The DAG nodes of the model graph. These need to be pointers because we don't actually know there specific type. We own these.
        DagNodeMap                                                  nodesMap;                                                       //!< Map between original nodes and own copy.
        std::vector<const DagNode*>                                 sources;                                                        //!< Set of source nodes for the model graph.
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const Model& x);                                //!< Overloaded output operator

    
}

#endif
