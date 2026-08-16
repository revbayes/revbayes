#ifndef DagNodeTypeUtilities_H
#define DagNodeTypeUtilities_H

#include "DagNode.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TypedDagNode.h"

namespace RevBayesCore {

    /** Return the node with the requested value type, or null if its value has another type. */
    template <typename valueType>
    TypedDagNode<valueType>* checkNodeOfType(DagNode *node)
    {
        return dynamic_cast<TypedDagNode<valueType>*>( node );
    }

    /** Return the node with the requested value type, or null if its value has another type. */
    template <typename valueType>
    const TypedDagNode<valueType>* checkNodeOfType(const DagNode *node)
    {
        return dynamic_cast<const TypedDagNode<valueType>*>( node );
    }

    /** Recover the node's known value type, checking the invariant in debug builds. */
    template <typename valueType>
    TypedDagNode<valueType>* assumeNodeHolds(DagNode *node)
    {
#ifndef NDEBUG
        TypedDagNode<valueType> *typedNode = checkNodeOfType<valueType>( node );
        if ( typedNode == NULL )
        {
            throw RbException("A DAG node does not hold the expected value type.");
        }
        return typedNode;
#else
        return static_cast<TypedDagNode<valueType>*>( node );
#endif
    }

    /** Recover the node's known value type, checking the invariant in debug builds. */
    template <typename valueType>
    const TypedDagNode<valueType>* assumeNodeHolds(const DagNode *node)
    {
#ifndef NDEBUG
        const TypedDagNode<valueType> *typedNode = checkNodeOfType<valueType>( node );
        if ( typedNode == NULL )
        {
            throw RbException("A DAG node does not hold the expected value type.");
        }
        return typedNode;
#else
        return static_cast<const TypedDagNode<valueType>*>( node );
#endif
    }

    /** Return the value of a node whose value type is established by the caller. */
    template <typename valueType>
    valueType& getValueOfType(DagNode *node)
    {
        return assumeNodeHolds<valueType>( node )->getValue();
    }

    /** Return the value of a node whose value type is established by the caller. */
    template <typename valueType>
    const valueType& getValueOfType(const DagNode *node)
    {
        return assumeNodeHolds<valueType>( node )->getValue();
    }

    /** Return the requested stochastic node, or null if the node has another type or behavior. */
    template <typename valueType>
    StochasticNode<valueType>* checkStochastic(DagNode *node)
    {
        return dynamic_cast<StochasticNode<valueType>*>( node );
    }

    /** Return the requested stochastic node, or null if the node has another type or behavior. */
    template <typename valueType>
    const StochasticNode<valueType>* checkStochastic(const DagNode *node)
    {
        return dynamic_cast<const StochasticNode<valueType>*>( node );
    }

    /** Require a stochastic node of the requested value type, or report an invalid model. */
    template <typename valueType>
    StochasticNode<valueType>* requireStochastic(DagNode *node)
    {
        StochasticNode<valueType> *stochasticNode = checkStochastic<valueType>( node );
        if ( stochasticNode == NULL )
        {
            throw RbException("A DAG node is not stochastic with the expected value type.");
        }
        return stochasticNode;
    }

    /** Require a stochastic node of the requested value type, or report an invalid model. */
    template <typename valueType>
    const StochasticNode<valueType>* requireStochastic(const DagNode *node)
    {
        const StochasticNode<valueType> *stochasticNode = checkStochastic<valueType>( node );
        if ( stochasticNode == NULL )
        {
            throw RbException("A DAG node is not stochastic with the expected value type.");
        }
        return stochasticNode;
    }

    /** Recover a known stochastic-node type, checking the invariant in debug builds. */
    template <typename valueType>
    StochasticNode<valueType>* assumeStochastic(DagNode *node)
    {
#ifndef NDEBUG
        return requireStochastic<valueType>( node );
#else
        return static_cast<StochasticNode<valueType>*>( node );
#endif
    }

    /** Recover a known stochastic-node type, checking the invariant in debug builds. */
    template <typename valueType>
    const StochasticNode<valueType>* assumeStochastic(const DagNode *node)
    {
#ifndef NDEBUG
        return requireStochastic<valueType>( node );
#else
        return static_cast<const StochasticNode<valueType>*>( node );
#endif
    }

    /** Replace a stored node pointer with a compatible node, checking the invariant in debug builds. */
    template <typename TargetNode, typename SourceNode>
    void replaceNodeReference(TargetNode *&node, SourceNode *replacement)
    {
#ifndef NDEBUG
        TargetNode *typedReplacement = dynamic_cast<TargetNode*>( replacement );
        if ( typedReplacement == NULL )
        {
            throw RbException("Cannot replace a DAG-node reference with an incompatible node.");
        }
        node = typedReplacement;
#else
        node = static_cast<TargetNode*>( replacement );
#endif
    }

}

#endif
