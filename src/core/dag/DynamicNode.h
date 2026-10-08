#ifndef DynamicNode_H
#define DynamicNode_H

#include "DynamicNodeBase.h"
#include "TypedDagNode.h"

namespace RevBayesCore {

    template<class valueType>
    class DynamicNode : public TypedDagNode<valueType>, protected DynamicNodeBase {

    public:
        DynamicNode(const std::string &n);
        DynamicNode(const DynamicNode &n);
        virtual                                            ~DynamicNode(void) = default;                                                //!< Virtual destructor

        // pure virtual methods
        virtual DynamicNode<valueType>*                    clone(void) const = 0;

        // public methods
        virtual DagNode*                                   cloneDAG(DagNodeMap &nodesMap, std::map<std::string, const DagNode*> &names) const; //!< Clone the entire DAG which is connected to this node

        // this function provided for derived classes used in the language layer, which need to override it
        virtual const std::string&                         getRevTypeOfValue(void);                                                     //!< Get Rev language type of value

    protected:
        virtual void                                       keepMe(const DagNode* affecter);                                             //!< Keep value of this and affected nodes
        virtual void                                       restoreMe(const DagNode *restorer);                                          //!< Restore value of this node
        virtual void                                       touchMe(const DagNode *toucher, bool touchAll);                              //!< Tell affected nodes value is reset
    };

}


/** Initialize the typed and type-independent portions of a dynamic node. */
template<class valueType>
RevBayesCore::DynamicNode<valueType>::DynamicNode(const std::string &n) : TypedDagNode<valueType>( n ),
    DynamicNodeBase()
{
    // nothing to do here
}


/** Copy the typed node while starting the dynamic state as touched. */
template<class valueType>
RevBayesCore::DynamicNode<valueType>::DynamicNode(const DynamicNode<valueType> &n) : TypedDagNode<valueType>( n ),
    DynamicNodeBase( n )
{
    // nothing to do here
}


/** Delegate type-independent graph cloning to the non-template implementation base. */
template<class valueType>
RevBayesCore::DagNode* RevBayesCore::DynamicNode<valueType>::cloneDAG(DagNodeMap &nodesMap, std::map<std::string, const DagNode*> &names) const
{
    return DynamicNodeBase::cloneDAG( *this, nodesMap, names );
}


/** Return the default type-erased error for core dynamic nodes. */
template<class valueType>
const std::string& RevBayesCore::DynamicNode<valueType>::getRevTypeOfValue(void)
{
    return DynamicNodeBase::getRevTypeOfValue();
}


/** Forward keep-state bookkeeping to the implementation base. */
template<class valueType>
void RevBayesCore::DynamicNode<valueType>::keepMe(const DagNode* /*affecter*/)
{
    DynamicNodeBase::keepMe();
}


/** Forward restore-state bookkeeping to the implementation base. */
template<class valueType>
void RevBayesCore::DynamicNode<valueType>::restoreMe(const DagNode* /*restorer*/)
{
    DynamicNodeBase::restoreMe();
}


/** Forward touch-state bookkeeping to the implementation base. */
template<class valueType>
void RevBayesCore::DynamicNode<valueType>::touchMe(const DagNode* /*toucher*/, bool /*touchAll*/)
{
    DynamicNodeBase::touchMe();
}

#endif
