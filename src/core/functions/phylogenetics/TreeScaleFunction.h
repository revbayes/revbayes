#ifndef TreeScaleFunction_H
#define TreeScaleFunction_H

#include <vector>

#include "Tree.h"
#include "TypedFunction.h"

namespace RevBayesCore {
class DagNode;
template <class valueType> class TypedDagNode;

/**
     * @brief The Tree Scale Function
     *
     * A core function for scaling the length of the tree. This function scales internal edges of the tree and adjusts the tip ages to the age specified in the parameter TipAges
     * The class also has the following member variables:
     *  @param tau the phylogenetic tree to be scaled
     *  @param scale the value to multiply internal node ages by
     *  @param tipAges the new ages of the tips after scaling
     *  @param scaleLimit The scale parameter should always be less than the max of the tipAges, denoted here
     */
    class TreeScaleFunction : public TypedFunction<Tree> {
        
    public:
        TreeScaleFunction(const TypedDagNode<Tree> *t, const TypedDagNode<double> *b, std::vector<double> m);
        virtual                                            ~TreeScaleFunction(void);                                                         //!< Virtual destructor
        
        // public member functions
        TreeScaleFunction*                                  clone(void) const;                                                                  //!< Create an independent clone
        void                                                reInitialized(void);                                                                //!< The arguments have been re-initialized
        void                                                update(void);
        
    protected:
        void                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode<Tree>*                           tau;
        const TypedDagNode<double>*                         scale;
        std::vector<double>                                 tipAges;
        double                                              scaleLimit;
    };
    
}

#endif
