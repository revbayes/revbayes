#ifndef TajimasPiFunction_H
#define TajimasPiFunction_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "TypedDagNode.h"
#include "TypedFunction.h"

#include <vector>

namespace RevBayesCore {
    
    /**
     * @brief Computing Tajima's pi function.
     *
     * This function computes Tajima's pi, which is an estimator of theta = 4*N*mu.
     * This function can be used to compute the per site or per locus estimator.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since Version 1.0, 2015-04-30
     *
     */
    class TajimasPiFunction : public TypedFunction<double> {
        
    public:
        TajimasPiFunction(const TypedDagNode<AbstractHomologousDiscreteCharacterData> *a, bool ps, bool excl);
        virtual                                                        ~TajimasPiFunction(void);                                                    //!< Virtual destructor
        
        // public member functions
        TajimasPiFunction*                                              clone(void) const;                                                              //!< Create an independent clone
        void                                                            update(void);
        
    protected:
        void                                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                        //!< Implementation of swaping parameters
        
    private:
        
        // members
        const TypedDagNode< AbstractHomologousDiscreteCharacterData >*  alignment;
        bool                                                            per_site;
        bool                                                            exclude_ambiguous_sites;
    };
    
}

#endif
