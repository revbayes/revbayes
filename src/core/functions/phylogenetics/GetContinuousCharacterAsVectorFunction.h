#ifndef GetContinuousCharacterAsVectorFunction_H
#define GetContinuousCharacterAsVectorFunction_H

#include <cstddef>
#include <vector>
#include <iosfwd>

#include "RbVector.h"
#include "TypedFunction.h"
#include "Taxon.h"
#include "TopologyNode.h"

namespace RevBayesCore {
class ContinuousCharacterData;
class DagNode;
class Tree;
template <class valueType> class TypedDagNode;

    class GetContinuousCharacterAsVectorFunction : public TypedFunction< RbVector<double> > {

    public:
        enum                                                                VECTOR_ORDER { ALPHABETICAL };
        GetContinuousCharacterAsVectorFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t>* s, VECTOR_ORDER ord );
        virtual                                                             ~GetContinuousCharacterAsVectorFunction(void);                                                         //!< Virtual destructor

        // public member functions
        GetContinuousCharacterAsVectorFunction*                             clone(void) const;                                                                  //!< Create an independent clone
        void                                                                update(void);

    protected:
        double                                                              getContinuousCharacter(const std::string &n, size_t site_index);
        std::vector<std::string>                                            getAlphabeticalSpeciesNames(void);
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                                reset(void);

    private:

        // members
        const TypedDagNode<ContinuousCharacterData>*                        data;
        const TypedDagNode<std::int64_t>*                                   site_index;
        VECTOR_ORDER                                                        order_by;

        std::vector<double>                                                 continuous_character;

    };

}

#endif
