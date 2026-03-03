#ifndef ComputeWithinSpeciesVarianceFromCharacterDataFunction_H
#define ComputeWithinSpeciesVarianceFromCharacterDataFunction_H

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

    class ComputeWithinSpeciesVarianceFromCharacterDataFunction : public TypedFunction< RbVector<double> > {

    public:
        enum                                                                MISSING_TREATMENT { MEAN, MEDIAN, NONE };
        ComputeWithinSpeciesVarianceFromCharacterDataFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t>* vs, const TypedDagNode<std::int64_t>* ns, MISSING_TREATMENT mtr );
        virtual                                                ~ComputeWithinSpeciesVarianceFromCharacterDataFunction(void);                                                         //!< Virtual destructor

        // public member functions
        ComputeWithinSpeciesVarianceFromCharacterDataFunction*                           clone(void) const;                                                                  //!< Create an independent clone
        void                                                                      update(void);

    protected:
        double                                                              computeWithinSpeciesVariance(const std::string &n, size_t v_site, size_t n_site);
        double                                                              computeMeanWithinSpeciesVariance(void);
        double                                                              computeMedianWithinSpeciesVariance(void);
        double                                                              getNumberOfSamplesForSpecies(const std::string &n, size_t site);
        std::vector<std::string>                                            getAlphabeticalSpeciesNames(void);
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                                resetWithinSpeciesVariances(void);

    private:

        // members
        const TypedDagNode<ContinuousCharacterData>*                        data;
        const TypedDagNode<std::int64_t>*                                   variance_site;
        const TypedDagNode<std::int64_t>*                                   num_sample_site;
        MISSING_TREATMENT                                                   missing_var_treatment;

        std::vector<double>                                                 within_species_variance;

    };

}

#endif
