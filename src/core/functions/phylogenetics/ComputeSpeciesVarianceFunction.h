#ifndef ComputeSpeciesVarianceFunction_H
#define ComputeSpeciesVarianceFunction_H

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

    class ComputeSpeciesVarianceFunction : public TypedFunction< RbVector<double> > {

    public:
        enum                                                                MISSING_TREATMENT { MEAN, MEDIAN, NONE };
        ComputeSpeciesVarianceFunction(const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t>* s, const std::vector<Taxon> &ta, MISSING_TREATMENT mtr, bool err );
        virtual                                                ~ComputeSpeciesVarianceFunction(void);                                                         //!< Virtual destructor

        // public member functions
        ComputeSpeciesVarianceFunction*                              clone(void) const;                                                                  //!< Create an independent clone
        void                                                                update(void);

    protected:
        double                                                              computeMeanForSpecies(const std::string &n, size_t i);
        double                                                              computeTipErrorOrVarianceForSpecies(const std::string &n, size_t i);
        double                                                              computeMeanErrorOrVarianceAcrossSpecies(void);
        double                                                              computeMedianErrorOrVarianceAcrossSpecies(void);
        double                                                              getNumberOfSamplesForSpecies(const std::string &n);
        std::vector<std::string>                                            getAlphabeticalSpeciesNames(void);
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                                reset(void);

    private:

        // members
        const TypedDagNode<ContinuousCharacterData>*                        data;
        const TypedDagNode<std::int64_t>*                                   site;
        std::vector<Taxon>                                                  taxa;
        MISSING_TREATMENT                                                   missing_var_treatment;

        std::vector<double>                                                 within_species_variance;
        const bool                                                          compute_VarOfMean;
    };

}

#endif
