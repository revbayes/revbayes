#ifndef ComputeEmpiricalWithinSpeciesVariancesFunction_H
#define ComputeEmpiricalWithinSpeciesVariancesFunction_H

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

    class ComputeEmpiricalWithinSpeciesVariancesFunction : public TypedFunction< RbVector<double> > {

    public:
        ComputeEmpiricalWithinSpeciesVariancesFunction(const TypedDagNode<Tree> *t, const TypedDagNode<ContinuousCharacterData> *d, const TypedDagNode<std::int64_t>* s, const std::vector<Taxon> &ta, bool lt, RevBayesCore::TypedDagNode<double>* dv );
        virtual                                                ~ComputeEmpiricalWithinSpeciesVariancesFunction(void);                                                         //!< Virtual destructor

        // public member functions
        ComputeEmpiricalWithinSpeciesVariancesFunction*                           clone(void) const;                                                                  //!< Create an independent clone
        void                                                                      update(void);

    protected:
        double                                                              computeMeanForSpecies(const std::string &n, size_t i);
        double                                                              computeWithinSpeciesVariance(const std::string &n, size_t i, bool lt);
        double                                                              getNumberOfSamplesForSpecies(const std::string &n);
        void                                                                swapParameterInternal(const DagNode *oldP, const DagNode *newP);                    //!< Implementation of swaping parameters
        void                                                                resetWithinSpeciesVariances(void);

    private:

        // members
        const TypedDagNode<Tree>*                                           tau;
        const TypedDagNode<ContinuousCharacterData>*                        data;
        const TypedDagNode<std::int64_t>*                                   site;
        bool                                                                log_transformed;
        std::vector<Taxon>                                                  taxa;
        RevBayesCore::TypedDagNode<double>*                                 default_var;

        std::vector<double>                                                 within_species_variance;

    };

}

#endif
