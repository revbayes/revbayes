#ifndef SITE_MIXTURE_MODEL_H
#define SITE_MIXTURE_MODEL_H

#include <stddef.h>
#include <vector>
#include <tuple>
#include <optional>

//#include "MatrixReal.h"
#include "MemberObject.h"  // For member functions.
#include "Cloneable.h"
#include "TransitionProbabilityMatrix.h"
#include "SiteModel.h"

namespace RevBayesCore {

    class Tree;

    /**
     * @brief SiteModel mixture class.
     *
     * Derived models that depend on the branch lengths will additionally need to
     * incorporate a reference to the tree.
     *
     * @copyright Copyright 2022-
     * @author Benjamin Redelings
     * @since 2022-09-22, version 1.0
     */

    class SiteMixtureModel: public Cloneable,
			    public MemberObject< RbVector<RbVector<RbVector<double>>> >,
			    public MemberObject< double >
    {
        std::vector<std::shared_ptr<const SiteModel>> components;
        std::vector<double> fractions;

    public:
        virtual                          ~SiteMixtureModel() = default;                                                  //!< Destructor
        SiteMixtureModel*                clone() const;

	const SiteModel&                 getComponent(int m) const;
        const std::vector<double>&       componentProbs() const;
	int                              size() const;
        int                              getNumberOfComponents() const;

        int                              getNumberOfStates() const;

        std::vector<TransitionProbabilityMatrix>
                                         calculateTransitionProbabilities(const Tree& t, int node, double rate) const;

	std::optional<double>            rate() const;
        void                             scale(double f);
        void                             setRate(double r);

        virtual void                     executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<RbVector<RbVector<double>>> &retValue) const;       //!< Execute the member-method
        virtual void                     executeMethod( const std::string &n, const std::vector<const DagNode*> &args, double &retValue) const;       //!< Execute the member-method

        // This is a hack to satisfy ModelVector<T>, which incorrectly assumes that these exist for all T.
        bool                             operator==(const SiteMixtureModel&) const {return false;}
        bool                             operator!=(const SiteMixtureModel&) const {return true;}
        bool                             operator<=(const SiteMixtureModel&) const {return false;}
        bool                             operator<(const SiteMixtureModel&) const {return false;}

        // virtual std::vector<int>            get_emitted_letters() const;                                          //!<Find out what alphet letter each state emits, for markov modulated models.

	SiteMixtureModel& operator=(const SiteMixtureModel&) = default;
	SiteMixtureModel& operator=(SiteMixtureModel&&) = default;

	SiteMixtureModel(const SiteMixtureModel&) = default;
	SiteMixtureModel(SiteMixtureModel&&) noexcept = default;

	SiteMixtureModel() = default; // This is required by RbVectorImpl::initFromString because the class is not abstract.
	SiteMixtureModel(const std::vector<std::shared_ptr<const SiteModel>>& c, const std::vector<double>& f);
	SiteMixtureModel(std::vector<std::shared_ptr<const SiteModel>>&& c, std::vector<double>&& f);
    };

    // We need this for TypedDagNode<SiteMixtureModel> for some reason...
    std::ostream&                                       operator<<(std::ostream& o, const SiteMixtureModel& x);

    template <typename T>    
    std::vector<std::shared_ptr<const T>> scale_models(const std::vector<std::shared_ptr<const T>>& models, const std::vector<double>& rates)
    {
	assert(models.size() == rates.size());

	std::vector<std::shared_ptr<const T>> scaled_models;

	for(int i=0; i<models.size(); i++)
	{
	    auto m = std::shared_ptr<T>(models[i]->clone());
	    m->scale(rates[i]);
	    scaled_models.push_back(m);
	}

	return scaled_models;
    }

    std::shared_ptr<const SiteMixtureModel> scaled_mixture(const SiteMixtureModel& sub_model, const std::vector<double>& fractions, const std::vector<double>& rates);
    std::shared_ptr<const SiteMixtureModel> mix_mixture(const std::vector<std::shared_ptr<const SiteMixtureModel>>& submodel, const std::vector<double>& fractions);
}



#endif
