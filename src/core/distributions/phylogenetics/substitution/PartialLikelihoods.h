#ifndef PartialLikelihoods_H
#define PartialLikelihoods_H

#include <vector>

constexpr double scale_factor = 115792089237316195423570985008687907853269984665640564039457584007913129639936e0;
constexpr double scale_min = 1.0/scale_factor;
constexpr double log_scale_factor = 177.445678223345999210811423093293201427328034396225345054e0;

template <typename T>
struct no_init_allocator {
    using value_type = T;

    no_init_allocator() = default;

    template <class U>
    no_init_allocator(const no_init_allocator<U>&) noexcept {}

    bool operator==(const no_init_allocator<T>&) const = default;

    T* allocate(std::size_t n) {
        return std::allocator<T>{}.allocate(n);
    }
    void deallocate(T* p, std::size_t n) {
        std::allocator<T>{}.deallocate(p, n);
    }

    // Called with no args during value-init → skip zeroing
    template <typename U>
    void construct(U* p) noexcept {
        ::new (static_cast<void*>(p)) U; // default-init, NOT value-init
    }

    // All other construct calls (copy, move, etc.) behave normally
    template <typename U, typename... Args>
        void construct(U* p, Args&&... args) {
        ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
    }
};

template <typename T>
using no_init_vector = std::vector<T, no_init_allocator<T>>;

class PartialLikelihoods
{
public:
    // With MPI, these might not be ALL the site patterns, but they are all the ones stored in this process.
    struct Dims
    {
        std::size_t num_site_mixtures;
        std::size_t num_patterns;
        std::size_t num_states;

        size_t size() const {return num_site_mixtures * num_patterns * num_states;}

        bool operator==(const Dims&) const = default;

    Dims(std::size_t nm, std::size_t np, std::size_t ns)
        :num_site_mixtures(nm),
            num_patterns(np),
            num_states(ns)
            { }
    };

private:
    // We can't change this without reallocating.
    Dims dims_;

public:

    const Dims& dims() const {return dims_;}

    // This shows how the entries are laid out inside the linear array.
    double& likelihood(int m, int p, int s)       {return likelihoods[s + dims_.num_states*(p + dims_.num_patterns*m)];}
    double  likelihood(int m, int p, int s) const {return likelihoods[s + dims_.num_states*(p + dims_.num_patterns*m)];}

    no_init_vector<double> likelihoods; // per mixture * pattern * state
    no_init_vector<int> scale; // per site

    PartialLikelihoods& operator=(const PartialLikelihoods&) = default;
    PartialLikelihoods& operator=(PartialLikelihoods&&) noexcept = default;

    PartialLikelihoods(const PartialLikelihoods&) = default;
    PartialLikelihoods(PartialLikelihoods&&) noexcept = default;

PartialLikelihoods()
    :dims_{0,0,0}
    {}

PartialLikelihoods(const Dims& d)
    :dims_(d),
     likelihoods(dims_.size()),
     scale(dims_.num_patterns)
    {
    }
};

#endif
