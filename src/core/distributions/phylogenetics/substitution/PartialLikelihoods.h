#ifndef PartialLikelihoods_H
#define PartialLikelihoods_H

#include <vector>

constexpr double scale_factor = 115792089237316195423570985008687907853269984665640564039457584007913129639936e0;
constexpr double scale_min = 1.0/scale_factor;
constexpr double log_scale_factor = 177.445678223345999210811423093293201427328034396225345054e0;

template <typename T, std::size_t Alignment = alignof(T)>
struct no_init_allocator {
    using value_type = T;

    static constexpr std::align_val_t alignment{
        std::max(Alignment, alignof(T))
    };

    no_init_allocator() = default;

    template <class U>
    no_init_allocator(const no_init_allocator<U, Alignment>&) noexcept {}

    template <class U>
    bool operator==(const no_init_allocator<U, Alignment>&) const noexcept { return true; }

    T* allocate(std::size_t n) {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_array_new_length();
        return static_cast<T*>(::operator new(n * sizeof(T), alignment));
    }

    void deallocate(T* p, std::size_t /*n*/) noexcept {
        ::operator delete(p, alignment);
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

    // Required for allocator rebinding
    template <class U>
    struct rebind {
        using other = no_init_allocator<U, Alignment>;
    };
};

template <typename T, std::size_t Alignment = alignof(T)>
using no_init_vector = std::vector<T, no_init_allocator<T, Alignment>>;

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
    double& likelihood(int m, int p, int s)
    {
        assert(0 <= m and m < dims_.num_site_mixtures);
        assert(0 <= p and p < dims_.num_patterns);
        assert(0 <= s and s < dims_.num_states);
        return likelihoods[s + dims_.num_states*(p + dims_.num_patterns*m)];
    }
    double  likelihood(int m, int p, int s) const
    {
        assert(0 <= m and m < dims_.num_site_mixtures);
        assert(0 <= p and p < dims_.num_patterns);
        assert(0 <= s and s < dims_.num_states);
        return likelihoods[s + dims_.num_states*(p + dims_.num_patterns*m)];
    }

    no_init_vector<double,32> likelihoods; // per mixture * pattern * state
    no_init_vector<int> scale; // per site

    PartialLikelihoods& operator=(const PartialLikelihoods&) = delete;
    PartialLikelihoods& operator=(PartialLikelihoods&&) noexcept = default;

    PartialLikelihoods(const PartialLikelihoods&) = delete;
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
