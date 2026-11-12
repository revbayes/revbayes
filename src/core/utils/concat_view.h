#ifndef RB_CONCAT_VIEW_H
#define RB_CONCAT_VIEW_H

#include <concepts>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

namespace RevBayesCore {

namespace detail {

// Limit concatenation to forward ranges of the same value type with ordinary, compatible references.
template <class... R>
concept concat_compatible_ranges = sizeof...(R) > 0 && (std::ranges::forward_range<R> && ...) &&
    requires { typename std::common_reference_t<std::ranges::range_reference_t<R>...>; } &&
    (std::same_as<std::ranges::range_value_t<R>,
                  std::common_type_t<std::ranges::range_value_t<R>...>> && ...) &&
    (std::is_lvalue_reference_v<std::ranges::range_reference_t<R>> && ...) &&
    (std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<R>>,
                  std::ranges::range_value_t<R>> && ...) &&
    std::is_lvalue_reference_v<std::common_reference_t<std::ranges::range_reference_t<R>...>> &&
    (std::convertible_to<std::ranges::range_reference_t<R>,
                         std::common_reference_t<std::ranges::range_reference_t<R>...>> && ...);

}

// A C++23 substitute for the forward-traversal part of C++26 concat_view:
// https://eel.is/c++draft/range.concat
// Adapt inputs with views::all: borrow lvalue containers, own movable temporary containers, and
// copy or move existing views. Owning a non-owning view does not extend its underlying data's
// lifetime. Input iterator-invalidation rules still apply; moving or destroying this view
// invalidates its iterators because they refer back to it. Traversal does not copy elements or allocate.
// Unlike the standard view, this omits proxy references, reverse traversal, and indexing.
// Replace it when C++26 concat is available.
template <std::ranges::view... Views>
    requires detail::concat_compatible_ranges<Views...>
class concat_view : public std::ranges::view_interface<concat_view<Views...>> {
    std::tuple<Views...> ranges_;

public:
    template <bool Const>
    class iterator {
        friend class concat_view;

        template <class T>
        using maybe_const = std::conditional_t<Const, const T, T>;

        maybe_const<concat_view>* parent_ = nullptr;
        std::variant<std::ranges::iterator_t<maybe_const<Views>>...> current_;
        static constexpr std::size_t last = sizeof...(Views) - 1;

        // Start in the first range, skipping it and any following empty ranges.
        explicit iterator(maybe_const<concat_view>* parent)
            : parent_(parent), current_(std::in_place_index<0>, std::ranges::begin(std::get<0>(parent->ranges_)))
        {
            satisfy<0>();
        }

        // The variant index identifies the active range, including when iterator types repeat.
        // On exhaustion, visit successive ranges in order until one has an element or the final end
        // is reached. This visits every element once, retaining duplicates and skipping empty ranges.
        template <std::size_t I>
        void satisfy()
        {
            if constexpr (I < last)
                if (std::get<I>(current_) == std::ranges::end(std::get<I>(parent_->ranges_)))
                {
                    current_.template emplace<I + 1>(std::ranges::begin(std::get<I + 1>(parent_->ranges_)));
                    satisfy<I + 1>();
                }
        }

        // Dispatch by index so each iterator is paired with its own range's end, even for equal types.
        template <std::size_t I = 0>
        void advance()
        {
            if (current_.index() == I)
            {
                ++std::get<I>(current_);
                satisfy<I>();
            }
            else if constexpr (I < last)
                advance<I + 1>();
        }

    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::common_type_t<std::ranges::range_value_t<maybe_const<Views>>...>;
        using difference_type = std::common_type_t<std::ranges::range_difference_t<maybe_const<Views>>...>;
        using reference = std::common_reference_t<std::ranges::range_reference_t<maybe_const<Views>>...>;

        iterator() = default;

        // Preserve element references while allowing the inputs to differ in constness.
        reference operator*() const
        {
            return std::visit([](const auto& it) -> reference { return *it; }, current_);
        }

        // Advance the active iterator, crossing range boundaries as needed.
        iterator& operator++()
        {
            advance();
            return *this;
        }

        // Return an independent copy at the previous position, as required for forward iteration.
        iterator operator++(int)
        {
            auto previous = *this;
            ++*this;
            return previous;
        }

        bool operator==(const iterator&) const = default;

        // Earlier ends are never stopping points: satisfy() has already crossed them.
        bool operator==(std::default_sentinel_t) const
        {
            return current_.index() == last &&
                   std::get<last>(current_) == std::ranges::end(std::get<last>(parent_->ranges_));
        }
    };

    explicit concat_view(Views... ranges) : ranges_(std::move(ranges)...) {}

    iterator<false> begin() { return iterator<false>(this); }

    // Const borrowing views can expose mutable elements; const owning views generally cannot.
    // Only instantiate the const iterator when every const input still meets our range requirements.
    auto begin() const requires detail::concat_compatible_ranges<const Views...>
    {
        return iterator<true>(this);
    }

    std::default_sentinel_t end() const { return {}; }
};

template <std::ranges::viewable_range... R>
concat_view(R&&...) -> concat_view<std::views::all_t<R>...>;

}

#endif
