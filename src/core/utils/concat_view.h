#ifndef RB_CONCAT_VIEW_H
#define RB_CONCAT_VIEW_H

#include <concepts>
#include <iterator>
#include <memory>
#include <ranges>
#include <tuple>
#include <type_traits>
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
// Borrow one or more lvalue ranges; neither containers nor elements are copied. All ranges must
// outlive traversal, and their iterator-invalidation rules apply. Unlike the standard view, this
// omits temporary ownership, proxy references, reverse traversal, and indexing. Replace it when
// C++26 concat is available.
template <class... R>
    requires detail::concat_compatible_ranges<R...>
class concat_view : public std::ranges::view_interface<concat_view<R...>> {
    std::tuple<R*...> ranges_;

public:
    class iterator {
        friend class concat_view;

        const concat_view* parent_ = nullptr;
        std::variant<std::ranges::iterator_t<R>...> current_;
        static constexpr std::size_t last = sizeof...(R) - 1;

        // Start in the first range, skipping it and any following empty ranges.
        explicit iterator(const concat_view* parent)
            : parent_(parent), current_(std::in_place_index<0>, std::ranges::begin(*std::get<0>(parent->ranges_)))
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
                if (std::get<I>(current_) == std::ranges::end(*std::get<I>(parent_->ranges_)))
                {
                    current_.template emplace<I + 1>(std::ranges::begin(*std::get<I + 1>(parent_->ranges_)));
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
        using value_type = std::common_type_t<std::ranges::range_value_t<R>...>;
        using difference_type = std::common_type_t<std::ranges::range_difference_t<R>...>;
        using reference = std::common_reference_t<std::ranges::range_reference_t<R>...>;

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
                   std::get<last>(current_) == std::ranges::end(*std::get<last>(parent_->ranges_));
        }
    };

    // Deducing lvalue arguments rejects temporaries even with explicit const range template arguments.
    template <class... Args>
        requires std::same_as<std::tuple<Args...>, std::tuple<R...>>
    explicit concat_view(Args&... ranges) : ranges_(std::addressof(ranges)...) {}

    iterator begin() const { return iterator(this); }
    std::default_sentinel_t end() const { return {}; }
};

template <class... R>
concat_view(R&...) -> concat_view<R...>;

}

#endif
