#ifndef RB_CONCAT_VIEW_H
#define RB_CONCAT_VIEW_H

#include <concepts>
#include <iterator>
#include <memory>
#include <ranges>
#include <type_traits>
#include <variant>

namespace RevBayesCore {

// A C++23 substitute for the forward-traversal part of C++26 concat_view:
// https://eel.is/c++draft/range.concat
// Borrow two lvalue ranges with ordinary element references; neither containers nor elements
// are copied. Both ranges must outlive traversal, and their iterator-invalidation rules apply.
// Unlike the standard view, this supports exactly two ranges and no temporary ownership,
// proxy references, reverse traversal, or indexing. Replace it when C++26 concat is available.
template <std::ranges::forward_range R1, std::ranges::forward_range R2>
    requires std::same_as<std::ranges::range_value_t<R1>, std::ranges::range_value_t<R2>> &&
             std::is_lvalue_reference_v<std::ranges::range_reference_t<R1>> &&
             std::is_lvalue_reference_v<std::ranges::range_reference_t<R2>> &&
             std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<R1>>,
                          std::ranges::range_value_t<R1>> &&
             std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<R2>>,
                          std::ranges::range_value_t<R2>> &&
             std::common_reference_with<std::ranges::range_reference_t<R1>,
                                        std::ranges::range_reference_t<R2>> &&
             std::is_lvalue_reference_v<std::common_reference_t<std::ranges::range_reference_t<R1>,
                                                               std::ranges::range_reference_t<R2>>>
class concat_view : public std::ranges::view_interface<concat_view<R1, R2>> {
    R1* first_;
    R2* second_;

public:
    class iterator {
        friend class concat_view;

        const concat_view* parent_ = nullptr;
        std::variant<std::ranges::iterator_t<R1>, std::ranges::iterator_t<R2>> current_;

        // Start in the first range, switching immediately if it is empty.
        explicit iterator(const concat_view* parent)
            : parent_(parent), current_(std::in_place_index<0>, std::ranges::begin(*parent->first_))
        {
            satisfy();
        }

        // The variant index records the active range. Once the first range is exhausted,
        // switch to the second; its end represents completion, including an empty second range.
        // Thus each element is visited once in input order, including duplicates in either range.
        void satisfy()
        {
            if (current_.index() == 0 && std::get<0>(current_) == std::ranges::end(*parent_->first_))
                current_.template emplace<1>(std::ranges::begin(*parent_->second_));
        }

    public:
        using iterator_concept = std::forward_iterator_tag;
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::ranges::range_value_t<R1>;
        using difference_type = std::common_type_t<std::ranges::range_difference_t<R1>,
                                                   std::ranges::range_difference_t<R2>>;
        using reference = std::common_reference_t<std::ranges::range_reference_t<R1>,
                                                  std::ranges::range_reference_t<R2>>;

        iterator() = default;

        // Preserve element references while allowing the inputs to differ in constness.
        reference operator*() const
        {
            if (current_.index() == 0)
                return *std::get<0>(current_);
            return *std::get<1>(current_);
        }

        // Advance the active iterator and cross the boundary only after exhausting the first range.
        iterator& operator++()
        {
            if (current_.index() == 0)
            {
                ++std::get<0>(current_);
                satisfy();
            }
            else
                ++std::get<1>(current_);
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

        // The first range's end is never a stopping point: satisfy() has already crossed it.
        bool operator==(std::default_sentinel_t) const
        {
            return current_.index() == 1 && std::get<1>(current_) == std::ranges::end(*parent_->second_);
        }
    };

    concat_view(R1& first, R2& second) : first_(std::addressof(first)), second_(std::addressof(second)) {}

    // Reject temporaries even when explicitly supplied template arguments make a range const.
    concat_view(R1&&, R2&) = delete;
    concat_view(R1&, R2&&) = delete;
    concat_view(R1&&, R2&&) = delete;

    iterator begin() const { return iterator(this); }
    std::default_sentinel_t end() const { return {}; }
};

template <class R1, class R2>
concat_view(R1&, R2&) -> concat_view<R1, R2>;

}

#endif
