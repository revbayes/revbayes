#ifndef IndexedCache_H
#define IndexedCache_H

#include <cassert>
#include <cstddef>
#include <optional>
#include <vector>

namespace RevBayesCore {

    template <typename T>
    class IndexedCache {
        struct ItemState {
            unsigned active : 1;   // which slot: 0 or 1
            unsigned dirty  : 1;   // does this slot need recomputation?
        };

        size_t num_items;
        std::vector<T> slots;               // size = 2 * num_items
        std::vector<ItemState> current_state;
        std::optional<std::vector<ItemState>> prev_state;

        // Return the storage slot for an item and active-slot index.
        T& slot(size_t item, unsigned s) {
            return slots[s * num_items + item];
        }

    public:
        // Construct an indexed two-slot cache with all items initially dirty.
        IndexedCache(size_t n)
            : num_items(n),
              slots(2 * n),
              current_state(n, {0, 1})   // all start as slot 0, dirty
            {}

        size_t size() const { return num_items;}

        // Return whether the current slot for an item needs recomputation.
        bool is_dirty(size_t item) const
        {
            assert( item >= 0 and item < num_items );

            return current_state[item].dirty;
        }

        // Return read-only access to the current clean slot for an item.
        const T& operator[](size_t item) const
        {
            assert( not is_dirty(item) );

            return slots[current_state[item].active * num_items + item];
        }

        // Return mutable access to the current clean slot for an item.
        T& get_mutable_item(size_t item)
        {
            assert( not is_dirty(item) );

            return slots[current_state[item].active * num_items + item];
        }

        // Mark an item dirty and detach its active slot from the saved state when needed.
        void mark_dirty(size_t item)
        {
            assert( item >= 0 and item < num_items );

            // Should we assume that we only mark things dirty after touching?
            // assert( is_touched() );

            // Here is where we unshare with the previous state.
            if (prev_state and not (*prev_state)[item].dirty and not is_dirty(item))
            {
                // If we mark something dirty multiple times, we don't want to flip the bit twice!
                current_state[item].active = (*prev_state)[item].active ^ 1;
            }

            current_state[item].dirty = 1;
        }

        bool is_touched() const
        {
            return prev_state.has_value();
        }

        // Mark every indexed item dirty.
        void mark_all_dirty()
        {
            for(size_t i = 0; i < num_items; i++)
                mark_dirty(i);
        }

        // Get a mutable reference to dirty slot and mark it clean.
        T& init_for_writing(size_t item)
        {
            // Don't over-write something that hasn't been invalidated.
            assert(is_dirty(item));

            // Mark the item clean.
            current_state[item].dirty = 0;

            // Access the now-clean item.
            return get_mutable_item(item);
        }

        // Accept the current cache metadata and discard the saved rollback state.
        void keep()
        {
            // We should only call keep if we proposed a new state and are accepting it.
            assert(prev_state);

            prev_state.reset();
        }

        // Restore the saved cache metadata after rejecting a proposal.
        void restore()
        {
            // Moving back to the previous state makes no sense if there is no previous state.
            assert(prev_state);

            current_state = *prev_state;
            prev_state.reset();
        }

        // Save active-slot and dirty metadata for proposal rollback.
        void touch()
        {
            if (not prev_state)
                prev_state = current_state;
        }

        // Resize after tree topology changes and reset all items dirty.
        void resize(size_t n)
        {
            num_items = n;
            slots.resize(2*n);
            current_state.assign(n, {0, 1});

            // if we change the number of nodes, how would we restore?
            if (prev_state)
                prev_state = current_state;
        }

        // Remove all cached slots and reset metadata to an empty saved state.
        void clear()
        {
            slots.clear();
            current_state.clear();
            num_items = 0;
            prev_state = current_state;
        }
    };

}

#endif
