#ifndef RbHelpJson_H
#define RbHelpJson_H

#include <string>
#include <vector>

namespace RevBayesCore {

    /**
     * Render every help entry, with its arguments, types and references, as a JSON array.
     * The global workspace must be initialized first.
     */
    std::string renderHelpJson(void);

    /**
     * Malformed help pages in help/md, one message each.
     */
    std::vector<std::string> checkHelpPages(void);

}

#endif
