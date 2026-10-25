#pragma once

#include "combinatorial_hex.hpp"
#include <string>

namespace jw::grid {
    // Display projection lives here; the patch and substitution core stay integer-only.
    void to_svg(const std::string& path,const patch& value,bool show_arrows=false);
    void rules_svg(const std::string& path,bool show_arrows=false);
}
