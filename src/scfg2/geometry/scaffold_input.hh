#ifndef CPARTY_SCFG2_GEOMETRY_SCAFFOLD_INPUT_HH
#define CPARTY_SCFG2_GEOMETRY_SCAFFOLD_INPUT_HH

#include <cstddef>
#include <string>
#include <string_view>

namespace cparty::scfg2::geometry {

inline bool validate_scaffold_input(std::string_view structure,
                                    std::string *error) {
  std::size_t depth = 0;
  for (const char symbol : structure) {
    if (symbol == '(') {
      ++depth;
    } else if (symbol == ')') {
      if (!depth) {
        if (error) *error = "unbalanced scaffold round brackets";
        return false;
      }
      --depth;
    } else if (symbol != '.' && symbol != '_' && symbol != 'x' && symbol != 'X') {
      if (error) {
        *error = "unsupported scaffold symbol: use round brackets '()' for "
                 "a noncrossing scaffold, '.'/'_' for free bases, or 'x'/'X' "
                 "for forced unpaired bases";
      }
      return false;
    }
  }
  if (depth) {
    if (error) *error = "unbalanced scaffold round brackets";
    return false;
  }
  return true;
}

}  // namespace cparty::scfg2::geometry

#endif
