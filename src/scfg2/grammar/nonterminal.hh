#ifndef CPARTY_SCFG2_GRAMMAR_NONTERMINAL_HH
#define CPARTY_SCFG2_GRAMMAR_NONTERMINAL_HH

#include <string_view>

namespace cparty::scfg2::grammar {

enum class NonTerminal : unsigned char {
  W,
  WI,
  V,
  VM,
  WMv,
  WMp,
  WM,
  WIP,
  VPL,
  VPR,
  VP,
  WMBW,
  WMBP,
  WMB,
  BE,
};

inline constexpr std::string_view nonterminal_name(NonTerminal nt) {
  switch (nt) {
    case NonTerminal::W: return "W";
    case NonTerminal::WI: return "WI";
    case NonTerminal::V: return "V";
    case NonTerminal::VM: return "VM";
    case NonTerminal::WMv: return "WMv";
    case NonTerminal::WMp: return "WMp";
    case NonTerminal::WM: return "WM";
    case NonTerminal::WIP: return "WIP";
    case NonTerminal::VPL: return "VPL";
    case NonTerminal::VPR: return "VPR";
    case NonTerminal::VP: return "VP";
    case NonTerminal::WMBW: return "WMBW";
    case NonTerminal::WMBP: return "WMBP";
    case NonTerminal::WMB: return "WMB";
    case NonTerminal::BE: return "BE";
  }
  return "UNKNOWN";
}

inline constexpr bool parse_nonterminal_name(std::string_view name, NonTerminal *out) {
  if (!out) {
    return false;
  }
  if (name == "W") *out = NonTerminal::W;
  else if (name == "WI") *out = NonTerminal::WI;
  else if (name == "V") *out = NonTerminal::V;
  else if (name == "VM") *out = NonTerminal::VM;
  else if (name == "WMv") *out = NonTerminal::WMv;
  else if (name == "WMp") *out = NonTerminal::WMp;
  else if (name == "WM") *out = NonTerminal::WM;
  else if (name == "WIP") *out = NonTerminal::WIP;
  else if (name == "VPL") *out = NonTerminal::VPL;
  else if (name == "VPR") *out = NonTerminal::VPR;
  else if (name == "VP") *out = NonTerminal::VP;
  else if (name == "WMBW") *out = NonTerminal::WMBW;
  else if (name == "WMBP") *out = NonTerminal::WMBP;
  else if (name == "WMB") *out = NonTerminal::WMB;
  else if (name == "BE") *out = NonTerminal::BE;
  else return false;
  return true;
}

}  // namespace cparty::scfg2::grammar

#endif  // CPARTY_SCFG2_GRAMMAR_NONTERMINAL_HH
