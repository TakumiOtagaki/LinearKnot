#ifndef CPARTY_SCFG2_CLI_CLI_COMMON_HH
#define CPARTY_SCFG2_CLI_CLI_COMMON_HH

// Input handling shared by the LinearKnot subcommands.

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace cparty::scfg2::cli {

struct Record {
  std::string name;
  std::string sequence;
  std::string scaffold;
};

inline bool parse_lengths(const std::string &raw, std::vector<int> *out) {
  std::stringstream ss(raw);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (token.empty()) return false;
    char *end = nullptr;
    const long value = std::strtol(token.c_str(), &end, 10);
    if (*end != '\0' || value < 1) return false;
    out->push_back(static_cast<int>(value));
  }
  return !out->empty();
}

inline std::string trim(const std::string &s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

inline bool read_records(const std::string &path, std::vector<Record> *out,
                  std::string *error) {
  std::ifstream in(path);
  if (!in) {
    *error = "cannot open " + path;
    return false;
  }
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(in, line)) {
    line = trim(line);
    if (!line.empty() && line[0] != '#') lines.push_back(line);
  }
  for (std::size_t at = 0; at < lines.size();) {
    Record record;
    if (lines[at][0] == '>') {
      record.name = trim(lines[at].substr(1));
      ++at;
    } else {
      record.name = "seq" + std::to_string(out->size() + 1);
    }
    if (at + 1 >= lines.size()) {
      *error = "record '" + record.name + "' lacks a sequence or scaffold line";
      return false;
    }
    record.sequence = lines[at++];
    record.scaffold = lines[at++];
    for (auto &c : record.sequence) {
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      if (c == 'T') c = 'U';
    }
    out->push_back(record);
  }
  return true;
}

}  // namespace cparty::scfg2::cli

#endif  // CPARTY_SCFG2_CLI_CLI_COMMON_HH
