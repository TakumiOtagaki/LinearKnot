#include "scfg2/cli/commands.hh"

#include <cstdlib>
#include <iostream>
#include <string_view>

namespace {

struct Command {
  std::string_view name;
  const char *description;
  int (*run)(int, char **);
};

constexpr Command commands[] = {
    {"fold", "Minimum-energy structure, ensemble energy and pair probabilities",
     cparty::scfg2::cli::run_fold},
    {"access", "Window accessibility and opening energies",
     cparty::scfg2::cli::run_access},
    {"profile", "Per-base structural context probabilities",
     cparty::scfg2::cli::run_profile},
    {"sample", "Structures sampled from the ensemble",
     cparty::scfg2::cli::run_sample},
    {"subopt", "Ranked suboptimal structures",
     cparty::scfg2::cli::run_subopt},
    {"mea", "Maximum expected accuracy prediction",
     cparty::scfg2::cli::run_mea},
    {"centroid", "Gamma-centroid prediction",
     cparty::scfg2::cli::run_centroid},
};

void usage(std::ostream &os) {
  os << "LinearKnot: hierarchical RNA folding and ensemble analysis\n\n"
        "Usage: linearknot COMMAND [options]\n\n"
        "Commands:\n";
  for (const auto &command : commands)
    os << "  " << command.name << "\t" << command.description << '\n';
  os << "\nUse 'linearknot COMMAND --help' for command options.\n"
        "  -h, --help  Show this help\n"
        "  --version   Show the LinearKnot version\n";
}

const Command *find_command(std::string_view name) {
  for (const auto &command : commands)
    if (command.name == name) return &command;
  return nullptr;
}

int unknown_command(std::string_view name) {
  std::cerr << "error: unknown command: " << name << '\n';
  usage(std::cerr);
  return EXIT_FAILURE;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    usage(std::cerr);
    return EXIT_FAILURE;
  }
  const std::string_view name = argv[1];
  if (name == "--help" || name == "-h") {
    usage(std::cout);
    return EXIT_SUCCESS;
  }
  if (name == "--version") {
    std::cout << "LinearKnot " << LINEARKNOT_VERSION << '\n';
    return EXIT_SUCCESS;
  }
  if (name == "help") {
    if (argc == 2) {
      usage(std::cout);
      return EXIT_SUCCESS;
    }
    if (argc != 3) {
      std::cerr << "error: usage: linearknot help [COMMAND]\n";
      return EXIT_FAILURE;
    }
    const auto *command = find_command(argv[2]);
    if (!command) return unknown_command(argv[2]);
    char help[] = "--help";
    char *arguments[] = {argv[2], help, nullptr};
    return command->run(2, arguments);
  }
  const auto *command = find_command(name);
  if (!command) return unknown_command(name);
  return command->run(argc - 1, argv + 1);
}
