#ifndef CPARTY_SCFG2_BEAM_MEMORY_OBSERVATION_HH
#define CPARTY_SCFG2_BEAM_MEMORY_OBSERVATION_HH

#include "scfg2/grammar/nonterminal.hh"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace cparty::scfg2::beam {

// Opt-in observations only: no allocation policy, pruning, or weights change.
// Byte counts cover stored values, vector capacities and hash bucket pointers;
// they exclude hash-node links, allocator headers, alignment and fragmentation.
struct ContainerMemory {
  std::string name;
  std::size_t entries = 0;
  std::size_t buckets = 0;
  std::size_t value_bytes = 0;
  std::size_t vector_size_bytes = 0;
  std::size_t vector_capacity_bytes = 0;

  std::size_t accounted_bytes() const {
    return value_bytes + buckets * sizeof(void *) + vector_capacity_bytes;
  }

  void add(const ContainerMemory &other) {
    entries += other.entries;
    buckets += other.buckets;
    value_bytes += other.value_bytes;
    vector_size_bytes += other.vector_size_bytes;
    vector_capacity_bytes += other.vector_capacity_bytes;
  }
};

enum class MemoryPhase { Scaffold, Pilot, Outer, Main };
inline const char *memory_phase_name(MemoryPhase phase) {
  switch (phase) {
    case MemoryPhase::Scaffold: return "scaffold";
    case MemoryPhase::Pilot: return "pilot";
    case MemoryPhase::Outer: return "outer";
    case MemoryPhase::Main: return "main";
  }
  return "unknown";
}

struct FamilyMemory {
  grammar::NonTerminal nonterminal;
  std::size_t live_items = 0;
};

struct MemoryObservation {
  MemoryPhase phase = MemoryPhase::Main;
  const char *event = "";
  std::size_t width = 0;
  std::size_t column = 0;
  std::size_t chart_items = 0;
  std::vector<ContainerMemory> containers;
  // Populated only at forward completion, after public pruning.
  std::vector<FamilyMemory> families;
};
using MemoryObserver = std::function<void(const MemoryObservation &)>;

template <class Vector>
ContainerMemory vector_memory(const char *name, const Vector &values) {
  using Cell = typename Vector::value_type;
  return {name, values.size(), 0, 0,
          values.size() * sizeof(Cell), values.capacity() * sizeof(Cell)};
}

template <class Vector>
ContainerMemory nested_vector_memory(const char *name, const Vector &values) {
  auto out = vector_memory(name, values);
  for (const auto &value : values) out.add(vector_memory("", value));
  return out;
}

template <class Map>
ContainerMemory hash_memory(const char *name, const Map &values) {
  return {name, values.size(), values.bucket_count(),
          values.size() * sizeof(typename Map::value_type), 0, 0};
}

template <class Map>
ContainerMemory hash_vector_memory(const char *name, const Map &values) {
  auto out = hash_memory(name, values);
  for (const auto &[key, vector] : values) {
    (void)key;
    using Cell = typename Map::mapped_type::value_type;
    out.vector_size_bytes += vector.size() * sizeof(Cell);
    out.vector_capacity_bytes += vector.capacity() * sizeof(Cell);
  }
  return out;
}

template <class Chart>
MemoryObservation chart_memory(MemoryPhase phase, const char *event,
                               std::size_t width, std::size_t column,
                               const Chart &chart, bool families = false) {
  MemoryObservation out;
  out.phase = phase;
  out.event = event;
  out.width = width;
  out.column = column;
  out.chart_items = chart.values().size();
  out.containers.push_back(hash_memory("chart", chart.values()));
  if (families) {
    const auto count = static_cast<std::size_t>(grammar::NonTerminal::BE) + 1;
    out.families.resize(count);
    for (std::size_t i = 0; i < count; ++i)
      out.families[i].nonterminal = static_cast<grammar::NonTerminal>(i);
    for (const auto &[item, value] : chart.values()) {
      (void)value;
      ++out.families[static_cast<std::size_t>(item.nonterminal)].live_items;
    }
  }
  return out;
}

} // namespace cparty::scfg2::beam
#endif
