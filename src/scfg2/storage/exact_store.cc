#include "scfg2/storage/exact_store.hh"

#include <unordered_map>
#include <unordered_set>

namespace cparty::scfg2::storage {
namespace {

template <typename T>
void hash_combine(std::size_t &seed, const T &value) {
  seed ^= std::hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct ItemSpanKey {
  NonTerminal nt = NonTerminal::W;
  cand_pos_t i = -1;
  cand_pos_t j = -1;

  bool operator==(const ItemSpanKey &other) const {
    return nt == other.nt && i == other.i && j == other.j;
  }
};

struct ItemSpanKeyHash {
  std::size_t operator()(const ItemSpanKey &key) const {
    std::size_t seed = 0;
    hash_combine(seed, static_cast<unsigned char>(key.nt));
    hash_combine(seed, key.i);
    hash_combine(seed, key.j);
    return seed;
  }
};

using TotalByItemMap = std::unordered_map<ExactItemKey, pf_t, ExactItemKeyHash>;
using ItemSet = std::unordered_set<ExactItemKey, ExactItemKeyHash>;
using ItemSetBySpanMap =
    std::unordered_map<ItemSpanKey, ItemSet, ItemSpanKeyHash>;

thread_local TotalByItemMap total_by_item;
thread_local ItemSetBySpanMap items_by_span;

ItemSpanKey item_span_key(const ExactItemKey &item) {
  return {item.nt, item.i, item.j};
}

}  // namespace

bool ExactItemKey::operator==(const ExactItemKey &other) const {
  return nt == other.nt && i == other.i && j == other.j && ip == other.ip &&
         jp == other.jp && aux0 == other.aux0 && aux1 == other.aux1 &&
         aux2 == other.aux2 && aux3 == other.aux3;
}

std::size_t ExactItemKeyHash::operator()(const ExactItemKey &key) const {
  std::size_t seed = 0;
  hash_combine(seed, static_cast<unsigned char>(key.nt));
  hash_combine(seed, key.i);
  hash_combine(seed, key.j);
  hash_combine(seed, key.ip);
  hash_combine(seed, key.jp);
  hash_combine(seed, key.aux0);
  hash_combine(seed, key.aux1);
  hash_combine(seed, key.aux2);
  hash_combine(seed, key.aux3);
  return seed;
}

namespace exact_store {

void reset() {
  total_by_item.clear();
  items_by_span.clear();
}

void set_total(const ExactItemKey &item, pf_t value) {
  if (value == 0) {
    erase(item);
    return;
  }
  total_by_item[item] = value;
  items_by_span[item_span_key(item)].insert(item);
}

void erase(const ExactItemKey &item) {
  total_by_item.erase(item);
  const ItemSpanKey span_key = item_span_key(item);
  const auto found = items_by_span.find(span_key);
  if (found == items_by_span.end()) return;
  found->second.erase(item);
  if (found->second.empty()) items_by_span.erase(found);
}

pf_t total(const ExactItemKey &item) {
  const auto found = total_by_item.find(item);
  return found == total_by_item.end() ? 0 : found->second;
}

void for_each_item_span(NonTerminal nt,
                        cand_pos_t i,
                        cand_pos_t j,
                        const VisitFn &visit) {
  const auto found = items_by_span.find({nt, i, j});
  if (found == items_by_span.end()) return;
  for (const auto &item : found->second) {
    const auto total_it = total_by_item.find(item);
    if (total_it != total_by_item.end()) visit(item, total_it->second);
  }
}

}  // namespace exact_store
}  // namespace cparty::scfg2::storage
