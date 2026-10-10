// EHeap: a binary heap of elements with keys, smallest key on top, as
// GtsEHeap. Elements keep their place in the heap through handles, for
// remove() and decrease_key(). The heap's layout follows GtsEHeap's exactly,
// so that elements with equal keys come off in the same order as in GTS for
// the same sequence of operations. GtsEHeap's randomized mode is not ported.
//
// Copyright (C) 1999 Stéphane Popinet (GTS: eheap.c)
// Copyright (C) 2026 Matevz Tadel
// Rewritten in C++ from GTS on 2026-10-09.
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace gts {

template <class T>
class EHeap
{
public:
  using Handle  = std::uint32_t;
  using KeyFunc = std::function<double(const T&)>;

  explicit EHeap(KeyFunc key = {}) : m_key(std::move(key)) {}

  // Inserts p with the key from the key function, or with the given key.
  Handle insert(T p) { const double k = m_key(p); return insert(std::move(p), k); }
  Handle insert(T p, double key);

  // Removes the element on top and returns it with its key.
  std::optional<std::pair<T, double>> remove_top();
  // The element on top and its key, without removing it.
  std::optional<std::pair<T, double>> top() const;

  // Removes the element of h. h is invalid afterwards.
  T    remove(Handle h);
  // Lowers the key of h's element; new_key must not be larger than the key.
  void decrease_key(Handle h, double new_key);

  const T& data(Handle h) const { return m_pairs[h].data; }
  double   key(Handle h) const  { return m_pairs[h].key; }
  double   key_of(const T& p) const { return m_key(p); }

  std::size_t size() const { return m_heap.size(); }
  bool        empty() const { return m_heap.empty(); }

  // While frozen, insertions do not keep the heap ordered; thaw() orders it
  // in O(n). This builds a heap from n elements in O(n).
  void freeze() { m_frozen = true; }
  void thaw();
  // Recomputes every key with the key function and reorders.
  void update();

  // Calls f(element) for each element, in the order of the heap's array.
  template <class F> void foreach(F&& f) const { for (Handle h : m_heap) f(m_pairs[h].data); }

private:
  struct Pair
  {
    T             data;
    double        key;
    std::uint32_t pos; // 1-based position in m_heap
  };

  Pair& at(std::uint32_t pos) { return m_pairs[m_heap[pos - 1]]; }
  void  place(std::uint32_t pos, Handle h) { m_heap[pos - 1] = h; m_pairs[h].pos = pos; }
  void  release(Handle h) { m_free.push_back(h); }
  static std::uint32_t parent(std::uint32_t i) { return i >= 2 ? i / 2 : 0; }

  void sift_up(std::uint32_t i);
  void sift_down(std::uint32_t i);

  KeyFunc             m_key;
  std::vector<Pair>   m_pairs;
  std::vector<Handle> m_free;
  std::vector<Handle> m_heap;
  bool                m_frozen = false;
};

template <class T>
typename EHeap<T>::Handle EHeap<T>::insert(T p, double key)
{
  Handle h;
  if (m_free.empty()) { h = Handle(m_pairs.size()); m_pairs.push_back({ std::move(p), key, 0 }); }
  else                { h = m_free.back(); m_free.pop_back(); m_pairs[h] = { std::move(p), key, 0 }; }
  m_heap.push_back(h);
  m_pairs[h].pos = std::uint32_t(m_heap.size());
  if ( ! m_frozen) sift_up(std::uint32_t(m_heap.size()));
  return h;
}

template <class T>
void EHeap<T>::sift_up(std::uint32_t i)
{
  const Handle child = m_heap[i - 1];
  const double key = m_pairs[child].key;
  std::uint32_t p;
  while ((p = parent(i)))
  {
    const Handle par = m_heap[p - 1];
    if (m_pairs[par].key > key)
    {
      place(p, child);
      place(i, par);
      i = p;
    }
    else
      break;
  }
}

template <class T>
void EHeap<T>::sift_down(std::uint32_t i)
{
  const std::uint32_t len = std::uint32_t(m_heap.size());
  const Handle par = m_heap[i - 1];
  const double key = m_pairs[par].key;
  for (;;)
  {
    const std::uint32_t lc = 2 * i, rc = 2 * i + 1;
    if (lc > len) break;
    std::uint32_t c = lc;
    if (rc <= len && ! (at(lc).key < at(rc).key)) c = rc;
    const Handle child = m_heap[c - 1];
    if ( ! (key > m_pairs[child].key)) break;
    place(i, child);
    place(c, par);
    i = c;
  }
}

template <class T>
std::optional<std::pair<T, double>> EHeap<T>::remove_top()
{
  if (m_heap.empty()) return std::nullopt;
  const Handle h = m_heap.front();
  std::pair<T, double> r(std::move(m_pairs[h].data), m_pairs[h].key);
  release(h);
  const Handle last = m_heap.back();
  m_heap.pop_back();
  if ( ! m_heap.empty())
  {
    place(1, last);
    sift_down(1);
  }
  return r;
}

template <class T>
std::optional<std::pair<T, double>> EHeap<T>::top() const
{
  if (m_heap.empty()) return std::nullopt;
  const Pair& p = m_pairs[m_heap.front()];
  return std::pair<T, double>(p.data, p.key);
}

// As gts_eheap_remove: moves the element to the top, then removes the top.
template <class T>
T EHeap<T>::remove(Handle h)
{
  std::uint32_t i = m_pairs[h].pos, p;
  while ((p = parent(i)))
  {
    const Handle par = m_heap[p - 1];
    place(p, h);
    place(i, par);
    i = p;
  }
  return std::move(remove_top()->first);
}

template <class T>
void EHeap<T>::decrease_key(Handle h, double new_key)
{
  m_pairs[h].key = new_key;
  if ( ! m_frozen) sift_up(m_pairs[h].pos);
}

template <class T>
void EHeap<T>::thaw()
{
  if ( ! m_frozen) return;
  for (std::uint32_t i = std::uint32_t(m_heap.size()) / 2; i > 0; --i)
    sift_down(i);
  m_frozen = false;
}

template <class T>
void EHeap<T>::update()
{
  m_frozen = true;
  for (Handle h : m_heap) m_pairs[h].key = m_key(m_pairs[h].data);
  thaw();
}

} // namespace gts
