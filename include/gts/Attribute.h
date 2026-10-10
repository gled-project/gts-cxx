// Attributes: data attached to vertices, edges or faces by index.
//
// An attribute is an array parallel to the vertex, edge or face array of a
// Geometrium. It grows with the elements and is renumbered with them by
// collect(). Hooks say how a value is made when an algorithm makes an element
// from others; this takes the place of subclassing the elements. Each type
// has default hooks, which a specialization of DefaultHooks can change for
// all attributes of that type, and Hooks can change for one attribute.
//
// Copyright (C) 2026 Matevz Tadel
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "gts/Elements.h"

#include <concepts>
#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <typeindex>
#include <vector>

namespace gts {

// Values that can be interpolated linearly: floating-point numbers, and
// classes with a + t * (b - a), such as Vec3. Integers are not: they are
// usually labels, counts or flags.
template <class T>
concept Affine = (std::floating_point<T> || std::is_class_v<T>) &&
  requires(const T& a, const T& b, double t) { { a + t * (b - a) } -> std::convertible_to<T>; };

// The hooks used for a type when Hooks<T> leaves one unset. An affine value
// is interpolated linearly and merged to the middle; any other value takes
// the first source. A piece of a split element copies its parent.
template <class T>
struct GenericHooks
{
  static T interpolate(const T& a, const T& b, double t)
  {
    if constexpr (Affine<T>) return T(a + t * (b - a));
    else                     { (void) b; (void) t; return a; }
  }
  static T merge(const T& a, const T& b) { return interpolate(a, b, 0.5); }
  static T split(const T& parent)        { return parent; }
};

// The customization point: specialize DefaultHooks for a type of your own,
// deriving from GenericHooks to keep the hooks you do not replace.
//
//   template <> struct gts::DefaultHooks<Colour> : gts::GenericHooks<Colour>
//   { static Colour merge(const Colour& a, const Colour& b) { return mix(a, b); } };
template <class T>
struct DefaultHooks : GenericHooks<T> {};

// How a value is made when an algorithm makes an element from others, for
// one attribute. A hook that is not set comes from DefaultHooks<T>.
//
//   gts::Hooks<float>().on_interpolate([](float a, float b, double t) { return a + t * (b - a); })
template <class T>
struct Hooks
{
  std::function<T(const T& a, const T& b, double t)> interpolate; // an element between a and b
  std::function<T(const T& a, const T& b)>           merge;       // a and b collapse into one
  std::function<T(const T& parent)>                  split;       // a piece of a split element

  Hooks&& on_interpolate(decltype(interpolate) f) && { interpolate = std::move(f); return std::move(*this); }
  Hooks&& on_merge(decltype(merge) f) &&             { merge = std::move(f);       return std::move(*this); }
  Hooks&& on_split(decltype(split) f) &&             { split = std::move(f);       return std::move(*this); }
};

// A typed handle on one attribute. Valid as long as its Geometrium.
template <class IdT, class T>
class Attribute
{
  static_assert( ! std::is_same_v<T, bool>,
                 "gts: std::vector<bool> packs bits and cannot hand out bool&; use std::uint8_t");
public:
  Attribute() = default;
  explicit Attribute(std::vector<T>* v) : m_v(v) {}

  T&       operator[](IdT id)       { return (*m_v)[id.i]; }
  const T& operator[](IdT id) const { return (*m_v)[id.i]; }

private:
  std::vector<T>* m_v = nullptr;
};

namespace detail {

// The type-erased storage of one attribute, kept by the Geometrium.
struct Column
{
  std::string           name;
  std::type_index       type = typeid(void);
  std::shared_ptr<void> data; // a std::vector<T>

  std::function<void(std::size_t)>                            resize;
  std::function<void(const std::vector<Index>&, std::size_t)> renumber;    // old -> new, kNone dropped
  std::function<void(Index dst, Index a, Index b, double t)>  interpolate;
  std::function<void(Index dst, Index a, Index b)>            merge;
  std::function<void(Index dst, Index parent)>                split;       // dst may be parent
  std::function<void(Index dst, const Column& src, Index s)>  copy_from;   // from another Geometrium

  template <class T>
  static Column make(std::string name, std::size_t n, Hooks<T> hooks)
  {
    static_assert( ! std::is_same_v<T, bool>,
                   "gts: std::vector<bool> packs bits and cannot hand out bool&; use std::uint8_t");
    auto data = std::make_shared<std::vector<T>>(n);
    std::vector<T>* v = data.get();
    Column c;
    c.name   = std::move(name);
    c.type   = typeid(T);
    c.data   = data;
    c.resize = [v](std::size_t m) { v->resize(m); };
    c.renumber = [v](const std::vector<Index>& map, std::size_t m) {
      std::vector<T> w(m);
      for (std::size_t i = 0; i < map.size(); ++i)
        if (map[i] != kNone) w[map[i]] = std::move((*v)[i]);
      v->swap(w);
    };
    using D = DefaultHooks<T>;
    if (hooks.interpolate)
      c.interpolate = [v, h = std::move(hooks.interpolate)](Index d, Index a, Index b, double t) { (*v)[d] = h((*v)[a], (*v)[b], t); };
    else
      c.interpolate = [v](Index d, Index a, Index b, double t) { (*v)[d] = D::interpolate((*v)[a], (*v)[b], t); };
    if (hooks.merge)
      c.merge = [v, h = std::move(hooks.merge)](Index d, Index a, Index b) { (*v)[d] = h((*v)[a], (*v)[b]); };
    else
      c.merge = [v](Index d, Index a, Index b) { (*v)[d] = D::merge((*v)[a], (*v)[b]); };
    if (hooks.split)
      c.split = [v, h = std::move(hooks.split)](Index d, Index s) { (*v)[d] = h((*v)[s]); };
    else
      c.split = [v](Index d, Index s) { (*v)[d] = D::split((*v)[s]); };
    c.copy_from = [v](Index d, const Column& src, Index s) {
      (*v)[d] = (*static_cast<const std::vector<T>*>(src.data.get()))[s];
    };
    return c;
  }

  template <class T>
  std::vector<T>* vector() const { return static_cast<std::vector<T>*>(data.get()); }
};

} // namespace detail
} // namespace gts
