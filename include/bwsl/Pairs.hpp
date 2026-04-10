//===-- Pairs.hpp ----------------------------------------------*- C++ -*-===//
//
//                       BeagleWarlord's Support Library
//
// Copyright 2016-2022 Guido Masella. All Rights Reserved.
// See LICENSE file for details
//
//===---------------------------------------------------------------------===//
///
/// @file
/// @author     Guido Masella (guido.masella@gmail.com)
/// @brief      Definitions for the Pairs
///
//===---------------------------------------------------------------------===//
#pragma once

// bwsl
#include <bwsl/MathUtils.hpp>

// std
#include <cmath>
#include <cstdlib>
#include <utility>

namespace bwsl::pairs {

///
/// Get the index of a pair from the values of the two single indices
///
inline auto
GetPairIndex(std::size_t a, std::size_t b, std::size_t maxval) -> std::size_t
{
  auto i = b * maxval + a;
  return i;
}

///
/// Get the index of a triple from the values of the three single indices
///
inline auto
GetTripleIndex(std::size_t a, std::size_t b, std::size_t c, std::size_t maxval) -> std::size_t
{
  auto i = c * maxval * maxval + b * maxval + a;
  return i;
}

///
/// Get the index of the pair from the pair of indices
///
inline auto
GetPairIndex(const std::pair<std::size_t, std::size_t> pair,
             const std::size_t maxval) -> std::size_t
{
  return GetPairIndex(pair.first, pair.second, maxval);
}


///
/// Get a triple of indices from the index of the triple
///
inline auto
GetTriple(std::size_t i, const std::size_t maxval)
  -> std::tuple<std::size_t, std::size_t, std::size_t>
{
  return std::make_tuple(i % maxval,
                         (i / maxval) % maxval,
                         i / (maxval * maxval));
}
///
/// Get a pair of indices from the index of the pair
///
inline auto
GetPair(std::size_t i, const std::size_t maxval)
  -> std::pair<std::size_t, std::size_t>
{
  return std::make_pair(i % maxval, i / maxval);
}

///
/// Get the number of indices of pairs of unsigned values
///
inline auto
GetNumPairs(std::size_t maxval) -> std::size_t
{
  return square(maxval);
}

///
/// Get the number of indices of triples of unsigned values
///
inline auto
GetNumTriples(std::size_t maxval) -> std::size_t
{
  return cube(maxval);
}


} // namespace bwsl::pairs

// vim: set ft=cpp ts=2 sts=2 et sw=2 tw=80: //
