//===-- ObservableGroup.hpp ------------------------------------*- C++ -*-===//
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
/// @brief      Definitions for the ObservableGroup Class
///
//===---------------------------------------------------------------------===//
#pragma once

// bwsl
#include <bwsl/Accumulators.hpp>

// fmt
#include <fmt/format.h>
#include <fmt/ostream.h>

// boost
#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/version.hpp>

// std
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace bwsl {

/// Description
template<typename Index_t>
class ObservableGroup
{
public:
  /// Default constructor
  ObservableGroup() = default;

  /// Copy constructor
  ObservableGroup(const ObservableGroup&) = default;

  /// Copy assignment operator
  auto operator=(const ObservableGroup&) -> ObservableGroup& = default;

  /// Move constructor
  ObservableGroup(ObservableGroup&&) noexcept = default;

  /// Move assignment operator
  auto operator=(ObservableGroup&&) noexcept -> ObservableGroup& = default;

  /// Default destructor
  virtual ~ObservableGroup() = default;

  /// Constructor
  ObservableGroup(std::string output_file);

  /// Constructor
  ObservableGroup(std::string output_file, std::vector<Index_t> const& indices);

  /// Measure an observable
  void Measure(Index_t idx, double val);

  /// Print the headers to the file
  void PrintHeaders() const;

  /// Print the results and reset all the accumulators
  void PrintAndReset(size_t precision = 10UL);

  auto AddObservable(Index_t key) -> ObservableGroup&;

protected:
private:
  /// Name of the associated output file
  std::string output_file_;

  /// Storage for the accumulators
  std::map<Index_t, accumulators::KahanAccumulator> accumulator_{};

  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive& ar, unsigned int version);
}; // class ObservableGroup

template<typename Index_t>
ObservableGroup<Index_t>::ObservableGroup(std::string output_file)
  : output_file_(std::move(output_file))
{
}

template<typename Index_t>
ObservableGroup<Index_t>::ObservableGroup(std::string output_file,
                                          std::vector<Index_t> const& indices)
  : output_file_(std::move(output_file))
{
  for (auto i : indices) {
    accumulator_.try_emplace(i);
  }
}

template<typename Index_t>
inline void
ObservableGroup<Index_t>::Measure(Index_t idx, double val)
{
  accumulator_.at(idx).Add(val);
}

template<typename Index_t>
inline void
ObservableGroup<Index_t>::PrintHeaders() const
{
  auto out = std::ofstream(output_file_.c_str(), std::ios::trunc);

  auto it = accumulator_.begin();
  fmt::print(out, "{}", it->first);
  while (++it != accumulator_.end()) {
    fmt::print(out, ",{}", it->first);
  }
  fmt::print(out, ",Count");
  fmt::print(out, "\n");
}

template<typename Index_t>
inline void
ObservableGroup<Index_t>::PrintAndReset(size_t precision)
{
  auto out = std::ofstream(output_file_.c_str(), std::ios::app);
  auto it = accumulator_.begin();
  auto count = it->second.Count();
  fmt::print(out, "{:.{}e}", it->second.Mean(), precision);
  it->second.Reset();
  while (++it != accumulator_.end()) {
    fmt::print(out, ",{:.{}e}", it->second.Mean(), precision);
    it->second.Reset();
  }
  fmt::print(out, ",{}", count);
  fmt::print(out, "\n");
}

template<typename Index_t>
inline auto
ObservableGroup<Index_t>::AddObservable(Index_t key)
  -> ObservableGroup<Index_t>&
{
  accumulator_.try_emplace(key);
  return *this;
}

template<typename Index_t>
template<class Archive>
void
ObservableGroup<Index_t>::serialize(Archive& ar,
                                    const unsigned int /* version */)
{
  // clang-format off
  ar & output_file_;
  ar & accumulator_;
  // clang-format on
}

} // namespace bwsl

// vim: set ft=cpp ts=2 sts=2 et sw=2 tw=80: //
