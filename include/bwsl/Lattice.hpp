//===-- Lattice.hpp --------------------------------------------*- C++ -*-===//
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
/// @brief      Definitions for the Lattice Class
///
//===---------------------------------------------------------------------===//
#pragma once

// bwsl
#include <bwsl/Approx.hpp>
#include <bwsl/Bravais.hpp>
#include <bwsl/HyperCubicGrid.hpp>
#include <bwsl/MathUtils.hpp>
#include <bwsl/Pairs.hpp>

// fmt
#include <fmt/format.h>
#include <fmt/ostream.h>

// std
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace bwsl {

///
/// Representation of a Lattice.
///
class Lattice : public HyperCubicGrid
{
public:
  /// Vector of offsets
  using vectorindex_t = std::vector<index_t>;

  /// Vector of neighbors
  using neighbors_t = std::vector<vectorindex_t>;

  /// Shorthand for real valued vectors
  using realvec_t = std::vector<double>;

  /// Default constructor
  Lattice() = default;

  /// Construct a lattice with given size from an infinite bravais lattice
  Lattice(Bravais const& bravais,
          gridsize_t const& size,
          boundaries_t boundaries = boundaries_t::Closed);

  /// Copy constructor
  Lattice(Lattice const& that) = default;

  /// Move constructor
  Lattice(Lattice&& that) = default;

  /// Copy assignment operator
  auto operator=(Lattice const& that) -> Lattice& = default;

  /// Move assignment operator
  auto operator=(Lattice&& that) -> Lattice& = default;

  /// Default destructor
  ~Lattice() override = default;

  /// Get nearest neighbors of site i
  [[nodiscard]] auto GetNeighbors(index_t i) const -> vectorindex_t const&
  {
    return neighbors_[i];
  }

  /// Distance betweeen two sites of the lattice
  [[nodiscard]] auto GetDistance(index_t a, index_t b) const -> double;

  /// Get the vector spawning from site @p a to site @p b .
  [[nodiscard]] auto GetVector(index_t a, index_t b) const -> realvec_t;

  /// Given an accumulator for the winding number return the total winding.
  /// The total winding is defined as the number of times we jump around the
  /// boundaries.
  [[nodiscard]] auto GetWinding(coords_t jump) const -> coords_t;

  /// Get the real space coordinates of site @p a .
  [[nodiscard]] auto GetPosition(index_t a) const -> realvec_t;

  /// Change a vector so that he will match the boundary conditions
  auto EnforceBoundaries(coords_t& coords) const -> void;

  /// Check if two sites are neighbors
  [[nodiscard]] auto AreNeighbors(index_t a, index_t b) const -> bool;

  /// Get the coordination number
  [[nodiscard]] auto GetCoordination(index_t a) const -> index_t;

  /// Get the coordination number
  [[nodiscard]] auto GetCoordination() const -> index_t;

  /// Get an allowed momentum
  [[nodiscard]] auto GetMomentum(size_t a) const -> realvec_t;

  /// Compute the structure factor given the occupations of the sites
  template<class T>
  auto AccumulateSk(std::vector<T> const& occupations,
                    realvec_t& sk,
                    double mult = 1.0) const -> void;

  /// Compute the structure factor given the occupations of the sites
  template<class T>
  [[nodiscard]] auto ComputeSk(std::vector<T> const& occupations,
                               double mult = 1.0) const -> realvec_t;

  /// Save the distances on a file
  auto SaveDistances(const std::string& fname) const -> void;

  /// Save the positions on a file
  auto SavePositions(const std::string& fname) const -> void;

  /// Save the momenta on a file
  auto SaveMomenta(const std::string& fname) const -> void;

  /// Save the distances on a file
  auto SavePairs(const std::string& fname) const -> void;

  /// Save the distances on a file
  auto SaveTriples(const std::string& fname) const -> void;

protected:
  /// Compute the positions of all the lattice points
  /// NOTE: the positions stored are in real space.
  [[nodiscard]] auto ComputePositions(Bravais const& bravais) const
    -> std::vector<realvec_t>;

  /// Compute the distance vectors.
  /// It is composed of vectors in real space between site 0 and site i.
  /// Here, i ∈ [0,1,...,N] where N is the largest site index.
  /// For closed boundary conditions each vector is chosen in such a manner
  /// that it represents the shortest distance between the first site and all
  /// the periodic images of the second one (or vice versa).
  /// [minimum distance convention]
  [[nodiscard]] auto ComputeVectors(Bravais const& bravais) const
    -> std::vector<realvec_t>;

  /// Compute the vector of distances respecting the minimum
  /// distance convention (if with closed boundaries).
  /// It is composed of magnitudes of distance vectors computed using
  /// ComputeVectors() above.
  [[nodiscard]] auto ComputeDistances(Bravais const& /*bravais*/) const
    -> realvec_t;

  /// Create the vector storing vector of neighbors for each lattice site.
  /// The vector of neighbors store the site indices of neighbors.
  /// It also takes into account the boundary conditions.
  [[nodiscard]] auto ComputeNeighbors(Bravais const& bravais) const
    -> neighbors_t;

  /// Compute the allowed momenta
  [[nodiscard]] auto ComputeMomenta(Bravais const& bravais) const
    -> std::vector<realvec_t>;

private:
  /// Positions of all the sites
  /// Assuming that the first site has position `(0,0)`
  std::vector<realvec_t> position_{};

  /// All the distance vectors between pairs of sites
  std::vector<realvec_t> vectors_{};

  /// All the distances on the lattice with minimum image convention
  realvec_t distance_{};

  /// vector of nearest neighbors
  neighbors_t neighbors_{};

  /// Allowed values momenta
  std::vector<realvec_t> momenta_{};
}; // class Lattice

inline Lattice::Lattice(Bravais const& bravais,
                        gridsize_t const& size,
                        boundaries_t boundaries)
  : HyperCubicGrid(size, boundaries)
  , position_(ComputePositions(bravais))
  , vectors_(ComputeVectors(bravais))
  , distance_(ComputeDistances(bravais))
  , neighbors_(ComputeNeighbors(bravais))
  , momenta_(ComputeMomenta(bravais))
{
}

inline void
Lattice::EnforceBoundaries(coords_t& coords) const
{
  for (auto i = 0UL; i < GetDim(); i++) {
    auto s = GetSize()[i];
    while (coords[i] < 0) {
      coords[i] += s;
    }
    while (coords[i] >= static_cast<long>(s)) {
      coords[i] -= s;
    }
  }
}

inline auto
Lattice::GetDistance(size_t a, size_t b) const -> double
{
  assert(IndexIsValid(a) && IndexIsValid(b));
  auto s = GetMappedSite(a, b);
  return distance_[s];
}

inline auto
Lattice::GetVector(size_t a, size_t b) const -> realvec_t
{
  assert(IndexIsValid(a) && IndexIsValid(b));
  auto s = GetMappedSite(a, b);
  return vectors_[s];
}

inline auto
Lattice::GetMomentum(size_t a) const -> realvec_t
{
  assert(IndexIsValid(a));
  return momenta_[a];
}

inline auto
Lattice::GetWinding(coords_t jumps) const -> coords_t
{
  assert(HasSameDimension(jumps));

  std::transform(jumps.begin(),
                 jumps.end(),
                 GetSize().cbegin(),
                 jumps.begin(),
                 std::divides<typename coords_t::value_type>());

  return jumps;
}

inline auto
Lattice::GetPosition(index_t a) const -> Lattice::realvec_t
{
  assert(IndexIsValid(a));
  return position_[a];
}

inline auto
Lattice::AreNeighbors(index_t a, index_t b) const -> bool
{
  assert(a < neighbors_.size());

  auto result = std::find(neighbors_[a].begin(), neighbors_[a].end(), b);

  return result != neighbors_[a].end();
}

inline auto
Lattice::ComputePositions(Bravais const& bravais) const
  -> std::vector<Lattice::realvec_t>
{
  auto p = std::vector<realvec_t>{};
  auto c0 = GetCoordinates(0);
  for (auto i = 0UL; i < GetNumSites(); i++) {
    auto x = bravais.GetVector(c0, GetCoordinates(i));
    p.push_back(x);
  }

  return p;
}

inline auto
Lattice::ComputeDistances(Bravais const& /*bravais*/) const
  -> Lattice::realvec_t
{
  auto p = realvec_t{};
  std::transform(vectors_.begin(),
                 vectors_.end(),
                 std::back_inserter(p),
                 [](realvec_t const& x) -> double {
                   return sqrt(sum_squared<realvec_t, double>(x));
                 });
  return p;
}

inline auto
Lattice::ComputeVectors(Bravais const& bravais) const
  -> std::vector<Lattice::realvec_t>
{
  auto p = std::vector<realvec_t>(GetNumSites(), realvec_t(GetDim(), 0.0));
  const auto imgsize = gridsize_t(GetDim(), 3);
  const auto nimg = accumulate_product(imgsize);

  const auto c0 = GetCoordinates(0UL);
  for (auto site = 0UL; site < GetNumSites(); site++) {
    auto cs = GetCoordinates(site);

    // minimum distance found
    auto [mindist, minvec] = bravais.GetDistanceVector(c0, cs);

    // search for the minimum distance across the first shell of
    // periodic images if we are with closed boundary condition
    if (HasClosedBoundaries()) {
      for (auto k = 0UL; k < nimg; k++) {
        auto img = index_to_array<coords_t, gridsize_t>(k, imgsize);
        auto csm = coords_t(cs);
        for (auto m = 0UL; m < GetDim(); m++) {
          csm[m] += (img[m] - 1) * GetSize()[m];
        }
        auto [dist, vec] = bravais.GetDistanceVector(c0, csm);

        if (dist < mindist) {
          mindist = dist;
          minvec = vec;
        }
      }
    }
    p[site] = minvec;
  }
  return p;
}

inline auto
Lattice::ComputeNeighbors(Bravais const& bravais) const -> Lattice::neighbors_t
{
  auto p = neighbors_t{};
  const auto gamma = bravais.GetGamma();

  for (auto i = 0UL; i < GetNumSites(); i++) {
    auto nn = gridsize_t{};
    auto ci = GetCoordinates(i);
    for (auto j = 0UL; j < gamma; j++) {
      auto cj = bravais.GetNeighbor(ci, j);

      // add the sites to the list of neighbors accordingly with the
      // boundary conditions set
      if (HasClosedBoundaries() || IsOnGrid(cj)) {
        EnforceBoundaries(cj);
        nn.push_back(GetIndex(cj));
      }
    }
    p.push_back(std::move(nn));
  }

  return p;
}

inline auto
Lattice::ComputeMomenta(Bravais const& bravais) const
  -> std::vector<Lattice::realvec_t>
{
  auto p = std::vector<realvec_t>{};

  // with open boundary conditions the momenta are not defined
  if (HasOpenBoundaries()) {
    return p;
  }

  // in a lattice with periodic boundary conditions the number
  // of reciprocal lattice vectors are the same number as the
  // sites of the direct lattice
  for (auto i = 0UL; i < GetNumSites(); i++) {
    auto ci = GetCoordinates(i);
    for (auto k = 0UL; k < ci.size(); k++) {
      ci[k] -= GetSize()[k] / 2UL;
    }
    auto kappa = bravais.GetReciprocalSpace(ci);
    for (auto k = 0UL; k < kappa.size(); k++) {
      kappa[k] /= GetSize()[k];
    }
    p.push_back(kappa);
  }

  return p;
}

template<class T>
inline auto
Lattice::AccumulateSk(std::vector<T> const& occupations,
                      realvec_t& sk,
                      double mult) const -> void
{
  auto n = GetNumSites();

  if (HasOpenBoundaries()) {
    return;
  }

  for (auto i = 0UL; i < GetNumSites(); i++) {
    auto const k = momenta_[i];
    auto im = 0.0;
    auto re = 0.0;

    for (auto j = 0UL; j < n; j++) {
      auto const x = GetVector(0, j);
      auto prod = 0.0;
      for (auto q = 0UL; q < GetDim(); q++) {
        prod += k[q] * x[q];
      }
      im += sin(prod) * occupations[j];
      re += cos(prod) * occupations[j];
    }
    sk[i] += mult * (square(im) + square(re)) / static_cast<double>(square(n));
  }
}

template<class T>
inline auto
Lattice::ComputeSk(std::vector<T> const& occupations, double mult) const
  -> realvec_t
{
  auto n = GetNumSites();
  auto sk = std::vector<double>(n, 0.0);

  if (!HasOpenBoundaries()) {
    AccumulateSk(occupations, sk, mult);
  }

  return sk;
}

inline auto
Lattice::SavePositions(std::string const& fname) const -> void
{
  auto out = std::ofstream{ fname.c_str() };

  fmt::print(out, "i");
  for (auto i = 0UL; i < GetDim(); i++) {
    fmt::print(out, ",x{}", i);
  }
  fmt::print(out, "\n");

  for (auto i = 0UL; i < GetNumSites(); i++) {
    auto vec = GetPosition(i);
    fmt::print(out, "{}", i);
    for (auto const& v : vec) {
      fmt::print(out, ",{}", v);
    }
    fmt::print(out, "\n");
  }
}

inline auto
Lattice::SaveDistances(std::string const& fname) const -> void
{
  auto out = std::ofstream{ fname.c_str() };

  fmt::print(out, "i");
  for (auto i = 0UL; i < GetDim(); i++) {
    fmt::print(out, ",d{}", i);
  }
  fmt::print(out, "\n");

  for (auto i = 0UL; i < GetNumSites(); i++) {
    auto vec = GetVector(0, i);
    fmt::print(out, "{}", i);
    for (auto const& v : vec) {
      fmt::print(out, ",{}", v);
    }
    fmt::print(out, "\n");
  }
}

inline void
Lattice::SaveMomenta(std::string const& fname) const
{
  auto out = std::ofstream{ fname.c_str() };

  fmt::print(out, "i");
  for (auto i = 0UL; i < GetDim(); i++) {
    fmt::print(out, ",k{}", i);
  }
  fmt::print(out, "\n");

  for (auto i = 0UL; i < GetNumSites(); i++) {
    fmt::print(out, "{}", i);
    for (auto const& k : momenta_[i]) {
      fmt::print(out, ",{}", k);
    }
    fmt::print(out, "\n");
  }
}

inline auto
Lattice::SavePairs(std::string const& fname) const -> void
{
  auto out = std::ofstream{ fname.c_str() };

  fmt::print(out, "i,a,b");

  auto npairs = pairs::GetNumPairs(GetNumSites());
  auto nsites = GetNumSites();

  for (auto i = 0UL; i < GetDim(); i++) {
    fmt::print(out, ",x{}", i);
  }
  for (auto i = 0UL; i < GetDim(); i++) {
    fmt::print(out, ",y{}", i);
  }
  for (auto i = 0UL; i < GetDim(); i++) {
    fmt::print(out, ",d{}", i);
  }
  fmt::print(out, "\n");

  for (auto i = 0UL; i < npairs; i++) {
    auto [a, b] = bwsl::pairs::GetPair(i, nsites);

    auto va = GetVector(0, a);
    fmt::print(out, "{},{},{}", bwsl::pairs::GetPairIndex(a, b, nsites), a, b);

    for (auto const& v : va) {
      fmt::print(out, ",{}", v);
    }

    auto vb = GetVector(0, b);
    for (auto const& v : vb) {
      fmt::print(out, ",{}", v);
    }

    auto vec = GetVector(a, b);
    for (auto const& v : vec) {
      fmt::print(out, ",{}", v);
    }

    fmt::print(out, "\n");
  }
}

inline auto
Lattice::SaveTriples(std::string const& fname) const -> void
{
    auto out = std::ofstream{ fname.c_str() };

    fmt::print(out, "i,a,b,c");

    auto nsites   = GetNumSites();
    auto ntriples = nsites * nsites * nsites;  // MUST match indexing scheme

    // Coordinates for each site
    for (auto i = 0UL; i < GetDim(); i++) fmt::print(out, ",ax{}", i);
    for (auto i = 0UL; i < GetDim(); i++) fmt::print(out, ",bx{}", i);
    for (auto i = 0UL; i < GetDim(); i++) fmt::print(out, ",cx{}", i);

    // Relative vectors a→b, a→c, b→c
    for (auto i = 0UL; i < GetDim(); i++) fmt::print(out, ",d_ab_{}", i);
    for (auto i = 0UL; i < GetDim(); i++) fmt::print(out, ",d_ac_{}", i);
    for (auto i = 0UL; i < GetDim(); i++) fmt::print(out, ",d_bc_{}", i);

    fmt::print(out, "\n");

    for (std::size_t idx = 0; idx < ntriples; idx++)
    {
        // Matches EXACTLY your GetTriple definition
        auto [a, b, c] = bwsl::pairs::GetTriple(idx, nsites);

        fmt::print(out, "{},{},{},{}", idx, a, b, c);

        // Coordinates at sites a, b, c
        auto va = GetVector(0, a);
        auto vb = GetVector(0, b);
        auto vc = GetVector(0, c);

        for (auto const& v : va) fmt::print(out, ",{}", v);
        for (auto const& v : vb) fmt::print(out, ",{}", v);
        for (auto const& v : vc) fmt::print(out, ",{}", v);

        // Relative vectors
        auto vab = GetVector(a, b);
        auto vac = GetVector(a, c);
        auto vbc = GetVector(b, c);

        for (auto const& v : vab) fmt::print(out, ",{}", v);
        for (auto const& v : vac) fmt::print(out, ",{}", v);
        for (auto const& v : vbc) fmt::print(out, ",{}", v);

        fmt::print(out, "\n");
    }
}

inline auto
Lattice::GetCoordination(size_t a) const -> size_t
{
  return neighbors_[a].size();
}

inline auto
Lattice::GetCoordination() const -> size_t
{
  return neighbors_[0].size();
}

} // namespace bwsl

// vim: set ft=cpp ts=2 sts=2 et sw=2 tw=80: //
