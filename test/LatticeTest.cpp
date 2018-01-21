//===-- LatticeTest.cpp ----------------------------------------*- C++ -*-===//
//
//                       BeagleWarlord's Support Library
//
// Copyright 2017 Guido Masella. All Rights Reserved.
// See LICENSE file for details
//
//===---------------------------------------------------------------------===//
///
/// \file
/// \version    0.1
/// \author     Guido Masella (guido.masella@gmail.com)
/// \date       May, 2017
/// \brief      Tests for the Lattice Class
///
//===---------------------------------------------------------------------===//
// bwsl
#include <bwsl/LatticeFactory.hpp>
#include <bwsl/Lattice.hpp>
#include <bwsl/SquareLattice.hpp>

// std
#include <vector>
#include <cmath>
#include <random>

// catch
#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>

using namespace bwsl;

// TDD example
TEST_CASE("Square Lattice", "[lattice]")
{

  auto structure = Lattice::CreateLattice("square", {3,3});

  SECTION("neighbors are correct")
  {
    REQUIRE(structure->AreNeighbors(4,1));
    REQUIRE(structure->AreNeighbors(1,4));
    REQUIRE(structure->AreNeighbors(0,1));
    REQUIRE(structure->AreNeighbors(0,3));
    REQUIRE(structure->AreNeighbors(0,6));
    REQUIRE(structure->AreNeighbors(0,2));
    REQUIRE(structure->AreNeighbors(8,5));
    REQUIRE(structure->AreNeighbors(8,7));
    REQUIRE(structure->AreNeighbors(8,2));
    REQUIRE(structure->AreNeighbors(8,6));
  }
  SECTION("distances are correct")
  {
    REQUIRE(structure->GetDistance(0,1) == Approx(1.0));
    REQUIRE(structure->GetDistance(0,2) == Approx(1.0));
    REQUIRE(structure->GetDistance(0,4) == Approx(sqrt(2)));
    REQUIRE(structure->GetDistance(0,8) == Approx(sqrt(2)));
    REQUIRE(structure->GetDistance(4,7) == Approx(1.0));
  }
  SECTION("windings are correct")
  {
    auto w = structure->GetWinding(0,1);
    auto w1 = structure->GetWinding(0,2);
    auto w3 = structure->GetWinding(2,8);
    REQUIRE(w.size() == 2);
    REQUIRE(w[0] == 0);
    REQUIRE(w[1] == 1);
    REQUIRE(w1[0] == 0);
    REQUIRE(w1[1] == -1);
    REQUIRE(w3[0] == -1);
    REQUIRE(w3[1] == 0);
  }
  SECTION("GetCoordinates and GetOffset invert each other")
  {
    for (auto i = 0ul; i < structure->GetNumSites(); i++) {
      auto offs = structure->GetOffset(structure->GetCoordinates(i));
      REQUIRE(offs == i);
    }
  }
  SECTION("GetRandomdistance is ok")
  {
    auto rng = std::mt19937_64{};
    auto dist = structure->GetRandomDistance(rng);
    REQUIRE(dist.size() == 2);
    REQUIRE(dist[0] >= -3);
    REQUIRE(dist[0] <= 3);
  }
}

TEST_CASE("Big Square Lattice")
{
  auto structure = Lattice::CreateLattice("square", {36,36});
  auto nsites = structure->GetNumSites();
  auto size = structure->GetSize();

  SECTION("distances are correct")
  {
    for (auto i = 0; i < nsites; i++) {
      auto coordi = structure->GetCoordinates(i);
      for (auto j = 0; j< nsites; j++) {
        if (i == j) continue;
        auto coordj = structure->GetCoordinates(j);
        auto distx = coordi[0] - coordj[0];
        auto disty = coordi[1] - coordj[1];
        if (distx > static_cast<long>(size[0]/2))
          distx -= static_cast<long>(size[0]);
        else if (distx < -static_cast<long>(size[0]/2))
          distx += static_cast<long>(size[0]);
        if (disty > static_cast<long>(size[1]/2))
          disty -= static_cast<long>(size[1]);
        else if (disty < -static_cast<long>(size[1]/2))
          disty += static_cast<long>(size[1]);
        auto dist = sqrt(static_cast<double>(bwsl::square(distx) + bwsl::square(disty)));

        REQUIRE(structure->GetDistance(i,j,0) == Approx(distx));
        REQUIRE(structure->GetDistance(i,j,1) == Approx(disty));
        REQUIRE(structure->GetDistance(i,j) == Approx(dist));
      }
    }
  }


}

// vim: set ft=cpp ts=2 sts=2 et sw=2 tw=80: //