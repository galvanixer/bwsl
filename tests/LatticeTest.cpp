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
#include <LatticeFactory.hpp>
#include <Lattice.hpp>
#include <SquareLattice.hpp>

// std
#include <vector>
#include <cmath>

// catch
#define CATCH_CONFIG_MAIN
#include <catch/catch.hpp>

using namespace bwsl;

// TDD example
TEST_CASE("Square Lattice", "[lattice]")
{

  Lattice* structure = Lattice::CreateLattice("square", {3,3});

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
    REQUIRE(structure->GetDistance(0,4) == Approx(sqrt(2)));
    REQUIRE(structure->GetDistance(0,8) == Approx(sqrt(2)));
    REQUIRE(structure->GetDistance(4,7) == Approx(1.0));
  }
}

// vim: set ft=cpp ts=2 sts=2 et sw=2 tw=80: //