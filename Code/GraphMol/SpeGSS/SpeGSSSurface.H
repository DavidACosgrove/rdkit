//
// Copyright (C) David Cosgrove 2026
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This is the interface for the class SpeGSSSurface which holds
// a surface for use with SpeGGS.

#ifndef RDKIT_SPEGSSSURFACE_H
#define RDKIT_SPEGSSSURFACE_H

#include <RDGeneral/export.h>
#include <GraphMol/SpeGSS/SpeGSSParams.h>

#include <cstdint>
#include <unordered_set>

#include <boost/container_hash/hash.hpp>

namespace RDKit {
class ROMol;

namespace SpeGSS {

using Edge = std::pair<std::uint64_t, std::uint64_t>;

class SpeGSSSurface {
 public:
  SpeGSSSurface() = default;
  SpeGSSSurface(const ROMol &mol, const SpeGSSParams &params = SpeGSSParams());
  // Create a surface just with the given subset of vertices.
  SpeGSSSurface(const SpeGSSSurface &surf,
                const std::vector<std::uint64_t> &vtxNums);
  SpeGSSSurface(const SpeGSSSurface &surf) = default;
  SpeGSSSurface(SpeGSSSurface &&surf) = default;
  SpeGSSSurface &operator=(const SpeGSSSurface &surf) = default;
  SpeGSSSurface &operator=(SpeGSSSurface &&surf) = default;
  ~SpeGSSSurface() = default;

  std::uint64_t getNumVertices() const { return d_vertices.size() / 3; }
  std::uint64_t getNumTriangles() const { return d_triangles.size() / 3; }
  std::uint64_t getNumEdges() const { return d_edges.size(); }
  double getVolume() const { return d_volume; }
  const std::vector<double> &getVertices() const { return d_vertices; }
  const std::vector<double> &getNormals() const { return d_normals; }
  const std::vector<std::uint64_t> &getTriangles() const { return d_triangles; }
  const std::unordered_set<Edge, boost::hash<Edge>> &getEdges() const {
    return d_edges;
  }

  // Returns true if the surface is a manifold.
  bool isManifold() const;

  void writeFile(const std::string &filename) const;
  void writeToStream(std::ostream &os) const;
  bool readFile(const std::string &filename);
  bool readFromStream(std::istream &is);

 private:
  // 3 doubles per vertex, x, y, z coords
  std::vector<double> d_vertices;
  // 3 doubles per normal, x, y and z components. One normal per
  // vertex, calculated by the mc33.
  std::vector<double> d_normals;
  // 3 ints per triangle, as supplied by mc33.
  std::vector<std::uint64_t> d_triangles;
  // d_edges are the edges that make up the triangles, 2 ints per edge.
  // Each edge is is stored only once, with large vertex number first.
  std::unordered_set<Edge, boost::hash<Edge>> d_edges;
  // The connection table of vertices.
  std::vector<std::vector<std::uint64_t>> d_vertConns;

  // Calculated from the grid, by counting the number of grid
  // points inside the contour value and multiplying by the grid
  // cube size.
  double d_volume{0.0};

  // Create the edges etc.
  void finishBuildingSurface();

  // Returns true if a vertex has < 2 connections, so is either
  // completely isolated or only in 1 edge.
  bool strayPointsAndEdges() const;
  // To be a manifold, all edges must be in no more than 2 triangles.
  // For a molecule, the surface should be closed so there should
  // be no edges in only 1 triangle, but that's not checked at
  // the moment.
  bool edgeCountsGoodForManifold() const;
  // Returns true if there's more than 1 piece to the surface. I.e.,
  // if starting from vertex 0 and moving from connected vertex to
  // connected vertex, all points aren't visited.
  bool oneComponentInMesh() const;
  // The triangles off each vertex must form a fan.  Returns true if so.
  // A fan is defined as having contiguous triangles, i.e. starting at one
  // of the neighbours of the vertex, and moving round in ascending angle,
  // each neighbour is connected to the next.
  bool vertexTrianglesAreFan() const;
};

}  // namespace SpeGSS
}  // namespace RDKit
#endif  // RDKIT_SPEGSSSURFACE_H
