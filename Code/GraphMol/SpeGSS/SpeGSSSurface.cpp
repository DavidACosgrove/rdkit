//
// Copyright (C) David Cosgrove 2026
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <algorithm>
#include <fstream>
#include <list>
#include <ranges>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "SpeGSSParams.h"

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SpeGSS/SpeGSSSurface.H>

#include <boost/dynamic_bitset.hpp>

// This is the Marching Cubes 33 implementation for generating the
// surface from the Gaussian density of a molecule.
#define mc33cpp_implementation
#include <MC33/mc33/mc33cpp.h>

namespace RDKit {
namespace SpeGSS {

// Bondi radii
// You can find more of these in Table 12 of this publication:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
const std::map<unsigned int, double> vdw_radii = {
    {0, 2.16},   // Dummy, same as Xe.
    {1, 1.10},   // H
    {2, 1.40},   // He
    {3, 1.81},   // Li
    {4, 1.53},   // Be
    {5, 1.92},   // B
    {6, 1.70},   // C
    {7, 1.55},   // N
    {8, 1.52},   // O
    {9, 1.47},   // F
    {10, 1.54},  // Ne
    {11, 2.27},  // Na
    {12, 1.73},  // Mg
    {13, 1.84},  // Al
    {14, 2.10},  // Si
    {15, 1.80},  // P
    {16, 1.80},  // S
    {17, 1.75},  // Cl
    {18, 1.88},  // Ar
    {19, 2.75},  // K
    {20, 2.31},  // Ca
    {31, 1.87},  // Ga
    {32, 2.11},  // Ge
    {33, 1.85},  // As
    {34, 1.90},  // Se
    {35, 1.83},  // Br
    {36, 2.02},  // Kr
    {37, 3.03},  // Rb
    {38, 2.49},  // Sr
    {49, 1.93},  // In
    {50, 2.17},  // Sn
    {51, 2.06},  // Sb
    {52, 2.06},  // Te
    {53, 1.98},  // I
    {54, 2.16},  // Xe
    {55, 3.43},  // Cs
    {56, 2.68},  // Ba
    {81, 1.96},  // Tl
    {82, 2.02},  // Pb
    {83, 2.07},  // Bi
    {84, 1.97},  // Po
    {85, 2.02},  // At
    {86, 2.20},  // Rn
    {87, 3.48},  // Fr
    {88, 2.83},  // Ra
};

namespace {
void getMoleculeExtremes(const ROMol &mol, int confId, double &xmin,
                         double &xmax, double &ymin, double &ymax, double &zmin,
                         double &zmax, double molPadding) {
  xmin = std::numeric_limits<double>::max();
  xmax = std::numeric_limits<double>::lowest();
  ymin = std::numeric_limits<double>::max();
  ymax = std::numeric_limits<double>::lowest();
  zmin = std::numeric_limits<double>::max();
  zmax = std::numeric_limits<double>::lowest();
  const auto &conf = mol.getConformer(confId);
  for (const auto &atom : mol.atoms()) {
    const auto &pos = conf.getAtomPos(atom->getIdx());
    if (pos.x < xmin) {
      xmin = pos.x;
    }
    if (pos.x > xmax) {
      xmax = pos.x;
    }
    if (pos.y < ymin) {
      ymin = pos.y;
    }
    if (pos.y > ymax) {
      ymax = pos.y;
    }
    if (pos.z < zmin) {
      zmin = pos.z;
    }
    if (pos.z > zmax) {
      zmax = pos.z;
    }
  }
  xmin -= molPadding;
  xmax += molPadding;
  ymin -= molPadding;
  ymax += molPadding;
  zmin -= molPadding;
  zmax += molPadding;
}

void calcNumGridPoints(double xrange, double yrange, double zrange,
                       double gridSpacing, unsigned int &numX,
                       unsigned int &numY, unsigned int &numZ) {
  numX = 1 + static_cast<int>(xrange / gridSpacing);
  numY = 1 + static_cast<int>(yrange / gridSpacing);
  numZ = 1 + static_cast<int>(zrange / gridSpacing);
}

double getStandardAtomRadius(const unsigned int atomicNum) {
  // Mostly they will be carbons, so just return that without lookup.
  if (atomicNum == 6) {
    return 1.7;
  }
  if (const auto rad = vdw_radii.find(static_cast<unsigned int>(atomicNum));
      rad != vdw_radii.end()) {
    return rad->second;
  }
  throw ValueErrorException("No VdW radius for atom with Z=" +
                            std::to_string(atomicNum));
}

void extractAlphasAndCoords(const ROMol &mol, int confId,
                            std::vector<double> &alphas,
                            std::vector<RDGeom::Point3D> &coords) {
  static double constexpr alpha_bit = std::numbers::pi / 1.3401;
  alphas.resize(mol.getNumAtoms());
  coords.resize(3 * mol.getNumAtoms());
  const auto &conf = mol.getConformer(confId);
  for (const auto &a : mol.atoms()) {
    const auto rad = getStandardAtomRadius(a->getAtomicNum());
    alphas[a->getIdx()] = alpha_bit / (rad * rad);
    coords[a->getIdx()] = conf.getAtomPos(a->getIdx());
  }
}

double calcDensity(const std::vector<double> &alphas,
                   const std::vector<RDGeom::Point3D> &coords, const double x,
                   const double y, const double z) {
  static double constexpr p =
      4.0 * std::numbers::pi / (3.0 * 1.5514);  // 4 * PI / 3 * lambda
  double dens = 0.0;
  const RDGeom::Point3D gp{x, y, z};
  for (size_t i = 0; i < alphas.size(); ++i) {
    const double alpha_dist = -alphas[i] * (coords[i] - gp).lengthSq();
    dens += p * exp(alpha_dist);
  }
  return dens;
}

void fillGrid(grid3d &grid, const ROMol &mol, double &volume,
              const SpeGSSParams &params, const int confId = -1) {
  double xmin, xmax, ymin, ymax, zmin, zmax;
  getMoleculeExtremes(mol, confId, xmin, xmax, ymin, ymax, zmin, zmax,
                      params.moleculePadding + params.gridSpacing);
  // std::cout << "x : " << xmin << " -> " << xmax << std::endl;
  // std::cout << "y : " << ymin << " -> " << ymax << std::endl;
  // std::cout << "z : " << zmin << " -> " << zmax << std::endl;
  unsigned int numX, numY, numZ;
  calcNumGridPoints(xmax - xmin, ymax - ymin, zmax - zmin, params.gridSpacing,
                    numX, numY, numZ);
  // std::cout << "Num grid points : " << numX << " x " << numY << " x " << numZ
  //           << std::endl;
  grid.set_grid_dimensions(numX, numY, numZ);
  grid.set_r0(xmin, ymin, zmin);
  grid.set_ratio_aspect(params.gridSpacing, params.gridSpacing,
                        params.gridSpacing);

  std::vector<double> alphas;
  std::vector<RDGeom::Point3D> coords;
  double maxDens = 0.0;
  double minDens = std::numeric_limits<double>::max();

  extractAlphasAndCoords(mol, confId, alphas, coords);
  unsigned int numInside = 0;
  for (unsigned int i = 0; i < numX; i++) {
    const double x = xmin + params.gridSpacing * i;
    for (unsigned int j = 0; j < numY; j++) {
      const double y = ymin + params.gridSpacing * j;
      for (unsigned int k = 0; k < numZ; k++) {
        const double z = zmin + params.gridSpacing * k;
        const auto dens = calcDensity(alphas, coords, x, y, z);
        if (dens > maxDens) {
          maxDens = dens;
        }
        if (dens < minDens) {
          minDens = dens;
        }
        if (dens > params.contourValue) {
          numInside++;
        }
        grid.set_grid_value(i, j, k, dens);
      }
    }
  }
  // std::cout << "Density range : " << minDens << " -> " << maxDens <<
  // std::endl;
  volume =
      numInside * params.gridSpacing * params.gridSpacing * params.gridSpacing;
}

void extractVertsAndNormals(surface &surf, std::vector<double> &vertices,
                            std::vector<double> &normals) {
  vertices.resize(surf.get_num_vertices() * 3);
  normals.resize(surf.get_num_vertices() * 3);
  for (unsigned int i = 0, j = 0; i < surf.get_num_vertices(); ++i, j += 3) {
    const auto v = surf.getVertex(i);
    vertices[j] = v[0];
    vertices[j + 1] = v[1];
    vertices[j + 2] = v[2];
    const auto n = surf.getNormal(i);
    normals[j] = n[0];
    normals[j + 1] = n[1];
    normals[j + 2] = n[2];
  }
}

void extractTriangles(surface &surf, std::vector<std::uint64_t> &triangles) {
  triangles.resize(surf.get_num_triangles() * 3);
  for (unsigned int i = 0, j = 0; i < surf.get_num_triangles(); ++i, j += 3) {
    const auto t = surf.getTriangle(i);
    triangles[j] = t[0];
    triangles[j + 1] = t[1];
    triangles[j + 2] = t[2];
  }
}

Edge makeEdge(std::uint64_t v1, std::uint64_t v2) {
  Edge edge;
  if (v1 > v2) {
    edge.first = v1;
    edge.second = v2;
  } else {
    edge.first = v2;
    edge.second = v1;
  }
  return edge;
}

void buildEdges(const std::vector<std::uint64_t> &triangles,
                std::unordered_set<Edge, boost::hash<Edge>> &edges) {
  auto insertEdge = [](std::uint64_t v1, std::uint64_t v2,
                       std::unordered_set<Edge, boost::hash<Edge>> &edges) {
    const auto edge = makeEdge(v1, v2);
    edges.insert(edge);
  };

  for (size_t i = 0; i < triangles.size(); i += 3) {
    insertEdge(triangles[i], triangles[i + 1], edges);
    insertEdge(triangles[i + 1], triangles[i + 2], edges);
    insertEdge(triangles[i + 2], triangles[i], edges);
  }
  // std::cout << "Number of edges : " << edges.size() << std::endl;
}

// find the angle between v1, v2, and v3.
double calcAngle(std::uint64_t v1, std::uint64_t v2, std::uint64_t v3,
                 const std::vector<double> &vertices) {
  const RDGeom::Point3D p1{vertices[3 * v1], vertices[3 * v1 + 1],
                           vertices[3 * v1 + 2]};
  const RDGeom::Point3D p2{vertices[3 * v2], vertices[3 * v2 + 1],
                           vertices[3 * v2 + 2]};
  const RDGeom::Point3D p3{vertices[3 * v3], vertices[3 * v3 + 1],
                           vertices[3 * v3 + 2]};
  const auto d1 = p2.directionVector(p1);
  const auto d2 = p3.directionVector(p1);
  return d1.signedAngleTo(d2);
}

void buildVertexConnections(
    const std::unordered_set<Edge, boost::hash<Edge>> &edges,
    std::vector<std::vector<std::uint64_t>> &vertexConnections) {
  for (const auto &[fst, snd] : edges) {
    vertexConnections[fst].push_back(snd);
    vertexConnections[snd].push_back(fst);
  }
  for (auto &vconn : vertexConnections) {
    vconn.shrink_to_fit();
  }
}

bool oneComponentInGraph(const std::vector<std::vector<std::uint64_t>> &conns) {
  boost::dynamic_bitset<> visited(conns.size());
  std::list<size_t> to_do;
  to_do.push_back(0);
  while (!to_do.empty()) {
    const auto p = to_do.front();
    to_do.pop_front();
    for (const auto c : conns[p]) {
      if (!visited[c]) {
        visited[c] = true;
        to_do.push_back(c);
      }
    }
  }
  if (visited.count() != visited.size()) {
    std::cout << "visited vertices : " << visited.count() << " of "
              << visited.size() << std::endl;
    for (size_t i = 0; i < visited.size(); ++i) {
      if (!visited[i]) {
        std::cout << "Not visited " << i << " : " << conns[i].size() << " : ";
        for (size_t j = 0; j < conns[i].size(); ++j) {
          std::cout << conns[i][j] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
  return visited.count() == visited.size();
}
}  // namespace

SpeGSSSurface::SpeGSSSurface(const ROMol &mol, const SpeGSSParams &params) {
  PRECONDITION(mol.getNumConformers(),
               "molecule must have at least 1 conformer");
  grid3d grid;
  fillGrid(grid, mol, d_volume, params);
  MC33 mc;
  mc.set_grid3d(grid);
  auto surf = std::make_unique<surface>();
  mc.calculate_isosurface(*surf, params.contourValue);
  // std::cout << "raw num vertices : " << surf->get_num_vertices() <<
  // std::endl; std::cout << "raw num triangles : " << surf->get_num_triangles()
  // << std::endl;

  extractVertsAndNormals(*surf, d_vertices, d_normals);
  extractTriangles(*surf, d_triangles);
  finishBuildingSurface();
}

SpeGSSSurface::SpeGSSSurface(const SpeGSSSurface &surf,
                             const std::vector<std::uint64_t> &vtxNums) {
  d_vertices.resize(vtxNums.size() * 3);
  d_normals.resize(vtxNums.size() * 3);
  std::unordered_map<std::uint64_t, std::uint64_t> oldToNew;
  std::uint64_t i = 0;
  for (auto vtxNum : vtxNums) {
    oldToNew[vtxNum] = i;
    d_vertices[3 * i] = surf.getVertices()[3 * vtxNum];
    d_vertices[3 * i + 1] = surf.getVertices()[3 * vtxNum + 1];
    d_vertices[3 * i + 2] = surf.getVertices()[3 * vtxNum + 2];
    d_normals[3 * i] = surf.getNormals()[3 * vtxNum];
    d_normals[3 * i + 1] = surf.getNormals()[3 * vtxNum + 1];
    d_normals[3 * i + 2] = surf.getNormals()[3 * vtxNum + 2];
    ++i;
    // std::cout << vtxNum << " - > " << oldToNew[vtxNum] << std::endl;
  }

  auto newFromOld = [&](std::uint64_t oldNum) -> std::uint64_t {
    const auto it = oldToNew.find(oldNum);
    if (it == oldToNew.end()) {
      return std::numeric_limits<std::uint64_t>::max();
    }
    return it->second;
  };

  auto wantTriangle = [&](const std::vector<std::uint64_t> &triangles,
                          unsigned int triangNum) -> bool {
    if (const auto newNum = newFromOld(triangles[triangNum]);
        newNum == std::numeric_limits<std::uint64_t>::max()) {
      return false;
    }
    if (const auto newNum = newFromOld(triangles[triangNum + 1]);
        newNum == std::numeric_limits<std::uint64_t>::max()) {
      return false;
    }
    if (const auto newNum = newFromOld(triangles[triangNum + 2]);
        newNum == std::numeric_limits<std::uint64_t>::max()) {
      return false;
    }
    return true;
  };
  for (size_t i = 0; i < surf.d_triangles.size(); i += 3) {
    if (wantTriangle(surf.d_triangles, i)) {
      // std::cout << "need " << surf.d_triangles[i] << " -> "
      //           << surf.d_triangles[i + 1] << " -> " << surf.d_triangles[i +
      //           2]
      //           << " vertices" << std::endl;
      d_triangles.push_back(newFromOld(surf.d_triangles[i]));
      d_triangles.push_back(newFromOld(surf.d_triangles[i + 1]));
      d_triangles.push_back(newFromOld(surf.d_triangles[i + 2]));
    }
  }
  finishBuildingSurface();
}

bool SpeGSSSurface::isManifold() const {
  if (strayPointsAndEdges()) {
    std::cout << "fails stray points and edges" << std::endl;
    return false;
  }
  if (!edgeCountsGoodForManifold()) {
    std::cout << "fails edge counts good for Manifold" << std::endl;
    return false;
  }
  if (!oneComponentInMesh()) {
    std::cout << "fails one component in Mesh" << std::endl;
    return false;
  }
  if (!vertexTrianglesAreFan()) {
    std::cout << "fails vertex triangles are fan" << std::endl;
    return false;
  }
  return true;
}

void SpeGSSSurface::writeFile(const std::string &filename) const {
  std::ofstream outs(filename);
  if (!outs.good()) {
    throw(std::string("Could not open file ") + filename);
  }
  writeToStream(outs);
  outs.close();
}

void SpeGSSSurface::writeToStream(std::ostream &os) const {
  os << "VERTICES" << std::endl << getNumVertices() << std::endl << std::endl;
  for (size_t i = 0; i < d_vertices.size(); i += 3) {
    os << d_vertices[i] << " " << d_vertices[i + 1] << " " << d_vertices[i + 2]
       << std::endl;
  }
  os << std::endl
     << "TRIANGLES" << std::endl
     << getNumTriangles() << std::endl
     << std::endl;
  for (size_t i = 0; i < d_triangles.size(); i += 3) {
    os << d_triangles[i] << " " << d_triangles[i + 1] << " "
       << d_triangles[i + 2] << std::endl;
  }
  os << std::endl << "NORMALS" << std::endl;
  for (size_t i = 0; i < d_normals.size(); i += 3) {
    os << d_normals[i] << " " << d_normals[i + 1] << " " << d_normals[i + 2]
       << std::endl;
  }
  os << std::endl << "END" << std::endl;
}

bool SpeGSSSurface::readFile(const std::string &filename) {
  std::ifstream ins(filename);
  if (!ins.good()) {
    throw(std::string("Could not open file ") + filename);
  }
  return readFromStream(ins);
}

bool SpeGSSSurface::readFromStream(std::istream &is) {
  std::string line;
  if (!std::getline(is, line)) {
    return false;
  }
  if (line != "VERTICES") {
    BOOST_LOG(rdErrorLog) << "Bad line " << line << std::endl;
    return false;
  }
  size_t numVerts;
  is >> numVerts;
  d_vertices.resize(numVerts * 3);
  for (size_t i = 0; i < numVerts * 3; i += 3) {
    is >> d_vertices[i] >> d_vertices[i + 1] >> d_vertices[i + 2];
  }
  // new line at end of last vertices line
  if (!std::getline(is, line)) {
    return false;
  }
  // blank line
  if (!std::getline(is, line)) {
    return false;
  }
  if (!std::getline(is, line)) {
    return false;
  }
  if (line != "TRIANGLES") {
    BOOST_LOG(rdErrorLog) << "Bad line " << line << std::endl;
    return false;
  }
  size_t numTriangles;
  is >> numTriangles;
  d_triangles.resize(numTriangles * 3);
  for (size_t i = 0; i < numTriangles * 3; i += 3) {
    is >> d_triangles[i] >> d_triangles[i + 1] >> d_triangles[i + 2];
  }
  if (!std::getline(is, line)) {
    return false;
  }
  if (!std::getline(is, line)) {
    return false;
  }
  if (!std::getline(is, line)) {
    return false;
  }
  if (line != "NORMALS") {
    BOOST_LOG(rdErrorLog) << "Bad line " << line << std::endl;
    return false;
  }
  d_normals.resize(numVerts * 3);
  for (size_t i = 0; i < numVerts * 3; i += 3) {
    is >> d_normals[i] >> d_normals[i + 1] >> d_normals[i + 2];
  }
  finishBuildingSurface();
  return true;
}

void SpeGSSSurface::finishBuildingSurface() {
  buildEdges(d_triangles, d_edges);
  d_vertConns.resize(d_vertices.size() / 3);
  buildVertexConnections(d_edges, d_vertConns);
}

bool SpeGSSSurface::strayPointsAndEdges() const {
  for (const auto &c : d_vertConns) {
    if (c.size() < 2) {
      return true;
    }
  }
  return false;
}

bool SpeGSSSurface::edgeCountsGoodForManifold() const {
  std::unordered_map<Edge, unsigned int, boost::hash<Edge>> edgeCounts;

  auto countEdge = [&edgeCounts](std::uint64_t v1, std::uint64_t v2) {
    const auto edge = makeEdge(v1, v2);
    if (const auto it = edgeCounts.find(edge); it == edgeCounts.end()) {
      edgeCounts[edge] = 1;
    } else {
      ++it->second;
      if (it->second == 3) {
        return false;
      }
    }
    return true;
  };
  for (size_t i = 0; i < d_triangles.size(); i += 3) {
    if (!countEdge(d_triangles[i], d_triangles[i + 1])) {
      return false;
    }
    if (!countEdge(d_triangles[i + 1], d_triangles[i + 2])) {
      return false;
    }
    if (!countEdge(d_triangles[i + 2], d_triangles[i])) {
      return false;
    }
  }
  unsigned int num0 = 0, num1 = 0, num2 = 0, numother = 0;
  for (const auto &val : edgeCounts | std::views::values) {
    if (val == 0) {
      num0++;
    } else if (val == 1) {
      num1++;
    } else if (val == 2) {
      num2++;
    } else {
      numother++;
    }
  }
  // std::cout << "num0 : " << num0 << std::endl;
  // std::cout << "num1 : " << num1 << std::endl;
  // std::cout << "num2 : " << num2 << std::endl;
  // std::cout << "numother : " << numother << std::endl;

  return true;
}

bool SpeGSSSurface::oneComponentInMesh() const {
  return oneComponentInGraph(d_vertConns);
}

bool SpeGSSSurface::vertexTrianglesAreFan() const {
  for (const auto &vconn : d_vertConns) {
    // For all the neighbours of this vertex, find all the edges between them
    // and make a connection table.  If there's a single component in the
    // graph so formed, then it's a fan as all the triangles are contiguous.
    std::vector<std::vector<std::uint64_t>> subConns(vconn.size());
    for (size_t i = 0; i < vconn.size(); ++i) {
      for (size_t j = 0; j < i; ++j) {
        const auto edge = makeEdge(vconn[i], vconn[j]);
        if (auto it = d_edges.find(edge); it != d_edges.end()) {
          subConns[i].push_back(j);
          subConns[j].push_back(i);
        }
      }
    }
    if (!oneComponentInGraph(subConns)) {
      return false;
    }
  }
  return true;
}
}  // namespace SpeGSS
}  // namespace RDKit