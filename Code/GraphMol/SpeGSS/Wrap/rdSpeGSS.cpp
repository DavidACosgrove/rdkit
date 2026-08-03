//
// Copyright (C) David Cosgrove 2026.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDBoost/python.h>
#include <RDBoost/Wrap.h>

#include <GraphMol/ROMol.h>
#include <GraphMol/SpeGSS/SpeGSSParams.h>
#include <GraphMol/SpeGSS/SpeGSSSurface.H>

namespace RDKit {
namespace {
boost::shared_ptr<SpeGSS::SpeGSSSurface> makeSurface1(
    const ROMol &mol, const SpeGSS::SpeGSSParams &params) {
  return boost::make_shared<SpeGSS::SpeGSSSurface>(mol, params);
}

template <class T>
python::list std_vector_to_py_list(const std::vector<T> &v) {
  python::object get_iter = python::iterator<std::vector<T>>();
  python::object iter = get_iter(v);
  python::list l(iter);
  return l;
}

python::list getVerticesHelper(const SpeGSS::SpeGSSSurface &surf) {
  const auto &verts = surf.getVertices();
  return std_vector_to_py_list<double>(verts);
}

python::list getNormalsHelper(const SpeGSS::SpeGSSSurface &surf) {
  const auto &norms = surf.getNormals();
  return std_vector_to_py_list<double>(norms);
}

python::list getTrianglesHelper(const SpeGSS::SpeGSSSurface &surf) {
  const auto &tris = surf.getTriangles();
  return std_vector_to_py_list<std::uint64_t>(tris);
}

}  // namespace

BOOST_PYTHON_MODULE(rdSpeGSS) {
  python::scope().attr("__doc__") =
      "Module containing implementation of SpeGSS - Spectral Geometry"
      " Shape Similarity."
      "  NOTE: This functionality is experimental and the API"
      " and/or results may change in future releases.";

  std::string docString = "SpeGSS parameters.";
  python::class_<SpeGSS::SpeGSSParams, boost::noncopyable>("SpeGSSParams",
                                                           docString.c_str())
      .def_readwrite(
          "moleculePadding", &SpeGSS::SpeGSSParams::moleculePadding,
          "Distance that grid extends beyond molecule when creating surface."
          "  Default=4.5, the same as the default distance cutoff for rdGaussianShape.")
      .def_readwrite(
          "gridSpacing", &SpeGSS::SpeGSSParams::gridSpacing,
          "Spacing between grid points when generating surface.  A smaller number"
          " will give a finer grid end hence a smoother surface with more points."
          "  Default=0.25Angstrom.")
      .def_readwrite("contourValue", &SpeGSS::SpeGSSParams::contourValue,
                     "Contour value for surface.  Default=0.5.")
      .def("__setattr__", &safeSetattr);

  docString = "SpeGSS Surface object.";
  python::class_<SpeGSS::SpeGSSSurface, boost::noncopyable>(
      "Surface", docString.c_str(), python::init<>())
      .def("__init__", python::make_constructor(
                           makeSurface1, python::default_call_policies(),
                           (python::arg("mol"), python::arg("params"))))
      .def("GetNumVertices", &SpeGSS::SpeGSSSurface::getNumVertices,
           "Get the number of vertices in the surface.")
      .def("GetNumTriangles", &SpeGSS::SpeGSSSurface::getNumTriangles,
           "Get the number of triangles in the surface.")
      .def("GetVolume", &SpeGSS::SpeGSSSurface::getVolume,
           "Return the grid-based volume enclosed by the surface.")
      .def("GetVertices", &getVerticesHelper,
           "Get the vertices for the surface as a list of floats.")
      .def("GetNormals", &getNormalsHelper,
           "Get the normals for the surface as a list of floats.")
      .def("GetTriangles", &getTrianglesHelper,
           "Get vertices for the triangles as a list of ints.")
      .def("IsManifold", &SpeGSS::SpeGSSSurface::isManifold,
           "Returns True is the surface is a manifold.")
      .def("WriteFile", &SpeGSS::SpeGSSSurface::writeFile,
           "Write surface to a file.")
      .def("ReadFile", &SpeGSS::SpeGSSSurface::readFile,
           "Read surface from a file.");
}
}  // namespace RDKit
