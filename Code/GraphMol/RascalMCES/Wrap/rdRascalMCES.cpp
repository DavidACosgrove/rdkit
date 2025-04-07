//
// Copyright (C) David Cosgrove 2023
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
#include <GraphMol/RascalMCES/RascalMCES.h>
#include <GraphMol/RascalMCES/RascalClusterOptions.h>
#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>

#define COMPARE_FUNC_NAME "__call__"
#define CALLBACK_FUNC_NAME "__call__"

namespace python = boost::python;

namespace {

python::list convertVecPairInt(const std::vector<std::pair<int, int>> &vec) {
  python::list pyres;
  for (const auto &p : vec) {
    python::tuple tup = python::make_tuple(p.first, p.second);
    pyres.append(tup);
  }
  return pyres;
}

python::list bondMatches(const RDKit::RascalMCES::RascalResult &self) {
  return convertVecPairInt(self.getBondMatches());
}
python::list atomMatches(const RDKit::RascalMCES::RascalResult &self) {
  return convertVecPairInt(self.getAtomMatches());
}

void largestFragmentOnly(RDKit::RascalMCES::RascalResult &self) {
  self.largestFragOnly();
}

struct RascalResult_wrapper {
  static void wrap() {
    std::string docString = "Used to return RASCAL MCES results.";
    python::class_<RDKit::RascalMCES::RascalResult>(
        "RascalResult", docString.c_str(), python::no_init)
        .def_readonly("smartsString",
                      &RDKit::RascalMCES::RascalResult::getSmarts,
                      "SMARTS string defining the MCES.")
        .def("bondMatches", bondMatches, python::args("self"),
             "A function returning a list of list "
             "of tuples, each inner list containing the matching bonds in the "
             "MCES as tuples of bond indices from mol1 and mol2")
        .def("atomMatches", atomMatches, python::args("self"),
             "Likewise for atoms.")
        .def(
            "largestFragmentOnly", largestFragmentOnly, python::args("self"),
            "Function that cuts the MCES down to the single largest frag.  This cannot be undone.")
        .def_readonly("similarity",
                      &RDKit::RascalMCES::RascalResult::getSimilarity,
                      "Johnson similarity between 2 molecules.")
        .def_readonly("numFragments",
                      &RDKit::RascalMCES::RascalResult::getNumFrags,
                      "Number of fragments in MCES.")
        .def_readonly("largestFragmentSize",
                      &RDKit::RascalMCES::RascalResult::getLargestFragSize,
                      "Number of atoms in largest fragment.")
        .def_readonly("tier1Sim", &RDKit::RascalMCES::RascalResult::getTier1Sim,
                      "The tier 1 similarity estimate.")
        .def_readonly("tier2Sim", &RDKit::RascalMCES::RascalResult::getTier2Sim,
                      "The tier 2 similarity estimate.")
        .def_readonly("timedOut", &RDKit::RascalMCES::RascalResult::getTimedOut,
                      "Whether it timed out.");
  }
};
}  // namespace

namespace RDKit {

python::list findMCESWrapper(const ROMol &mol1, const ROMol &mol2,
                             const python::object &py_opts) {
  RascalMCES::RascalOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RascalMCES::RascalOptions>(py_opts);
  }
  std::vector<RDKit::RascalMCES::RascalResult> results;
  {
    NOGIL gil;
    results = RascalMCES::rascalMCES(mol1, mol2, opts);
  }
  // If the search was cancelled, there will be at least an empty
  // result with the cancelled flag set.  There may be the best
  // MCES found so far as well, as with a timeout, but that's going
  // to be ignored.
  if (results.front().getCancelled()) {
    throw_runtime_error("Rascal cancelled.");
  }

  python::list pyres;
  for (auto &res : results) {
    pyres.append(res);
  }
  return pyres;
}

std::vector<std::shared_ptr<ROMol>> extractMols(python::object mols) {
  std::vector<std::shared_ptr<ROMol>> cmols;
  unsigned int nElems = python::extract<unsigned int>(mols.attr("__len__")());
  cmols.resize(nElems);
  for (unsigned int i = 0; i < nElems; ++i) {
    if (!mols[i]) {
      throw_value_error("molecule is None");
    }
    cmols[i] = python::extract<std::shared_ptr<ROMol>>(mols[i]);
  }
  return cmols;
}

python::list packOutputMols(
    const std::vector<std::vector<unsigned int>> &clusters) {
  python::list pyres;
  for (auto &clus : clusters) {
    python::list mols;
    for (auto &m : clus) {
      mols.append(m);
    }
    pyres.append(mols);
  }
  return pyres;
}

python::list rascalClusterWrapper(python::object mols,
                                  const python::object &py_opts) {
  RascalMCES::RascalClusterOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RascalMCES::RascalClusterOptions>(py_opts);
  }
  auto cmols = extractMols(mols);
  std::vector<RDKit::UINT_VECT> clusters;
  {
    NOGIL gil;
    clusters = RascalMCES::rascalCluster(cmols, opts);
  }
  return packOutputMols(clusters);
}

python::list rascalButinaClusterWrapper(python::object mols,
                                        const python::object &py_opts) {
  RascalMCES::RascalClusterOptions opts;
  if (!py_opts.is_none()) {
    opts = python::extract<RascalMCES::RascalClusterOptions>(py_opts);
  }
  auto cmols = extractMols(mols);
  std::vector<RDKit::UINT_VECT> clusters;
  {
    NOGIL gil;
    clusters = RascalMCES::rascalButinaCluster(cmols, opts);
  }
  return packOutputMols(clusters);
}

struct PyAtomCompareUserData {
  const void *userData;
  python::object pyUserData;
};

struct PyCompareFunctionWrapper
    : public boost::python::wrapper<PyCompareFunctionWrapper> {
  PyCompareFunctionWrapper() {}
  PyCompareFunctionWrapper(PyObject *py_obj) {
    PRECONDITION(py_obj, "PyObject* must not be NULL");
    d_pyObject.reset(
        new python::object(python::handle<>(python::borrowed(py_obj))));
  }
  virtual ~PyCompareFunctionWrapper() {}
  virtual const char *subclassName() const {
    throw std::invalid_argument(
        "subclassName() must be overridden in the "
        "derived class");
  }
  const python::object &pyObject() const { return *d_pyObject; }

 protected:
  void failedToExtractPyObject() const {
    std::stringstream ss;
    ss << "Failed to extract object from " << subclassName() << " subclass";
    PyErr_SetString(PyExc_RuntimeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  void errorNotOverridden() const {
    std::stringstream ss;
    ss << "The " COMPARE_FUNC_NAME "() method must be overridden in the rdFMCS."
       << subclassName() << " subclass";
    PyErr_SetString(PyExc_AttributeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  void errorNotDefined() const {
    // should never happen as the method is virtual but not pure in the C++
    // class
    std::stringstream ss;
    ss << "The " CALLBACK_FUNC_NAME
          "() method must be defined "
          "in the "
       << subclassName() << " subclass";
    PyErr_SetString(PyExc_AttributeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  void errorNotCallable() const {
    std::stringstream ss;
    ss << "The " COMPARE_FUNC_NAME " attribute in the " << subclassName()
       << " subclass is not a callable method";
    PyErr_SetString(PyExc_TypeError, ss.str().c_str());
    python::throw_error_already_set();
  }
  virtual bool hasPythonOverride(const char *attrName) const {
    auto obj = get_override(attrName);
    return PyCallable_Check(obj.ptr());
  }

  void extractPyMCSWrapper() {
    d_pyObjectExtractor.reset(
        new python::extract<PyCompareFunctionWrapper *>(*d_pyObject));
    if (d_pyObjectExtractor->check()) {
      PyObject *callable =
          PyObject_GetAttrString(d_pyObject->ptr(), CALLBACK_FUNC_NAME);
      if (!callable) {
        errorNotDefined();
      }
      if (!PyCallable_Check(callable)) {
        errorNotCallable();
      }
      if (!pyObjectExtract()->hasPythonOverride(CALLBACK_FUNC_NAME)) {
        errorNotOverridden();
      }
    } else {
      std::stringstream ss;
      ss << "expected an instance of the rdFMCS." << subclassName()
         << " subclass";
      PyErr_SetString(PyExc_TypeError, ss.str().c_str());
      python::throw_error_already_set();
    }
  }
  PyCompareFunctionWrapper *pyObjectExtract() const {
    return (*d_pyObjectExtractor)();
  }

 private:
  std::unique_ptr<python::object> d_pyObject;
  std::unique_ptr<python::extract<PyCompareFunctionWrapper *>>
      d_pyObjectExtractor;
};

struct PyAtomCompareFunction : PyCompareFunctionWrapper {
  PyAtomCompareFunction() {}
  PyAtomCompareFunction(PyObject *py_obj) : PyCompareFunctionWrapper(py_obj) {}
  ~PyAtomCompareFunction() {}

  PyAtomCompareFunction *extractPyObject() const {
    auto res = dynamic_cast<PyAtomCompareFunction *>(pyObjectExtract());
    if (!res) {
      failedToExtractPyObject();
    }
    return res;
  }

  const char *subclassName() const { return "PyAtomCompareFunction"; }
  virtual bool operator()(const ROMol &, unsigned int, const ROMol &,
                          unsigned) const {
    errorNotOverridden();
    return false;
  }
};

class PyRascalOptions : public boost::noncopyable {
 public:
  PyRascalOptions() : d_rascalOptions(new RascalMCES::RascalOptions()) {}
  ~PyRascalOptions() = default;

  const RascalMCES::RascalOptions *get() const { return d_rascalOptions.get(); }
  double getSimilarityThreshold() const {
    return d_rascalOptions->similarityThreshold;
  }
  void setSimilarityThreshold(double val) {
    d_rascalOptions->similarityThreshold = val;
  }
  bool getSingleLargestFrag() const {
    return d_rascalOptions->singleLargestFrag;
  }
  void setSingleLargestFrag(bool val) {
    d_rascalOptions->singleLargestFrag = val;
  }
  bool getCompleteAromaticRings() const {
    return d_rascalOptions->completeAromaticRings;
  }
  void setCompleteAromaticRings(bool val) {
    d_rascalOptions->completeAromaticRings = val;
  }
  bool getRingMatchesRingOnly() const {
    return d_rascalOptions->ringMatchesRingOnly;
  }
  void setRingMatchesRingOnly(bool val) {
    d_rascalOptions->ringMatchesRingOnly = val;
  }
  bool getCompleteSmallestRings() const {
    return d_rascalOptions->completeSmallestRings;
  }
  void setCompleteSmallestRings(bool val) {
    d_rascalOptions->completeSmallestRings = val;
  }
  bool getExactConnectionsMatch() const {
    return d_rascalOptions->exactConnectionsMatch;
  }
  void setExactConnectionsMatch(bool val) {
    d_rascalOptions->exactConnectionsMatch = val;
  }
  int getMinFragSize() const { return d_rascalOptions->minFragSize; }
  void setMinFragSize(const int val) { d_rascalOptions->minFragSize = val; }
  int getMaxFragSeparation() const {
    return d_rascalOptions->maxFragSeparation;
  }
  void setMaxFragSeparation(const int val) {
    d_rascalOptions->maxFragSeparation = val;
  }
  bool getAllBestMCESs() const { return d_rascalOptions->allBestMCESs; }
  void setAllBestMCESs(const bool val) { d_rascalOptions->allBestMCESs = val; }
  bool getReturnEmptyMCES() const { return d_rascalOptions->returnEmptyMCES; }
  void setReturnEmptyMCES(const bool val) {
    d_rascalOptions->returnEmptyMCES = val;
  }
  int getTimeout() const { return d_rascalOptions->timeout; }
  void setTimeout(const int val) { d_rascalOptions->timeout = val; }
  unsigned int getMaxBondMatchPairs() const {
    return d_rascalOptions->maxBondMatchPairs;
  }
  void setMaxBondMatchPairs(const unsigned int val) {
    d_rascalOptions->maxBondMatchPairs = val;
  }
  std::string getEquivalentAtoms() const {
    return d_rascalOptions->equivalentAtoms;
  }
  void setEquivalentAtoms(const std::string &val) {
    d_rascalOptions->equivalentAtoms = val;
  }
  bool getIgnoreBondOrders() const { return d_rascalOptions->ignoreBondOrders; }
  void setIgnoreBondOrders(const bool val) {
    d_rascalOptions->ignoreBondOrders = val;
  }
  bool ignoreAtomAromaticity() const {
    return d_rascalOptions->ignoreAtomAromaticity;
  }
  void setIgnoreAtomAromaticity(const bool val) {
    d_rascalOptions->ignoreAtomAromaticity = val;
  }
  unsigned int getMinCliqueSize() const {
    return d_rascalOptions->minCliqueSize;
  }
  void setMinCliqueSize(const unsigned int val) {
    d_rascalOptions->minCliqueSize = val;
  }
  python::object getAtomCompareFunction() const {
    python::object atomCompareFunction;
    if (!d_rascalOptions->atomCompareFunction) {
    }
    return atomCompareFunction;
  }
  void setAtomCompareFunction(PyObject *atomComp) {
    d_rascalOptions->atomCompareFunction = AtomComparePyFunction;
  }

 private:
  static bool AtomComparePyFunction(const ROMol &mol1, unsigned int atom1,
                                    const ROMol &mol2, unsigned int atom2,
                                    void *userData) {
    PRECONDITION(userData, "userData must not be NULL");
    bool res = false;
    res =
        python::call_method<bool>(nullptr, COMPARE_FUNC_NAME, boost::ref(mol1),
                                  atom1, boost::ref(mol2), atom2);
    return res;
  }
  std::unique_ptr<RascalMCES::RascalOptions> d_rascalOptions;
};

BOOST_PYTHON_MODULE(rdRascalMCES) {
  python::scope().attr("__doc__") =
      "Module containing implementation of RASCAL Maximum Common Edge Substructure algorithm.";
  RascalResult_wrapper::wrap();

  std::string docString = "RASCAL Options";
  python::class_<RDKit::RascalMCES::RascalOptions, boost::noncopyable>(
      "RascalOptions", docString.c_str())
      .def_readwrite(
          "similarityThreshold",
          &RDKit::RascalMCES::RascalOptions::similarityThreshold,
          "Threshold below which MCES won't be run.  Between 0.0 and 1.0, default=0.7.")
      .def_readwrite(
          "singleLargestFrag",
          &RDKit::RascalMCES::RascalOptions::singleLargestFrag,
          "Return the just single largest fragment of the MCES.  This is equivalent to running with allBestMCEs=True, finding the result with the largest largestFragmentSize, and calling its largestFragmentOnly method.")
      .def_readwrite(
          "completeAromaticRings",
          &RDKit::RascalMCES::RascalOptions::completeAromaticRings,
          "If True (default), partial aromatic rings won't be returned.")
      .def_readwrite(
          "ringMatchesRingOnly",
          &RDKit::RascalMCES::RascalOptions::ringMatchesRingOnly,
          "If True (default is False), ring bonds won't match non-ring bonds.")
      .def_readwrite(
          "completeSmallestRings",
          &RDKit::RascalMCES::RascalOptions::completeSmallestRings,
          "If True (default is False), only complete rings present in both input molecule's RingInfo will be returned. Implies completeAromaticRings and ringMatchesRingOnly.")
      .def_readwrite(
          "exactConnectionsMatch",
          &RDKit::RascalMCES::RascalOptions::exactConnectionsMatch,
          "If True (default is False), atoms will only match atoms if they have the same\n"
          " number of explicit connections.  E.g. the central atom of\n"
          " C(C)(C) won't match either atom in CC")
      .def_readwrite(
          "minFragSize", &RDKit::RascalMCES::RascalOptions::minFragSize,
          "Imposes a minimum on the number of atoms in a fragment that may be part of the MCES.  Default -1 means no minimum.")
      .def_readwrite(
          "maxFragSeparation",
          &RDKit::RascalMCES::RascalOptions::maxFragSeparation,
          "Maximum number of bonds between fragments in the MCES for both to be reported.  Default -1 means no maximum.  If exceeded, the smaller fragment will be removed.")
      .def_readwrite(
          "allBestMCESs", &RDKit::RascalMCES::RascalOptions::allBestMCESs,
          "If True, reports all MCESs found of the same maximum size.  Default False means just report the first found.")
      .def_readwrite(
          "returnEmptyMCES", &RDKit::RascalMCES::RascalOptions::returnEmptyMCES,
          "If the estimated similarity between the 2 molecules doesn't meet the similarityThreshold, no results are returned.  If you want to know what the"
          " estimates were, set this to True, and examine the tier1Sim and tier2Sim properties of the result then returned.")
      .def_readwrite(
          "timeout", &RDKit::RascalMCES::RascalOptions::timeout,
          "Maximum time (in seconds) to spend on an individual MCESs determination.  Default 60, -1 means no limit.")
      .def_readwrite(
          "maxBondMatchPairs",
          &RDKit::RascalMCES::RascalOptions::maxBondMatchPairs,
          "Too many matching bond (vertex) pairs can cause the process to run out of memory."
          "  The default of 1000 is fairly safe.  Increase with caution, as memory use increases"
          " with the square of this number.  ")
      .def_readwrite("equivalentAtoms",
                     &RDKit::RascalMCES::RascalOptions::equivalentAtoms,
                     "SMARTS strings defining atoms that should"
                     "be considered equivalent. e.g."
                     "[F,Cl,Br,I] so all halogens will match each other."
                     "Space-separated list allowing more than 1"
                     "class of equivalent atoms.")
      .def_readwrite("ignoreBondOrders",
                     &RDKit::RascalMCES::RascalOptions::ignoreBondOrders,
                     "If True, will treat all bonds as the same,"
                     " irrespective of order.  Default=False.")
      .def_readwrite("ignoreAtomAromaticity",
                     &RDKit::RascalMCES::RascalOptions::ignoreAtomAromaticity,
                     "If True, matches atoms solely on atomic number."
                     "  If False, will treat aromatic and aliphatic atoms"
                     " as different.  Default=True.")
      .def_readwrite("minCliqueSize",
                     &RDKit::RascalMCES::RascalOptions::minCliqueSize,
                     "Normally, the minimum clique size is specified"
                     " via the similarityThreshold.  Sometimes it's"
                     " more convenient to specify it directly.  If this"
                     " is > 0, it will over-ride the"
                     " similarityThreshold."
                     "  Note that this refers to the"
                     " minimum number of BONDS in the MCES. Default=0.")
      .def_readwrite(
          "atomCompareFunction", &RDKit::PyAtomCompareFunction::operator(),
          "Function to be used to compare atoms in 2 molecules.  Must"
          " be an instance of a user-defined subclass of rdRascalMCES.AtomCompareFunction.")
      .def_readwrite("atomCompareFunctionUserData",
                     &RDKit::RascalMCES::RascalOptions::atomCompareUserData,
                     "Optional user data for use in atomCompareFunction.");

  docString =
      "Find one or more MCESs between the 2 molecules given.  Returns a list of "
      "RascalResult objects."
      "- mol1"
      "- mol2 The two molecules for which to find the MCES"
      "- opts Optional RascalOptions object changing the default run mode."
      "";
  python::def("FindMCES", &RDKit::findMCESWrapper,
              (python::arg("mol1"), python::arg("mol2"),
               python::arg("opts") = python::object()),
              docString.c_str());

  docString =
      "RASCAL Cluster Options.  Most of these pertain to RascalCluster calculations.  Only similarityCutoff is used by RascalButinaCluster.";
  python::class_<RDKit::RascalMCES::RascalClusterOptions, boost::noncopyable>(
      "RascalClusterOptions", docString.c_str())
      .def_readwrite(
          "similarityCutoff",
          &RDKit::RascalMCES::RascalClusterOptions::similarityCutoff,
          "Similarity cutoff for molecules to be in the same cluster.  Between 0.0 and 1.0, default=0.7.")
      .def_readwrite(
          "minFragSize", &RDKit::RascalMCES::RascalClusterOptions::minFragSize,
          "The minimum number of atoms in a fragment for it to be included in the MCES.  Default=3.")
      .def_readwrite(
          "maxNumFrags", &RDKit::RascalMCES::RascalClusterOptions::maxNumFrags,
          "The maximum number of fragments allowed in the MCES for each pair of molecules. Default=2.  So that the MCES"
          " isn't a lot of small fragments scattered around the molecules giving an inflated estimate of similarity.")
      .def_readwrite(
          "numThreads", &RDKit::RascalMCES::RascalClusterOptions::numThreads,
          "Number of threads to use during clustering.  Default=-1 means all the hardware threads less one.")
      .def_readwrite(
          "a", &RDKit::RascalMCES::RascalClusterOptions::a,
          "The penalty score for each unconnected component in the MCES. Default=0.05.")
      .def_readwrite(
          "b", &RDKit::RascalMCES::RascalClusterOptions::a,
          "The weight of matched bonds over matched atoms. Default=2.")
      .def_readwrite(
          "minIntraClusterSim",
          &RDKit::RascalMCES::RascalClusterOptions::minIntraClusterSim,
          "Two pairs of molecules are included in the same cluster if the similarity between"
          " their MCESs is greater than this.  Default=0.9.")
      .def_readwrite(
          "clusterMergeSim",
          &RDKit::RascalMCES::RascalClusterOptions::clusterMergeSim,
          "Two clusters are merged if the fraction of molecules they have in common is greater than this.  Default=0.6.");

  docString =
      "Use the RASCAL MCES similarity metric to do fuzzy clustering.  Returns a list of lists "
      "of molecules, each inner list being a cluster.  The last cluster is all the "
      "molecules that didn't fit into another cluster (the singletons)."
      "- mols List of molecules to be clustered"
      "- opts Optional RascalOptions object changing the default run mode."
      "";
  python::def("RascalCluster", &RDKit::rascalClusterWrapper,
              (python::arg("mols"), python::arg("opts") = python::object()),
              docString.c_str());
  docString =
      "Use the RASCAL MCES similarity metric to do Butina clustering"
      " (Butina JCICS 39 747-750 (1999)).  Returns a list of lists of molecules,"
      " each inner list being a cluster.  The last cluster is all the"
      " molecules that didn't fit into another cluster (the singletons)."
      "- mols List of molecules to be clustered"
      "- opts Optional RascalOptions object changing the default run mode."
      "";
  python::def("RascalButinaCluster", &RDKit::rascalButinaClusterWrapper,
              (python::arg("mols"), python::arg("opts") = python::object()),
              docString.c_str());

  python::class_<RDKit::PyAtomCompareFunction, boost::noncopyable>(
      "AtomCompareFunction",
      "Base class.  Subclass and override AtomCompareFunction." COMPARE_FUNC_NAME
      "()"
      "to define custom atom comparison functions, then set"
      " RascalOptions.atomCompareFunction to an instance of the subclass.")
      .def(COMPARE_FUNC_NAME, &RDKit::PyAtomCompareFunction::operator(),
           (python::arg("self"), python::arg("mol1"), python::arg("atom1"),
            python::arg("mol2"), python::arg("atom2")),
           "override to implement custom atom comparison");
}

}  // namespace RDKit