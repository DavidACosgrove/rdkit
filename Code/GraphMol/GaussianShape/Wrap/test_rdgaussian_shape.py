import unittest

from rdkit import Chem
from rdkit.Chem import rdGaussianShape
from rdkit import RDConfig
from rdkit.Geometry import Point3D


datadir = RDConfig.RDBaseDir + '/Code/GraphMol/GaussianShape/test_data'


class TestCase(unittest.TestCase):

  def setUp(self):
    suppl = Chem.SDMolSupplier(datadir + '/test1.sdf')
    self.ref = suppl[0]
    self.probe = suppl[1]

  def test1_Defaults(self):
    tpl = rdGaussianShape.AlignMol(self.ref, self.probe)
    self.assertAlmostEqual(tpl[0], 0.497, places=3)
    self.assertAlmostEqual(tpl[1], 0.760, places=3)
    self.assertAlmostEqual(tpl[2], 0.233, places=3)
    

  def test2_NoColor(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    ovOpts.optimMode = rdGaussianShape.OptimMode.SHAPE_ONLY
    shpOpts=  rdGaussianShape.ShapeInputOptions()
    shpOpts.useColors = False
    tpl = rdGaussianShape.AlignMol(self.ref, self.probe , shpOpts, shpOpts, ovOpts)
    self.assertAlmostEqual(tpl[0], 0.760, places=3)
    self.assertAlmostEqual(tpl[1], 0.760, places=3)
    self.assertAlmostEqual(tpl[2], 0.0, places=3)

  def test3_FromShape(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    shpOpts=  rdGaussianShape.ShapeInputOptions()
    shp = rdGaussianShape.ShapeInput(self.ref, -1, shpOpts, ovOpts)
    self.assertAlmostEqual(shp.ShapeVolume, 591.058, places=3)
    self.assertAlmostEqual(shp.ColorVolume, 31.935, places=3)
    self.assertTrue(type(shp) == rdGaussianShape.ShapeInput)
    tpl = rdGaussianShape.AlignMol(shp, self.probe)
    self.assertAlmostEqual(tpl[0], 0.497, places=3)
    self.assertAlmostEqual(tpl[1], 0.760, places=3)
    self.assertAlmostEqual(tpl[2], 0.233, places=3)

    mol = shp.ShapeToMol()
    self.assertEqual(Chem.MolToSmiles(mol), "CC(=O)Oc1ccccc1C(=O)O")
    mol = shp.ShapeToMol(True)
    self.assertEqual(Chem.MolToSmiles(mol), "CC(=O)Oc1ccccc1C(=O)O.[Xe].[Xe].[Xe].[Xe].[Xe].[Xe]")
    mol = shp.ShapeToMol(False, False)
    self.assertEqual(Chem.MolToSmiles(mol), "C.C.C.C.C.C.C.C.C.C.C.C.C")
    mol = shp.ShapeToMol(True, False)
    self.assertEqual(Chem.MolToSmiles(mol), "C.C.C.C.C.C.C.C.C.C.C.C.C.[Xe].[Xe].[Xe].[Xe].[Xe].[Xe]")

  def test4_customFeatures(self):
    m1 = Chem.MolFromSmiles(
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|")
    opts = rdGaussianShape.ShapeInputOptions()
    opts.customFeatures = [[(1, Point3D(-1.75978, 0.148897, 0), 1.0),
                            (2, Point3D(1.7571, -0.120174, 0.1), 1.0)]]
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    shp = rdGaussianShape.ShapeInput(m1, -1, opts, ovOpts)
    self.assertEqual(shp.NumAtoms, 4)
    self.assertEqual(shp.NumFeatures, 2)
    m2 = Chem.Mol(m1)
    opts2 = rdGaussianShape.ShapeInputOptions()
    opts2.customFeatures = [[(2, Point3D(-1.75978, 0.148897, 0), 1.0, [1, 2, 3]),
                             (1, Point3D(1.7571, -0.120174, 0.1), 1.0, [4, 5, 6])]]
    shp2 = rdGaussianShape.ShapeInput(m2, -1, opts2, ovOpts)
    tpl = rdGaussianShape.AlignShapes(shp, shp2, ovOpts)
    self.assertAlmostEqual(tpl[0], 0.999, places=3)
    self.assertAlmostEqual(tpl[1], 1.000, places=3)
    self.assertAlmostEqual(tpl[2], 0.998, places=3)
    tf = tpl[3]
    self.assertGreater(0.0, tf[0])
    self.assertEqual(1.0, tf[15])

    # check the getter:
    cfs = opts2.customFeatures
    self.assertEqual(len(cfs), 1)
    self.assertEqual(len(cfs[0]), 2)
    self.assertEqual(len(cfs[0][0]), 4)
    self.assertEqual(cfs[0][0][0], 2)
    self.assertEqual(cfs[0][1][3][0], 4)
    self.assertEqual(cfs[0][1][3][1], 5)
    self.assertEqual(cfs[0][1][3][2], 6)

  def test5_customFeatures(self):
    m1 = Chem.MolFromSmiles(
      "O=CC=O |(-1.75978,0.148897,0;-0.621382,-0.394324,0;0.624061,0.3656,.1;1.7571,-0.120174,.1)|")
    opts = rdGaussianShape.ShapeInputOptions()
    opts.customFeatures = [[(1, Point3D(-1.75978, 0.148897,
                                       0), 1.0), (2, Point3D(1.7571, -0.120174, 0.1), 1.0)]]
    m2 = Chem.Mol(m1)
    opts2 = rdGaussianShape.ShapeInputOptions()
    opts2.customFeatures = [[(2, Point3D(-1.75978, 0.148897,
                                        0), 1.0), (1, Point3D(1.7571, -0.120174, 0.1), 1.0)]]
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    tpl = rdGaussianShape.AlignMol(m1, m2, opts, opts2, ovOpts)
    self.assertAlmostEqual(tpl[0], 0.999, places=3)
    self.assertAlmostEqual(tpl[1], 1.000, places=3)
    self.assertAlmostEqual(tpl[2], 0.998, places=3)

  def test6_FixedScore(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    # Just to make sure it's there and returns a value.
    opts = rdGaussianShape.ShapeInputOptions()
    tpl = rdGaussianShape.ScoreMol(self.ref, self.ref, opts, opts, ovOpts)
    self.assertAlmostEqual(tpl[0], 1.0, places=3)
    self.assertAlmostEqual(tpl[1], 1.0, places=3)
    self.assertAlmostEqual(tpl[2], 1.0, places=3)

    opts = rdGaussianShape.ShapeInputOptions()
    opts.useColors = False
    ovOpts.normalize = False
    shp = rdGaussianShape.ShapeInput(self.ref, -1, opts, ovOpts)
    tpl = rdGaussianShape.ScoreMol(shp, self.probe, opts)
    self.assertAlmostEqual(tpl[0], 0.0, places=3)
    self.assertAlmostEqual(tpl[1], 0.0, places=3)
    self.assertAlmostEqual(tpl[2], 0.0, places=3)
    
    opts.useColors = True
    shp1 = rdGaussianShape.ShapeInput(self.probe, -1, opts, ovOpts)
    shp2 = rdGaussianShape.ShapeInput(self.probe, -1, opts, ovOpts)
    tpl = rdGaussianShape.ScoreShape(shp1, shp2, ovOpts)
    self.assertAlmostEqual(tpl[0], 1.0, places=3)
    self.assertAlmostEqual(tpl[1], 1.0, places=3)
    self.assertAlmostEqual(tpl[2], 1.0, places=3)
    self.assertAlmostEqual(tpl[3], 751.0, places=1)
    self.assertAlmostEqual(tpl[4], 42.5, places=1)

  def test7_customAtomRadii(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    opts = rdGaussianShape.ShapeInputOptions()
    opts.allCarbonRadii = False
    opts.atomRadii = [(4, 1.9)]
    shp = rdGaussianShape.ShapeInput(self.ref, -1, opts, ovOpts)
    self.assertAlmostEqual(shp.ShapeVolume, 559.361, places=3)

  def test8_atomSubset(self):
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    opts = rdGaussianShape.ShapeInputOptions()
    opts.atomSubset = [4, 5, 6, 7, 8, 9]
    opts.useColors = False
    opts.allCarbonRadii = False
    shp = rdGaussianShape.ShapeInput(self.ref, -1, opts, ovOpts)
    self.assertAlmostEqual(shp.ShapeVolume, 259.144, places=3)
    self.assertEqual(shp.ColorVolume, 0.0)

  def test9_tversky(self):
    scores = rdGaussianShape.AlignMol(self.ref, self.probe)
    self.assertAlmostEqual(scores[0], 0.497, places=3)
    self.assertAlmostEqual(scores[1], 0.760, places=3)
    self.assertAlmostEqual(scores[2], 0.233, places=3)

    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    ovOpts.simAlpha = 0.95
    ovOpts.simBeta = 0.05
    ref_tversky = rdGaussianShape.AlignMol(self.ref, self.probe, overlayOpts=ovOpts)
    self.assertAlmostEqual(ref_tversky[0], 0.700, places=3)
    self.assertAlmostEqual(ref_tversky[1], 0.968, places=3)
    self.assertAlmostEqual(ref_tversky[2], 0.433, places=3)

    ovOpts.simAlpha = 0.05
    ovOpts.simBeta = 0.95
    fit_tversky = rdGaussianShape.AlignMol(self.ref, self.probe, overlayOpts=ovOpts)
    self.assertAlmostEqual(fit_tversky[0], 0.557, places=3)
    self.assertAlmostEqual(fit_tversky[1], 0.780, places=3)
    self.assertAlmostEqual(fit_tversky[2], 0.335, places=3)

  def test10_multipleConformers(self):

    def loadFile(fileName):
      mol = None
      with open(fileName, "r") as f:
        for line in f:
          line = line.strip()
          if not line:
            continue
          m = Chem.MolFromSmiles(line)
          if mol is None:
            mol = m
          else:
            mol.AddConformer(m.GetConformer())
      return mol
      
    esofile = datadir + "/esomeprazole_multi.smi"
    esomeprazole = loadFile(esofile)

    self.assertEqual(esomeprazole.GetNumConformers(), 10)
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    opts = rdGaussianShape.ShapeInputOptions()
    # Set the to default values, just to show that they can be set.
    opts.sortShapes = True
    opts.includeDummies = True
    opts.shapePruneThreshold = -1.0
    shapes1 = rdGaussianShape.ShapeInput(esomeprazole, -1, opts, ovOpts)
    self.assertEqual(shapes1.NumShapes, 10)

    ranfile = datadir + "/ranitidine_multi.smi"
    ranitidine = loadFile(ranfile)
    shapes2 = rdGaussianShape.ShapeInput(ranitidine, -1, opts, ovOpts)
    self.assertEqual(shapes2.NumShapes, 10)

    bestSim, best1, best2, xform = shapes1.BestSimilarity(shapes2)
    self.assertAlmostEqual(bestSim[0], 0.449, places=3)
    self.assertEqual(best1, 2)
    self.assertEqual(best2, 6)
    self.assertEqual(len(xform), 16)

    stuff = rdGaussianShape.ScoreMoleculeAllConformers(esomeprazole, ranitidine)
    self.assertEqual(len(stuff), 4)
    self.assertEqual(len(stuff[0]), 10)
    self.assertEqual(len(stuff[0][0]), 10)
    self.assertEqual(stuff[1], 8)
    self.assertEqual(stuff[2], 3)
    self.assertAlmostEqual(stuff[0][8][3], 0.449, places=3)

  def test11_subshapeStarts(self):
    ref = Chem.MolFromSmiles(r"[C@H](CCc1ccccc1)(C(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2)CC(C)C |(35.672,41.482,-5.722;34.516,40.842,-6.512;34.843,39.355,-6.7;33.819,38.475,-7.45;33.825,38.414,-8.838;32.951,37.553,-9.53;32.064,36.747,-8.81;32.096,36.799,-7.402;32.985,37.656,-6.73;35.934,42.778,-6.452;36.833,42.858,-7.316;35.175,43.735,-6.275;35.516,41.561,-4.218;36.707,42.096,-3.513;38.055,41.449,-3.859;39.11,42.138,-3.959;37.975,40.129,-3.983;39.134,39.277,-4.298;38.825,38.04,-5.133;37.649,37.934,-5.605;39.788,37.369,-5.652;39.985,38.945,-3.037;39.221,37.953,-2.164;37.934,37.961,-1.823;37.579,36.695,-1.314;38.63,35.975,-1.286;39.736,36.771,-1.642;41.052,36.341,-1.48;41.213,35.042,-0.964;40.095,34.215,-0.69;38.765,34.665,-0.855;36.506,41.966,-2.002;37.6,42.757,-1.31;37.546,44.225,-1.728;37.408,42.58,0.19),wD:0.0,wU:17.21,13.33|")
    probe = Chem.MolFromSmiles(r"N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2 |(37.935,40.394,-3.825;39.119,39.593,-4.13;38.758,38.486,-5.101;37.526,38.337,-5.395;39.716,37.852,-5.605;39.883,39.108,-2.906;39.086,38.098,-2.209;38.093,38.363,-1.34;37.565,37.179,-0.881;38.201,36.136,-1.471;39.193,36.684,-2.308;40.015,35.812,-3.036;39.846,34.441,-2.913;38.844,33.933,-2.075;38.015,34.752,-1.333),wU:1.0|")
    ovOpts = rdGaussianShape.ShapeOverlayOptions()
    ovOpts.startMode = rdGaussianShape.StartMode.SUBSHAPE
    ovOpts.maxSubshapeStarts = 4
    res = rdGaussianShape.AlignMol(ref, probe, overlayOpts=ovOpts)
    self.assertAlmostEqual(res[0], 0.331, places=3)

    subshapeOpts = rdGaussianShape.SubshapeOptions()
    subshapeOpts.clusterStartsThreshold = 1.0
    ovOpts.subshapeOptions = subshapeOpts
    
    subshapeOpts2 = ovOpts.subshapeOptions
    self.assertEqual(subshapeOpts2.clusterStartsThreshold, 1.0)

    
if __name__ == '__main__':
  unittest.main()

  
