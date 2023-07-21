#  Copyright (c) 2023 David Cosgrove and other RDKit contributors
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

# These tests are just to check that the Python wrappers are working
# ok.  The bulk of the tests are in the C++ code.
import unittest

from rdkit import Chem
from rdkit.Chem import rdRascalMCES

class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        mol1 = Chem.MolFromSmiles("c1ccccc1Cl")
        mol2 = Chem.MolFromSmiles("c1ccccc1F")
        opts = rdRascalMCES.RascalOptions()

        results = rdRascalMCES.FindMCES(mol1, mol2, opts)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].smartsString, 'c1ccccc1')
        self.assertEqual(len(results[0].bondMatches()), 6)
        self.assertEqual(len(results[0].atomMatches()), 6)

    def test2(self):
        # Test single largest fragment extraction
        ad1 = Chem.MolFromSmiles("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934")
        ad2 = Chem.MolFromSmiles("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336")

        opts = rdRascalMCES.RascalOptions()
        results = rdRascalMCES.FindMCES(ad1, ad2, opts)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].smartsString, 'N(-C)-c1ccc(-CC(=O)-NCCCCCCCCCC):cc1.NC12CC3CC(-C1)-CC(-C2)-C3')
        results[0].largestFragmentOnly()
        self.assertEqual(results[0].smartsString, 'N(-C)-c1ccc(-CC(=O)-NCCCCCCCCCC):cc1')

    def test3(self):
        # Test not specifying options
        mol1 = Chem.MolFromSmiles("c1ccccc1Cl")
        mol2 = Chem.MolFromSmiles("c1ccccc1F")

        results = rdRascalMCES.FindMCES(mol1, mol2)
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].smartsString, 'c1ccccc1')
        self.assertEqual(len(results[0].bondMatches()), 6)
        self.assertEqual(len(results[0].atomMatches()), 6)

    def test4(self):
        # Test setting non-default option
        mol1 = Chem.MolFromSmiles('Oc1cccc2C(=O)C=CC(=O)c12')
        mol2 = Chem.MolFromSmiles('O1C(=O)C=Cc2cc(OC)c(O)cc12')
        results = rdRascalMCES.FindMCES(mol1, mol2)
        self.assertEqual(len(results), 0)

        opts = rdRascalMCES.RascalOptions()
        opts.similarityThreshold = 0.5
        results = rdRascalMCES.FindMCES(mol1, mol2, opts)
        self.assertEqual(len(results), 1)
        

    def test5(self):
        # Check exposure of mol1() and mol2() in RascalResults.
        smi1 = "c1ccccc1C(C(=O)CC)(c1ccccc1)CC(C)N(C)C"
        mol1 = Chem.MolFromSmiles(smi1)
        smi2 = "c1ccccc1C1(CCN(C)CC1)C(=O)OCC"
        mol2 = Chem.MolFromSmiles(smi2)
        opts = rdRascalMCES.RascalOptions()
        opts.similarityThreshold = 0.6
        results = rdRascalMCES.FindMCES(mol1, mol2, opts)
        self.assertEqual(len(results), 1)
        outmol1 = results[0].mol1()
        self.assertEqual(Chem.MolToSmiles(outmol1), 'CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1')
        outmol2 = results[0].mol2()
        self.assertEqual(Chem.MolToSmiles(outmol2), 'CCOC(=O)C1(c2ccccc2)CCN(C)CC1')
        
