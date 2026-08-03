#  Copyright (c) 2026 David Cosgrove and other RDKit contributors
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

import os
import tempfile
import unittest

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdSpeGSS


class TestCase(unittest.TestCase):

    def testSurfaceCreation(self):
        params = rdSpeGSS.SpeGSSParams()
        benzene = Chem.MolFromSmiles("c1ccccc1 |(-1.39928,0.010791,-0.0463459;-0.701518,-1.17621,-0.047681;0.676993,-1.18937,-0.00208102;1.38771,-0.0146145,0.0458731;0.719034,1.19151,0.048576;-0.646501,1.17791,0.00285385)|")
        surf1 = rdSpeGSS.Surface(benzene, params)
        self.assertEqual(surf1.GetNumVertices(), 1806)
        self.assertEqual(surf1.GetNumTriangles(), 3608)

        params.gridSpacing = 0.5
        params.contourValue = 1.0
        surf2 = rdSpeGSS.Surface(benzene, params)
        self.assertEqual(surf2.GetNumVertices(), 344)
        self.assertEqual(surf2.GetNumTriangles(), 684)
        verts = surf2.GetVertices()
        self.assertEqual(len(verts), 3 * 344)
        norms = surf2.GetNormals()
        self.assertEqual(len(norms), 3 * 344)
        triangles = surf2.GetTriangles()
        self.assertEqual(len(triangles), 3 * 684)

        surf3 = rdSpeGSS.Surface()
        self.assertEqual(surf3.GetNumVertices(), 0)
        self.assertEqual(surf3.GetNumTriangles(), 0)

    def testSurfaceWriteRead(self):
        params = rdSpeGSS.SpeGSSParams()
        benzene = Chem.MolFromSmiles("c1ccccc1 |(-1.39928,0.010791,-0.0463459;-0.701518,-1.17621,-0.047681;0.676993,-1.18937,-0.00208102;1.38771,-0.0146145,0.0458731;0.719034,1.19151,0.048576;-0.646501,1.17791,0.00285385)|")
        surf1 = rdSpeGSS.Surface(benzene, params)
        self.assertEqual(surf1.GetNumVertices(), 1806)
        self.assertEqual(surf1.GetNumTriangles(), 3608)
        with tempfile.NamedTemporaryFile(delete_on_close=False) as fp:
            fp.close()
            surf1.WriteFile(fp.name)
            surf2 = rdSpeGSS.Surface()
            surf2.ReadFile(fp.name)
            self.assertEqual(surf1.GetNumVertices(), surf2.GetNumVertices())

        
if __name__ == "__main__":
  unittest.main()
