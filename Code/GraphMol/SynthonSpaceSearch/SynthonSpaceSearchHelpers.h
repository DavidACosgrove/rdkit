//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef SYNTHONSPACESEARCHHELPERS_H
#define SYNTHONSPACESEARCHHELPERS_H

class ShapeInput;

namespace RDKit::SynthonSpaceSearch {

using ShapeSet = std::vector<std::unique_ptr<ShapeInput>>;

}
#endif  // SYNTHONSPACESEARCHHELPERS_H
