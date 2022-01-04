//
//  Copyright (C) 2014-2021 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//
// This class is a helper used by DrawMol to draw annotation (atom and bond
// notes, for example) onto the molecule.
// It is not part of the public API.

#ifndef RDKIT_DRAWANNOTATION_H
#define RDKIT_DRAWANNOTATION_H

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class DrawText;
class MolDraw2D;

class DrawAnnotation {
  friend class DrawMol;

 public:
  ~DrawAnnotation() = default;

 private:
  DrawAnnotation(const std::string &note, const TextAlignType &align,
                 const std::string &cls, double relFontScale,
                 const Point2D &pos, const DrawColour &colour,
                 DrawText &textDrawer);
  DrawAnnotation(const DrawAnnotation &) = delete;
  DrawAnnotation(const DrawAnnotation &&) = delete;
  DrawAnnotation &operator=(const DrawAnnotation &) = delete;

  // expects xmin etc to be initialised to something sensible.
  void findExtremes(double &xmin, double &xmax,
                    double &ymin, double &ymax, double padding=0.0) const;
  void draw(MolDraw2D &molDrawer) const;
  // this is for debugging almost always.
  void drawRects(MolDraw2D &molDrawer) const;
  void scale(const Point2D &scaleFactor);
  void move(const Point2D &trans);
  bool doesRectClash(const StringRect &rect, double padding) const;

  std::string text_;
  TextAlignType align_;
  std::string class_; // for SVG output, most likely
  double relFontScale_; // fraction of current fontScale to use
  DrawText &textDrawer_;
  Point2D pos_ = Point2D(0.0, 0.0);
  DrawColour colour_ = DrawColour(0.0, 0.0, 0.0, 0.0);
  std::vector<std::shared_ptr<StringRect>> rects_;
}
 ;

} // namespace RDKit
#endif  // RDKIT_DRAWANNOTATION_H