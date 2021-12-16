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

// A set of shapes used in 2D drawing.d  Not part of the public API.

#ifndef RDKIT_DRAWSHAPE_H
#define RDKIT_DRAWSHAPE_H

#include <vector>

#include <Geometry/point.h>
#include <GraphMol/MolDraw2D/MolDraw2DHelpers.h>

namespace RDKit {

class MolDraw2D;
const DashPattern noDash;
const DashPattern dots{2.0, 6.0};
const DashPattern dashes{6,0, 6.0};
const DashPattern shortDashes{2.0, 2.0};

class DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  virtual ~DrawShape() = default;

 protected:
  DrawShape(const std::vector<Point2D> &points, int lineWidth = 2,
            bool scaleLineWidth = false,
            DrawColour lineColour = DrawColour(0, 0, 0), bool fill = false,
            int atom1 = -1, int atom2 = -1, int bond = -1);

  void draw(MolDraw2D &drawer);
  virtual void myDraw(MolDraw2D &drawer) const = 0;
  virtual void findExtremes(double &xmin, double &xmax,
                            double &ymin, double &ymax) const;
  virtual void scale(const Point2D &scale_factor);
  virtual void move(const Point2D &trans);

  std::vector<Point2D> points_;
  int lineWidth_;
  bool scaleLineWidth_;
  DrawColour lineColour_;
  bool fill_;
  int atom1_, atom2_, bond_;
};

class DrawShapeArrow: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeArrow() = default;

 protected:
  DrawShapeArrow(const std::vector<Point2D> &points, int lineWidth = 2,
                 bool scaleLineWidth = false,
                 DrawColour lineColour = DrawColour(0, 0, 0), bool fill = false,
                 double frac = 0.2, double angle = M_PI / 6);
  void myDraw(MolDraw2D &drawer) const;

  double frac_;
  double angle_;
};

class DrawShapeEllipse: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeEllipse() = default;

 protected:
  DrawShapeEllipse(const std::vector<Point2D> &points, int lineWidth = 2,
                   bool scaleLineWidth = false,
                   DrawColour lineColour = DrawColour(0, 0, 0),
                   bool fill = false);
  void myDraw(MolDraw2D &drawer) const;

};

class DrawShapePolyline: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapePolyline() = default;

 protected:
  DrawShapePolyline(const std::vector<Point2D> &points, int lineWidth = 2,
                    bool scaleLineWidth = false,
                    DrawColour lineColour = DrawColour(0, 0, 0),
                    bool fill = false, int atom1 = -1, int atom2 = -1,
                    int bond = -1, DashPattern dashPattern = noDash);
  void myDraw(MolDraw2D &drawer) const;

  DashPattern dashPattern_;
};

class DrawShapeSolidWedge : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeSolidWedge() = default;

 protected:
  DrawShapeSolidWedge(const std::vector<Point2D> points, const DrawColour &col1,
                      const DrawColour &col2, bool inverted, bool splitBonds,
                      int atom1 = -1, int atom2 = -1, int bond = -1);
  void buildTriangles();
  void myDraw(MolDraw2D &drawer) const;
  void scale(const Point2D &scale_factor);

  DrawColour col2_;
  bool inverted_;
  bool splitBonds_;
};

class DrawShapeDashedWedge : protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeDashedWedge() = default;

 protected:
  DrawShapeDashedWedge(const std::vector<Point2D> points,
                       const DrawColour &col1, const DrawColour &col2,
                       bool inverted, int atom1 = -1, int atom2 = -1,
                       int bond = -1);
  void buildLines();
  void myDraw(MolDraw2D &drawer) const;
  void scale(const Point2D &scale_factor);
  void move(const Point2D &trans);

  DrawColour col2_;
  bool inverted_;
  std::vector<DrawColour> lineColours_;
  // for when we re-create the lines when it gets too wide, this is
  // the initial points[0] from the c'tor.
  Point2D at1Cds_;
};

class DrawShapeWavyLine: protected DrawShape {
  friend class MolDraw2D;
  friend class DrawMol;

 public:
  ~DrawShapeWavyLine() = default;

 protected:
  DrawShapeWavyLine(const std::vector<Point2D> points, int lineWidth = 2,
                    bool scaleLineWidth = false,
                    const DrawColour &col1 = DrawColour(0, 0, 0),
                    const DrawColour &col2 = DrawColour(0, 0, 0),
                    int atom1 = -1, int atom2 = -1, int bond = -1);
  void myDraw(MolDraw2D &drawer) const;

  DrawColour col2_;
};

std::vector<Point2D> calcScaledWedgePoints(const Point2D &point,
                                           const Point2D &end1,
                                           const Point2D &end2,
                                           double widthsq);
} // namespace RDKit

#endif  // RDKIT_DRAWSHAPE_H