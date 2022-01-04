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

#include <GraphMol/MolDraw2D/DrawMolMCH.h>

namespace RDKit {

// ****************************************************************************
DrawMolMCH::DrawMolMCH(
    const ROMol &mol, const std::string &legend, int width, int height,
    MolDrawOptions &drawOptions, DrawText &textDrawer,
    const std::map<int, std::vector<DrawColour>> &highlight_atom_map,
    const std::map<int, std::vector<DrawColour>> &highlight_bond_map,
    const std::map<int, double> &highlight_radii,
    const std::map<int, int> &highlight_linewidth_multipliers,
    int confId)
    : DrawMol(mol, legend, width, height, drawOptions, textDrawer, nullptr,
              nullptr, nullptr, nullptr, nullptr, &highlight_radii, confId),
      mcHighlightAtomMap_(highlight_atom_map),
      mcHighlightBondMap_(highlight_bond_map),
  highlightLinewidthMultipliers_(highlight_linewidth_multipliers) {
  std::cout << "Top of DrawMolMCH c'tor" << std::endl;
}

// ****************************************************************************
void DrawMolMCH::extractHighlights() {
  DrawMol::extractHighlights();
  extractMCHighlights();
}

// ****************************************************************************
void DrawMolMCH::extractMCHighlights() {
  std::cout << "DrawMolMCH::extractHighlights" << std::endl;
  makeBondHighlights();
  std::cout << "after bonds : " << highlights_.size() << std::endl;
  makeAtomHighlights();
  std::cout << "after bonds and atoms : " << highlights_.size() << std::endl;
}

// ****************************************************************************
void DrawMolMCH::makeBondHighlights() {
  for (auto hb : mcHighlightBondMap_) {
    int bond_idx = hb.first;
    int lineWidth = drawOptions_.bondLineWidth;
    if (!drawOptions_.fillHighlights) {
      lineWidth = getHighlightBondWidth(drawOptions_, bond_idx,
                                        &highlightLinewidthMultipliers_);
    }
    auto bond = drawMol_->getBondWithIdx(bond_idx);
    int at1_idx = bond->getBeginAtomIdx();
    int at2_idx = bond->getEndAtomIdx();
    Point2D at1_cds = atCds_[at1_idx];
    Point2D at2_cds = atCds_[at2_idx];
    Point2D perp = calcPerpendicular(at1_cds, at2_cds);
    double rad = 0.7 * drawOptions_.highlightRadius;

    auto draw_adjusted_line = [&](Point2D p1, Point2D p2,
                                  const DrawColour &col) {
      adjustLineEndForHighlight(at1_idx, p2, p1);
      adjustLineEndForHighlight(at2_idx, p1, p2);
      std::vector<Point2D> pts{p1, p2};
      DrawShape *pl = new DrawShapeSimpleLine(
          pts, lineWidth, drawOptions_.scaleBondWidth, col, at1_idx, at2_idx,
          bond->getIdx(), noDash);
      highlights_.emplace_back(std::unique_ptr<DrawShape>(pl));
    };

    if (hb.second.size() < 2) {
      DrawColour col;
      if (hb.second.empty()) {
        col = drawOptions_.highlightColour;
      } else {
        col = hb.second.front();
      }
      if (drawOptions_.fillHighlights) {
        std::vector<Point2D> line_pts;
        line_pts.emplace_back(at1_cds + perp * rad);
        line_pts.emplace_back(at2_cds + perp * rad);
        line_pts.emplace_back(at2_cds - perp * rad);
        line_pts.emplace_back(at1_cds - perp * rad);
        DrawShape *pl = new DrawShapePolyLine(
            line_pts, lineWidth, drawOptions_.scaleBondWidth, col, true,
            at1_idx, at2_idx, bond->getIdx(), noDash);
        highlights_.emplace_back(std::unique_ptr<DrawShape>(pl));
      } else {
        draw_adjusted_line(at1_cds + perp * rad, at2_cds + perp * rad, col);
        draw_adjusted_line(at1_cds - perp * rad, at2_cds - perp * rad, col);
      }
    } else {
      double col_rad = 2.0 * rad / hb.second.size();
      std::cout << "rad = " << rad << "   col_rad = " << col_rad << std::endl;
      if (drawOptions_.fillHighlights) {
        Point2D p1 = at1_cds - perp * rad;
        Point2D p2 = at2_cds - perp * rad;
        std::vector<Point2D> line_pts;
        for (size_t i = 0; i < hb.second.size(); ++i) {
          line_pts.clear();
          line_pts.emplace_back(p1);
          line_pts.emplace_back(p1 + perp * col_rad);
          line_pts.emplace_back(p2 + perp * col_rad);
          line_pts.emplace_back(p2);
          DrawShape *pl = new DrawShapePolyLine(
              line_pts, lineWidth, drawOptions_.scaleBondWidth, hb.second[i],
              true, at1_idx, at2_idx, bond->getIdx(), noDash);
          highlights_.emplace_back(std::unique_ptr<DrawShape>(pl));
          p1 += perp * col_rad;
          p2 += perp * col_rad;
        }
      } else {
        std::cout << "not fill highlights" << std::endl;
        std::vector<DrawColour> cols{hb.second};
        if (cols.size() % 2) {
          draw_adjusted_line(at1_cds, at2_cds, cols[0]);
          cols.erase(cols.begin());
        }
        int step = 0;
        for (size_t i = 0; i < cols.size(); ++i) {
          std::cout << "drawing colour " << i << " : " << cols[i].r
                    << ", " << cols[i].g << ", " << cols[i].b << std::endl;
          // draw even numbers from the bottom, odd from the top
          Point2D offset = perp * (rad - step * col_rad);
          if (!(i % 2)) {
            draw_adjusted_line(at1_cds - offset, at2_cds - offset,
                               cols[i]);
          } else {
            draw_adjusted_line(at1_cds + offset, at2_cds + offset,
                               cols[i]);
            step++;
          }
        }
      }
    }
  }
}

// ****************************************************************************
void DrawMolMCH::makeAtomHighlights() {
  for (auto &ha : mcHighlightAtomMap_) {
    double xradius, yradius;
    Point2D centre;
    int lineWidth = getHighlightBondWidth(drawOptions_, -1, nullptr);
    calcSymbolEllipse(ha.first, centre, xradius, yradius);
    if (ha.second.size() == 1) {
      Point2D offset(xradius, yradius);
      Point2D p1 = centre - offset;
      Point2D p2 = centre + offset;
      std::vector<Point2D> pts{p1, p2};
      DrawShape *ell = new DrawShapeEllipse(
          pts, lineWidth, true, ha.second.front(), drawOptions_.fillHighlights,
          ha.first);
      highlights_.emplace_back(std::unique_ptr<DrawShape>(ell));
    } else {
      double arc_size = 360.0 / double(ha.second.size());
      std::cout << "making arcs for atom " << ha.first << " of size " << arc_size << std::endl;
      double arc_start = 270.0;
      for (size_t i = 0; i < ha.second.size(); ++i){
        double arc_stop = arc_start + arc_size;
        arc_stop = arc_stop >= 360.0 ? arc_stop - 360.0 : arc_stop;
        std::cout << "   " << i << " : " << arc_start << " : " << arc_stop << std::endl;
        std::vector<Point2D> pts{centre, Point2D(xradius, yradius)};
        DrawShape *arc = new DrawShapeArc(
            pts, arc_start, arc_stop, lineWidth, true, ha.second[i],
            drawOptions_.fillHighlights, ha.first);
        highlights_.emplace_back(std::unique_ptr<DrawShape>(arc));
        arc_start += arc_size;
        arc_start = arc_start >= 360.0 ? arc_start - 360.0 : arc_start;
      }
    }

  }
}

// ****************************************************************************
void DrawMolMCH::adjustLineEndForHighlight(int at_idx, Point2D p1,
                                           Point2D &p2) const {
  // this code is transliterated from
  // http://csharphelper.com/blog/2017/08/calculate-where-a-line-segment-and-an-ellipse-intersect-in-c/
  // which has it in C#
  double xradius, yradius;
  Point2D centre;
  calcSymbolEllipse(at_idx, centre, xradius, yradius);
  // cout << "ellipse is : " << centre.x << ", " << centre.y << " rads " <<
  // xradius << " and " << yradius << endl; cout << "p1 = " << p1.x << ", " <<
  // p1.y << endl << "p2 = " << p2.x << ", " << p2.y << endl;
  if (xradius < 1.0e-6 || yradius < 1.0e-6) {
    return;
  }

  // move everything so the ellipse is centred on the origin.
  p1 -= centre;
  p2 -= centre;
  double a2 = xradius * xradius;
  double b2 = yradius * yradius;
  double A =
      (p2.x - p1.x) * (p2.x - p1.x) / a2 + (p2.y - p1.y) * (p2.y - p1.y) / b2;
  double B = 2.0 * p1.x * (p2.x - p1.x) / a2 + 2.0 * p1.y * (p2.y - p1.y) / b2;
  double C = p1.x * p1.x / a2 + p1.y * p1.y / b2 - 1.0;

  auto t_to_point = [&](double t) -> Point2D {
    Point2D ret_val;
    ret_val.x = p1.x + (p2.x - p1.x) * t + centre.x;
    ret_val.y = p1.y + (p2.y - p1.y) * t + centre.y;
    return ret_val;
  };

  double disc = B * B - 4.0 * A * C;
  if (disc < 0.0) {
    // no solutions, leave things as they are.  Bit crap, though.
    return;
  } else if (fabs(disc) < 1.0e-6) {
    // 1 solution
    double t = -B / (2.0 * A);
    // cout << "t = " << t << endl;
    p2 = t_to_point(t);
  } else {
    // 2 solutions - take the one nearest p1.
    double disc_rt = sqrt(disc);
    double t1 = (-B + disc_rt) / (2.0 * A);
    double t2 = (-B - disc_rt) / (2.0 * A);
    // cout << "t1 = " << t1 << "  t2 = " << t2 << endl;
    double t;
    // prefer the t between 0 and 1, as that must be between the original
    // points.  If both are, prefer the lower, as that will be nearest p1,
    // so on the bit of the ellipse the line comes to first.
    bool t1_ok = (t1 >= 0.0 && t1 <= 1.0);
    bool t2_ok = (t2 >= 0.0 && t2 <= 1.0);
    if (t1_ok && !t2_ok) {
      t = t1;
    } else if (t2_ok && !t1_ok) {
      t = t2;
    } else if (t1_ok && t2_ok) {
      t = std::min(t1, t2);
    } else {
      // the intersections are both outside the line between p1 and p2
      // so don't do anything.
      return;
    }
    // cout << "using t = " << t << endl;
    p2 = t_to_point(t);
  }
  // cout << "p2 = " << p2.x << ", " << p2.y << endl;
}

// ****************************************************************************
void DrawMolMCH::calcSymbolEllipse(unsigned int atomIdx, Point2D &centre,
                                   double &xradius, double &yradius) const {
  centre = atCds_[atomIdx];
  xradius = drawOptions_.highlightRadius;
  yradius = xradius;
  if (highlightRadii_ &&
      highlightRadii_->find(atomIdx) != highlightRadii_->end()) {
    xradius = highlightRadii_->find(atomIdx)->second;
    yradius = xradius;
  }

  if (drawOptions_.atomHighlightsAreCircles || !atomLabels_[atomIdx] ||
      atomLabels_[atomIdx]->symbol_.empty()) {
    return;
  }

  double x_min, y_min, x_max, y_max;
  x_min = y_min = std::numeric_limits<double>::max();
  x_max = y_max = -std::numeric_limits<double>::max();
  atomLabels_[atomIdx]->findExtremes(x_min, x_max, y_min, y_max);

  static const double root_2 = sqrt(2.0);
  xradius = std::max(xradius, root_2 * 0.5 * (x_max - x_min));
  yradius = std::max(yradius, root_2 * 0.5 * (y_max - y_min));
  centre.x = 0.5 * (x_max + x_min);
  centre.y = 0.5 * (y_max + y_min);
}

}  // namespace RDKit