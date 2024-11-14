#ifndef WINDOW
#define WINDOW

// Qt headers
#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/DemosMainWindow.h>

// UI generated header
#include "ui_drawing_window_description.h"

// CGAL headers
#include "../include/Hyperbolic_Dirichlet_domain_2.h"
#include "../include/Anchored_hyperbolic_surface_triangulation_2.h"

typedef CGAL::Exact_rational                                          													NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>>	Kernel;

typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                         Traits;
typedef typename Traits::Complex                                                                                ComplexNumber;
typedef typename Traits::Hyperbolic_point_2                                                                     Point;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                                                     Isometry;
typedef CGAL::Hyperbolic_surface_triangulation_2<Traits, CGAL::Anchored_Combinatorial_Map_Attributes<Traits>>   Triangulation;

typedef CGAL::Anchored_hyperbolic_surface_triangulation_2<Traits>                                               Anchored_Triangulation;
typedef typename Triangulation::Anchor                                                                          Anchor;
typedef CGAL::Combinatorial_map<2,CGAL::Anchored_Combinatorial_Map_Attributes<Traits>>                          Combinatorial_Map;

class DemoWindowItem :
    public CGAL::Qt::GraphicsItem
{
    Q_OBJECT // Qt macro for Qt objects
    // (Q_OBJECT does not support templates)
private:
  typedef CGAL::Bbox_2                                Bbox_2; // "Bounding box" : just a box type used for drawing

  // Edges to draw
  std::vector<std::pair<Point,Point>>                 _edges;
  std::vector<std::pair<Point,Point>>                 _dirichlet_edges;

  // Pens for drawing
  QPen                                                _poincare_disk_pen;
  QPen                                                _edges_pen;
  QPen                                                _dirichlet_pen;

  // radius of the poincaré disk
  const int _poincare_disk_radius_in_pixels = 600;
  // Approximation treshold : used to decide when to simplify a computation (ex : draw a line instead of an arc if an hyperbolic segment is very small)
  const double computation_treshold = 0.001;
  const double computation_treshold_squared = computation_treshold*computation_treshold;

public:
  // Initializer
  DemoWindowItem();

  // Qt methods
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  QRectF boundingRect() const;
  void modelChanged();

  // Drawing methods
  void draw_unfold(Triangulation& triangulation);
  void draw_dirichlet(Domain& domain);
  void draw_triangles(std::vector<Anchor> anchors);
  void draw_triangulation(Anchored_Triangulation& triangulation);

private:
  // Sub-methods for drawing edges and vertices
  void draw_point(QPainter* painter, Point position);

  void draw_edge(QPainter* painter, Point source, Point target);
  void draw_line(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y);
  void draw_arc(QPainter* painter, double point_1_x, double point_1_y, double point_2_x, double point_2_y, double center_x, double center_y);

  double deg_angle(double x, double y);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DemoWindow :
    public CGAL::Qt::DemosMainWindow, public Ui::MainWindow
{
    Q_OBJECT // Qt macro for Qt objects
    // (Q_OBJECT does not support templates)
private:
  QGraphicsScene                  _scene;
  DemoWindowItem*                 _item;
  QString                         path;

public:
  DemoWindow();
  DemoWindowItem& item();

  // Events handling
  void keyPressEvent(QKeyEvent* event);

public Q_SLOTS:
  void on_actionSave_as_SVG_triggered();
  void on_actionSave_as_PNG_triggered();
};

#endif // WINDOW
