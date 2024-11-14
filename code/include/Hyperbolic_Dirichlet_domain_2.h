#ifndef CGAL_HYPERBOLIC_DIRICHLET_DOMAIN_2
#define CGAL_HYPERBOLIC_DIRICHLET_DOMAIN_2

#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include "../include/Anchored_hyperbolic_surface_triangulation_2.h"

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>
#include <CGAL/Circular_kernel_2.h>


typedef CGAL::Exact_rational																							NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>>	Kernel;

typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>         ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                    	Traits;
typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                       Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                 Isometry;
typedef typename Traits::FT                                                 Number;
typedef CGAL::Complex_without_sqrt<Number>                                  Complex_number;
typedef typename Traits::Point_2                                            Point;
typedef CGAL::Hyperbolic_surface_triangulation_2<Traits, CGAL::Anchored_Combinatorial_Map_Attributes<Traits>>   Triangulation;
typedef CGAL::Anchored_hyperbolic_surface_triangulation_2<Traits>                                               Anchored_Triangulation;
typedef typename Anchored_Triangulation::Anchor                                      							Anchor;
typedef CGAL::Combinatorial_map<2,CGAL::Anchored_Combinatorial_Map_Attributes<Traits>>                          CMap;
typedef typename CMap::Dart_const_handle					Dart_handle;
typedef typename CMap::Dart_const_handle					Dart_const_handle;
typedef typename CMap::Dart_const_range						Dart_const_range;

std::vector<std::tuple<Dart_const_handle,Point,Point,Point>> unfold(Triangulation& triangulation);
std::vector<CGAL::Circular_arc_point_2<Kernel>> compute_dirichlet_vertices(Domain& domain);

#endif  //CGAL_HYPERBOLIC_DIRICHLET_DOMAIN_2
