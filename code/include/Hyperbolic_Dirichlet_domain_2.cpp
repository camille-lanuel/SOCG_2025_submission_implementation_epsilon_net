#include "Hyperbolic_Dirichlet_domain_2.h"

// Input: triangulation with a single vertex v
// Output: lift of this triangulation with one lift of v mapped to the origin of the Poincaré disk and all its incident triangles arout it.
std::vector<std::tuple<Dart_const_handle,Point,Point,Point>> unfold(Triangulation& triangulation)
{
    Anchor& anchor = triangulation.anchor();
    CMap& cmap = triangulation.combinatorial_map();

    std::vector<std::tuple<Dart_const_handle,Point,Point,Point>> lifted_triangles;  // vector of lifted triangles (future output)
    std::map<Dart_const_handle, Point> positions;  // map that will contain the computed lift of the vertex of each dart

    // create a mark for visited darts
    size_t visited_darts_mark = cmap.get_new_mark();
    cmap.unmark_all(visited_darts_mark);

    // translate the first vertex of the anchor to the origin of the Poincaré disk, and add the positions of the translated vertices of the anchor
    Isometry center_the_drawing = CGAL::hyperbolic_translation<Traits>(anchor.vertices[0]);
    positions[anchor.dart] = center_the_drawing.evaluate(anchor.vertices[0]);
    positions[triangulation.const_ccw(anchor.dart)] = center_the_drawing.evaluate(anchor.vertices[1]);
    positions[triangulation.const_cw(anchor.dart)] = center_the_drawing.evaluate(anchor.vertices[2]);
    cmap.mark(anchor.dart, visited_darts_mark);

    // add the first triangle (the translated anchor) to the vector of triangles
    std::tuple<Dart_const_handle,Point,Point,Point> value = std::make_tuple(anchor.dart, positions[anchor.dart], positions[triangulation.const_ccw(anchor.dart)], positions[triangulation.const_cw(anchor.dart)]);
    lifted_triangles.push_back(value);

    // visit all the darts one by one by turning around the central vertex
    Dart_const_handle invader = anchor.dart;
    while( cmap.number_of_unmarked_darts(visited_darts_mark)>1 ){  // >1 because the first triangle appears twice
        Dart_const_handle invaded = triangulation.const_opposite(invader);

        // get the positions of the vertices of the invader's triangle
        Point a = positions[triangulation.const_ccw(invader)];
        Point b = positions[triangulation.const_cw(invader)];
        Point c = positions[invader];
        Complex_number cross_ratio = triangulation.get_cross_ratio(invader);

        // retieve the positions of the invaded's triangle
        positions[invaded] = a;
        positions[triangulation.const_ccw(invaded)] = c;
        Point d = triangulation.fourth_point_from_cross_ratio(a, b, c, cross_ratio);
        positions[triangulation.const_cw(invaded)] = d;

        // add the three vertices to the vector of lifted triangles
        value = std::make_tuple(invaded, Point(a), Point(c), Point(d));
        lifted_triangles.push_back(value);
        cmap.mark(invaded, visited_darts_mark);

        invader = triangulation.const_ccw(invaded);
    }

    cmap.free_mark(visited_darts_mark);

    return lifted_triangles;
}

// Input: Fundamental domain whose vertices are the same point on the surface
// Output: vertices of a Dirichlet domain centered at the origin of the Poincaré disk
std::vector<Traits::Hyperbolic_Voronoi_point_2> compute_dirichlet_vertices(Domain& domain)
{
    Triangulation triangulation = Triangulation(domain);
    triangulation.make_delaunay();

    std::vector<std::tuple<Dart_const_handle,Point,Point,Point>> realized_triangles = unfold(triangulation);
    std::vector<CGAL::Circular_arc_point_2<Kernel>> dirichlet_vertices;
    std::vector<int> dirichlet_pairings;

    ParentTraits gt;
    CGAL::internal::Construct_hyperbolic_circumcenter_CK_2<ParentTraits> chc(gt);

    for (std::tuple<Dart_const_handle,Point,Point,Point>& triangle : realized_triangles){
        CGAL::Circular_arc_point_2<Kernel> circumcenter = chc(std::get<1>(triangle), std::get<2>(triangle), std::get<3>(triangle));
        dirichlet_vertices.push_back(circumcenter);
        int n = dirichlet_vertices.size();
    }

    return(dirichlet_vertices);
}
