#ifndef CGAL_ANCHORED_HYPERBOLIC_SURFACE_TRIANGULATION_2
#define CGAL_ANCHORED_HYPERBOLIC_SURFACE_TRIANGULATION_2

#include <CGAL/Hyperbolic_surface_triangulation_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>
#include <CGAL/Timer.h>


namespace CGAL{

template<class Traits>
struct Anchored_Combinatorial_Map_Attributes{
        template<class CMap>
        struct Dart_wrapper{
                typedef Cell_attribute<CMap, Complex_without_sqrt<typename Traits::FT>> Edge_attrib;
                typedef Cell_attribute<CMap, typename Hyperbolic_surface_triangulation_2<Traits,Anchored_Combinatorial_Map_Attributes<Traits>>::Anchor> Face_attrib;
                typedef std::tuple<void,Edge_attrib,Face_attrib> Attributes;
                };
        };

enum Locate_type{
        OUTSIDE = -1,
        VERTEX,
        EDGE,
        FACE
};

template<class Traits>
class Anchored_hyperbolic_surface_triangulation_2: public Hyperbolic_surface_triangulation_2<Traits,Anchored_Combinatorial_Map_Attributes<Traits>>{
public:
        typedef Hyperbolic_surface_triangulation_2<Traits,Anchored_Combinatorial_Map_Attributes<Traits>>                     Base;
        typedef Combinatorial_map<2,Anchored_Combinatorial_Map_Attributes<Traits>>                                           Combinatorial_Map;
        typedef typename Hyperbolic_surface_triangulation_2<Traits,Anchored_Combinatorial_Map_Attributes<Traits>>::Anchor    Anchor;

        typedef typename Traits::FT                                                                     Number;
        typedef typename Traits::Complex                                                                ComplexNumber;
        typedef typename Traits::Hyperbolic_point_2                                                     Point;
        typedef typename Traits::Hyperbolic_Voronoi_point_2                                                                                             Voronoi_point;
        typedef typename Combinatorial_Map::Dart_handle                                                 Dart_handle;
        typedef typename Combinatorial_Map::Dart_range                                                  Dart_range;
        typedef typename Combinatorial_Map::template One_dart_per_cell_range<0>                         Vertex_range;
        typedef typename Combinatorial_Map::template One_dart_per_cell_range<1>                         Edge_range;
        typedef typename Combinatorial_Map::template One_dart_per_cell_range<2>                         Face_range;
        typedef typename Combinatorial_Map::Dart_const_handle                                           Dart_const_handle;
        typedef typename Combinatorial_Map::Dart_const_range                                            Dart_const_range;
        typedef typename Combinatorial_Map::template One_dart_per_cell_const_range<1>                   Edge_const_range;
        typedef typename Combinatorial_Map::template One_dart_per_cell_const_range<2>                   Face_const_range;
        typedef Hyperbolic_isometry_2<Traits>                                                                   Isometry;

        const int NB_SIDES = 3;
        const int NULL_INDEX = -1;

        //---------- constructor
        Anchored_hyperbolic_surface_triangulation_2(Combinatorial_Map& cmap, typename Base::Anchor& anch);

        //---------- utilities
        Anchor& anchor(const Dart_handle dart);
        Anchor create_anchor(const Dart_handle dart, const Point& a, const Point& b, const Point& c) const;
        int index_of_dart_in_anchor(const Dart_handle dart);
        bool is_in_anchor(const Dart_handle dart, const Anchor& anch);
        bool is_in_anchor(const Dart_handle dart);
        Dart_handle ith_dart_of_anchor(const int i, const Anchor& anch);
        void display_anchor_vertices(const Anchor& anch, bool round = true) const;
        bool are_triangles_equal(const Anchor& anchor1, const Anchor& anchor2);
        bool is_valid();

        //---------- location and insertion
        std::tuple<Locate_type, int> lies_in_anchor(const Point& query, const Anchor& anch) const;
        Orientation hyperbolic_orientation_2(const Point& p, const Point& q, const Point& r) const;
        bool are_segments_intersecting(const Point& p1, const Point& p2, const Point& q1, const Point& q2) const;
        std::tuple<Anchor, int> locate_visibility_walk(const Point& query, const Anchor& anch);
        std::tuple<Anchor, int> locate_visibility_walk(const Point& query);
        std::tuple<Anchor, int> locate_straight_walk(const Point& query, const Anchor& anch);

        std::vector<Anchor> insert(const Point& query, Anchor& anch);
        std::vector<Anchor> insert(const Point& query);

        //---------- Delaunay related methods
        void flip(Dart_handle dart);
        int make_delaunay();
        std::tuple<int, std::vector<Dart_handle>> delaunay_insert(const Point& query, Anchor& anch);
        std::tuple<int, std::vector<Dart_handle>> delaunay_insert(const Point& query);

        //---------- eps-net methods
        void epsilon_net(double epsilon);
        
        double shortest_loop();

private:
        //---------- constructor
        void set_anchors();

        //---------- utilities
        Anchor& set_new_anchor(Dart_handle dart, const Point& r, const Point& s, const Point& t, const Point& query);
        
        //---------- location and insertion
        std::vector<Anchor> insert_in_face(const Point& query, Anchor& anch);
        std::vector<Anchor> insert_on_edge(const Point& query, Anchor& anch, const int);

        //---------- Delaunay related methods
        void push_flippable_edge(const Dart_handle dart, std::list<Dart_handle>& darts_to_flip);
        std::tuple<int, std::vector<Dart_handle>> restore_delaunay(std::list<Dart_handle>& darts_to_flip);

        //---------- eps-net methods
        Number delta(const Point& u, const Point& v) const;
        Point approx_circumcenter(const Anchor& anch) const;
        Number delta_min(const Anchor& anch, const Point& approx_center) const;
        Number delta_max(const Anchor& anch, const Point& approx_center) const;
        void push_unmarked_triangle(const Dart_handle dart, std::list<Dart_handle>& large_triangles, size_t& large_triangles_mark);
        bool is_epsilon_covering(const double BOUND);
        bool is_epsilon_packing(const double BOUND);
        bool is_epsilon_net(const double BOUND);
};



//---------- constructor


template<class Traits>
Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchored_hyperbolic_surface_triangulation_2(Combinatorial_Map& cmap, typename Base::Anchor& anch)
: Base(cmap, anch)
{       
        set_anchors();
}

// only used in constructor
template<class Traits>
void Anchored_hyperbolic_surface_triangulation_2<Traits>::set_anchors()
{
        std::map<Dart_handle, Point> positions;

        // create a mark for visited darts
        size_t visited_darts_mark = this->_combinatorial_map.get_new_mark();
        this->_combinatorial_map.unmark_all(visited_darts_mark);

        Point a = this->_anchor.vertices[0];
        Point b = this->_anchor.vertices[1];
        Point c = this->_anchor.vertices[2];

        positions[this->_anchor.dart] = a;
        positions[Base::ccw(this->_anchor.dart)] = b;
        positions[Base::cw(this->_anchor.dart)] = c;
        
        Anchor current_anchor = this->_anchor;

        this->_combinatorial_map.template set_attribute<2>(this->_anchor.dart,
                this->_combinatorial_map.template create_attribute<2>(current_anchor));
        this->_combinatorial_map.mark(this->_anchor.dart, visited_darts_mark);
        this->_combinatorial_map.mark(Base::ccw(this->_anchor.dart), visited_darts_mark);
        this->_combinatorial_map.mark(Base::cw(this->_anchor.dart), visited_darts_mark);

        // visit each triangle
        Dart_handle invader = this->_anchor.dart;
        while (this->_combinatorial_map.number_of_unmarked_darts(visited_darts_mark) > 0) {
                Dart_handle invaded = Base::opposite(invader);

                if (not this->_combinatorial_map.is_marked(invaded, visited_darts_mark)) {
                        this->_combinatorial_map.mark(invaded, visited_darts_mark);
                        this->_combinatorial_map.mark(Base::ccw(invaded), visited_darts_mark);
                        this->_combinatorial_map.mark(Base::cw(invaded), visited_darts_mark);

                        // get the information of the invader's triangle
                        a = positions[Base::ccw(invader)];
                        b = positions[Base::cw(invader)];
                        c = positions[invader];
                        ComplexNumber cross_ratio = Base::get_cross_ratio(invader);

                        // retieve the positions of the invaded's triangle
                        positions[invaded] = a;
                        positions[Base::ccw(invaded)] = c;
                        Point d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
                        positions[Base::cw(invaded)] = d;

                        // add the anchor in the face attibute of the cmap
                        current_anchor = create_anchor(invaded, a, c, d);
                        this->_combinatorial_map.template set_attribute<2>(invaded,
                                this->_combinatorial_map.template create_attribute<2>(current_anchor));
                }
                invader = Base::ccw(invaded);
        }

        this->_combinatorial_map.free_mark(visited_darts_mark);
}



//---------- utilities


template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor& Anchored_hyperbolic_surface_triangulation_2<Traits>::anchor(const Dart_handle dart)
{
        return this->_combinatorial_map.template info<2>(dart);
}

template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor Anchored_hyperbolic_surface_triangulation_2<Traits>::create_anchor(const Dart_handle dart, const Point& a, const Point& b, const Point& c) const
{
        Anchor anch;
        anch.dart = dart;
        anch.vertices[0] = a;
        anch.vertices[1] = b;
        anch.vertices[2] = c;
        return anch;
}

template<class Traits>
int Anchored_hyperbolic_surface_triangulation_2<Traits>::index_of_dart_in_anchor(const Dart_handle dart)
{       
        int index = NULL_INDEX;
        Anchor& anch = anchor(dart);
        Dart_handle current_dart = anch.dart;
        for (int i = 0; i < NB_SIDES; i++) {
                if (current_dart == dart) {
                        index = i;
                        break;
                }
                current_dart = Base::ccw(current_dart);
        }
        CGAL_assertion(index != NULL_INDEX);
        return index;
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::is_in_anchor(const Dart_handle dart, const Anchor& anch)
{
        return (anchor(dart).dart == anch.dart);
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::is_in_anchor(const Dart_handle dart)
{
        return is_in_anchor(dart, this->_anchor);
}

template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Dart_handle Anchored_hyperbolic_surface_triangulation_2<Traits>::ith_dart_of_anchor(const int i, const Anchor& anch)
{       
        Dart_handle dart = anch.dart;
        for (int j = 0; j < i; j++) {
                dart = Base::ccw(dart);
        }
        return dart;
}

template<class Traits>
void Anchored_hyperbolic_surface_triangulation_2<Traits>::display_anchor_vertices(const Anchor& anch, bool round) const
{       
        if (round) {
                for (int i = 0; i < NB_SIDES; i++) {
                        std::cout << "vertex " << i << " : ";
                        std::cout << "(" << to_double(anch.vertices[i].x()) << "," << to_double(anch.vertices[i].y()) <<")" << std::endl;
                }
        } else {
                for (int i = 0; i < NB_SIDES; i++) {
                        std::cout << "vertex " << i << " : ";
                        std::cout << "(" << anch.vertices[i].x() << "," << anch.vertices[i].y() <<")" << std::endl;
                }
        }
}

// Output: true iff the triangles described by the anchors are equal (same triangle in the cmap and same vertices in H^2)
template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::are_triangles_equal(const Anchor& anchor1, const Anchor& anchor2)
{
        Anchor& triangle1 = anchor(anchor1.dart);
        Anchor& triangle2 = anchor(anchor2.dart);
        
        // check if anchors correspond to same triangle in the cmap
        if (triangle1.dart != triangle2.dart) {
                return false;
        }

        // check vertices
        bool res = true;
        for (int i = 0; i < NB_SIDES; i++) {
                bool found = false;
                for (int j = 0; j < NB_SIDES; j++) {
                        if (anchor1.vertices[i] ==  anchor2.vertices[j]) {
                                found = true;
                        }
                }
                res = res && found;
        }
        return res;
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::is_valid()
{
        CGAL_assertion(Base::is_valid());
        if (not Base::is_valid()) {
                return false;
        }

        for (typename Face_range::iterator it = this->_combinatorial_map.template one_dart_per_cell<2>().begin();
                it != this->_combinatorial_map.template one_dart_per_cell<2>().end(); ++it) {
                Anchor& current = anchor(it);
                Dart_handle current_dart = current.dart;

                for (int i = 0; i < NB_SIDES; i++) {
                        Dart_handle opposite_dart = Base::opposite(current_dart);
                        Point c1 = current.vertices[i];
                        Point a1 = current.vertices[(i + 1) % NB_SIDES];
                        Point b1 = current.vertices[(i + 2) % NB_SIDES];
                        ComplexNumber cross_ratio = Base::get_cross_ratio(current_dart);
                        Point d1 = Base::fourth_point_from_cross_ratio(a1, b1, c1, cross_ratio);

                        int j = index_of_dart_in_anchor(opposite_dart);
                        Anchor& neighbor = anchor(opposite_dart);
                        Point a2 = neighbor.vertices[j];
                        Point c2 = neighbor.vertices[(j + 1) % NB_SIDES];
                        Isometry pair_sides = isometry_pairing_the_sides<Traits>(a2, c2, a1, c1);
                        CGAL_assertion(pair_sides.evaluate(a2) == a1);
                        CGAL_assertion(pair_sides.evaluate(c2) == c1);

                        Point d2 = pair_sides.evaluate(anchor(opposite_dart).vertices[(j + 2) % NB_SIDES]);
                        CGAL_assertion(d2 == d1);
                        if (d2 != d1) {
                                return false;
                        }

                        current_dart = Base::ccw(current_dart);
                }
        }
        return true;
}



//---------- location and insertion

//TODO put this out of this class
template<class Traits>
Orientation Anchored_hyperbolic_surface_triangulation_2<Traits>::hyperbolic_orientation_2(const Point& p, const Point& q, const Point& r) const
{
        Point origin = Point(0, 0);
        Orientation orientation_to_origin = orientation(p, origin, q);
        if (orientation_to_origin == COLLINEAR) {
                return orientation(p, q, r);
        }

        Traits gt;
        CGAL::internal::Side_of_oriented_hyperbolic_segment_2 orientation_test = gt.side_of_oriented_hyperbolic_segment_2_object();
        Oriented_side orientation_to_disk = orientation_test(p, q, r);
        if (orientation_to_disk == ON_POSITIVE_SIDE) {
                return orientation_to_origin;
        } else if (orientation_to_disk == ON_NEGATIVE_SIDE) {
                if (orientation_to_origin == LEFT_TURN) {
                        return RIGHT_TURN;
                } else{
                        return LEFT_TURN;
                }
        } else{
                return COLLINEAR;
        }
}

// Output: The locate type lt of query relative to the anchor, and an index corresponding to:
// - if lt == FACE: NULL_INDEX (= -1),
// - if lt == EDGE: index of the edge on which query lies,
// - if lt == VERTEX: index of the vertex on which query lies,
// - if lt == OUTSIDE: index of the first edge such that query and the third point of the triangle lies on different sides.
template<class Traits>
std::tuple<Locate_type, int> Anchored_hyperbolic_surface_triangulation_2<Traits>::lies_in_anchor(const Point& query, const Anchor& anch) const
{
        Locate_type lt = FACE;
        int index = NULL_INDEX;
        for (int i = 0; i < NB_SIDES; i++) {
                Orientation ori_query = hyperbolic_orientation_2(anch.vertices[i], anch.vertices[(i + 1) % NB_SIDES], query);
                if (ori_query == RIGHT_TURN) {
                        lt = OUTSIDE;
                        index = i;
                        break;
                }
                if (ori_query == COLLINEAR) {
                        lt = EDGE;
                        index = i;
                        if (hyperbolic_orientation_2(anch.vertices[(i + 1) % NB_SIDES], anch.vertices[(i + 2) % NB_SIDES], query) == COLLINEAR) {
                                lt = VERTEX;
                                index = (i + 1) % NB_SIDES;
                                break;
                        }
                }
        }
        return std::make_tuple(lt, index);
}

// Output: an anchor of the triangle in which query lies,
// and an int corresponding to the number of traversed triangles to find it.
template<class Traits>
std::tuple<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor, int> Anchored_hyperbolic_surface_triangulation_2<Traits>::locate_visibility_walk(const Point& query, const Anchor& anch)
{
        ComplexNumber z_query (query.x(), query.y());
        CGAL_precondition(z_query.squared_modulus() < Number(1));
        CGAL_expensive_precondition(Base::is_delaunay()); 

        // initialisation
        auto [lt, index] = lies_in_anchor(query, anch);
        if (lt != OUTSIDE){
                return std::tuple(anch, 0);
        }

        bool found = false;
        int count = 0;

        Dart_handle dart = anch.dart;
        for (int i = 0; i < index; i++) {
                dart = Base::ccw(dart);
        }
        Point c = anch.vertices[index % NB_SIDES];
        Point a = anch.vertices[(index + 1) % NB_SIDES];
        Point b = anch.vertices[(index + 2) % NB_SIDES];
        Point d;

        // visibility walk
        while (not found) {
                dart = Base::opposite(dart);
                ComplexNumber cross_ratio = Base::get_cross_ratio(dart);
                d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
                if (hyperbolic_orientation_2(c, d, query) == RIGHT_TURN){
                        b = a;
                        a = d;
                        dart = Base::ccw(dart);
                } else if (hyperbolic_orientation_2(d, a, query) == RIGHT_TURN) {
                        b = c;
                        c = d;
                        dart = Base::cw(dart);
                } else {
                        found = true;
                }
                count++;
        }
        Anchor res = create_anchor(dart, a, c, d);
        CGAL_assertion(std::get<0>(lies_in_anchor(query, res)) != OUTSIDE);
        return std::tuple(res, count);
}

template<class Traits>
std::tuple<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor, int> Anchored_hyperbolic_surface_triangulation_2<Traits>::locate_visibility_walk(const Point& query)
{
        return locate_visibility_walk(query, this->_anchor);
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::are_segments_intersecting(const Point& p1, const Point& p2, const Point& q1, const Point& q2) const
{
        Traits gt;
        CGAL::internal::Side_of_oriented_hyperbolic_segment_2 orientation_test = gt.side_of_oriented_hyperbolic_segment_2_object();

        Oriented_side ori_q1 = orientation_test(p1, p2, q1);
        Oriented_side ori_q2 = orientation_test(p1, p2, q2);
        if (ori_q1 * ori_q2 == 1) {
                return false;
        }

        Oriented_side ori_p1 = orientation_test(q1, q2, p1);
        Oriented_side ori_p2 = orientation_test(q1, q2, p2);
        if (ori_p1 * ori_p2 == 1) {
                return false;
        }

        return true;
}

template<class Traits>
std::tuple<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor, int> Anchored_hyperbolic_surface_triangulation_2<Traits>::locate_straight_walk(const Point& query, const Anchor& anch)
{
        ComplexNumber z_query (query.x(), query.y());
        CGAL_precondition(z_query.squared_modulus() < Number(1));

        Dart_handle dart = anch.dart;
        Point p = anch.vertices[0];
        Point r = anch.vertices[1];
        Point l = anch.vertices[2];

        // initialisation
        int count = 0;
        if (hyperbolic_orientation_2(r, p, query) < 0) {
                while (hyperbolic_orientation_2(l, p, query) < 0) {
                        dart = Base::opposite(Base::cw(dart));
                        Point old_l = l;
                        l = Base::fourth_point_from_cross_ratio(p, r, l, Base::get_cross_ratio(dart));
                        r = old_l;
                        count++;
                }
        } else {
                while (hyperbolic_orientation_2(r, p, query) >= 0) {
                        Point old_r = r;
                        r = Base::fourth_point_from_cross_ratio(r, l, p, Base::get_cross_ratio(dart));
                        l = old_r;
                        dart = Base::ccw(Base::opposite(dart));
                        count++;
                }
        }

        // straight walk
        dart = Base::ccw(dart);
        Point p0 = p;
        while (hyperbolic_orientation_2(query, r, l) < 0) {
                ComplexNumber cross_ratio = Base::get_cross_ratio(dart);
                Point s = Base::fourth_point_from_cross_ratio(l, p0, r, cross_ratio);
                if (hyperbolic_orientation_2(s, p, query) < 0) {
                        p0 = r;
                        r = s;
                        dart = Base::cw(Base::opposite(dart));
                } else {
                        p0 = l;
                        l = s;
                        dart = Base::ccw(Base::opposite(dart));
                }
                count++;
        }
        Anchor res = create_anchor(dart, r, l, p0);
        CGAL_assertion(std::get<0>(lies_in_anchor(query, res)) != OUTSIDE);
        return std::tuple(res, count);
}

// Input: dart whose cross-ratio is those of the edge [r, t] and is computed with s and a fourth vertex.
// Output: modified anch such that its vertices are {r, query, t} and its dart represents the edge between r and t.
// The cross-ratio of the edge [r, t] is updated with its value where query is replaced by s.
// Warning: the other cross-ratios are not updated, they will be in insert_on_edge and insert_in_face
template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor& Anchored_hyperbolic_surface_triangulation_2<Traits>::set_new_anchor(Dart_handle dart, const Point& r, const Point& s, const Point& t, const Point& query)
{       
        bool is_triangulation_anchor = is_in_anchor(dart);

        // compute new cross-ratio
        ComplexNumber old_cr = Base::get_cross_ratio(dart);
        Point u = Base::fourth_point_from_cross_ratio(r, s, t, old_cr);
        ComplexNumber new_cr = Base::cross_ratio(r, query, t, u);

        // modify anch and cross-ratio
        Anchor& anch = anchor(dart);
        anch.dart = dart;
        this->_combinatorial_map.template info<1>(dart) = new_cr;
        CGAL_assertion(this->_combinatorial_map.template info<1>(anch.dart)==new_cr);
        anch.vertices[0] = t;
        anch.vertices[1] = r;
        anch.vertices[2] = query;
        this->_combinatorial_map.template info<2>(dart) = anch;

        if (is_triangulation_anchor) {
                this->_anchor = anch;
        }

        return anch;
}

// Output: the three anchors incident to query after its insertion inside anch.
// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
// Note to self: No need to manage triangulation's anchor as it is done in set_new_anchor.
template<class Traits>
std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor> Anchored_hyperbolic_surface_triangulation_2<Traits>::insert_in_face(const Point& query, Anchor& anch)
{       
        std::tuple<Locate_type, int> locate_res = lies_in_anchor(query, anch);
        Locate_type lt = std::get<0>(locate_res);
        int index = std::get<1>(locate_res);
        CGAL_precondition(lt == FACE);
        
        Dart_handle current_dart = anch.dart;
        this->_combinatorial_map.insert_cell_0_in_cell_2(anch.dart);
        std::vector<Anchor> new_anchors;
        for (int i = 0; i < NB_SIDES; i++) {
                Point c = anch.vertices[i];
                Point a = anch.vertices[(i + 1) % NB_SIDES];
                Point b = anch.vertices[(i + 2) % NB_SIDES];
                new_anchors.push_back(set_new_anchor(current_dart, a, b, c, query));
                this->_combinatorial_map.template set_attribute<1>(Base::ccw(current_dart),
                        this->_combinatorial_map.template create_attribute<1>(Base::cross_ratio(query, c, a, b)));
                current_dart = Base::ccw(Base::opposite(Base::ccw(current_dart)));
        }

        return new_anchors;
}

// Output: the four anchors incident to query after its insertion on one edge of anch.
// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
template<class Traits>
std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor> Anchored_hyperbolic_surface_triangulation_2<Traits>::insert_on_edge(const Point& query, Anchor& anch, const int index)
{
        CGAL_precondition(std::get<0>(lies_in_anchor(query, anch)) == EDGE);

        // find dart on which we insert query
        Dart_handle insertion_dart = ith_dart_of_anchor(index, anch);

        // gather information
        bool need_to_change_t_anchor = (is_in_anchor(insertion_dart) || is_in_anchor(Base::opposite(insertion_dart)));
        std::vector<Anchor> new_anchors;
        Point c = anch.vertices[index];
        Point a = anch.vertices[(index + 1) % NB_SIDES];
        Point b = anch.vertices[(index + 2) % NB_SIDES];
        ComplexNumber cross_ratio = Base::get_cross_ratio(insertion_dart);
        Point d = Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);

        Dart_handle dart_ab = Base::ccw(insertion_dart);
        Dart_handle dart_bc = Base::ccw(dart_ab);
        Dart_handle dart_cd = Base::ccw(Base::opposite(insertion_dart));
        Dart_handle dart_da = Base::ccw(dart_cd);

        // insert vertex on edge in the cmap and create new triangles
        this->_combinatorial_map.insert_cell_0_in_cell_1(insertion_dart);
        this->_combinatorial_map.insert_cell_1_in_cell_2(dart_bc, Base::cw(dart_ab));
        this->_combinatorial_map.insert_cell_1_in_cell_2(dart_da, Base::cw(dart_cd));
        CGAL_assertion(Base::ccw(dart_ab) == Base::opposite(Base::cw(dart_bc)));
        CGAL_assertion(Base::ccw(dart_bc) == Base::opposite(Base::cw(dart_cd)));
        CGAL_assertion(Base::ccw(dart_cd) == Base::opposite(Base::cw(dart_da)));
        CGAL_assertion(Base::ccw(dart_da) == Base::opposite(Base::cw(dart_ab)));

        // set new anchors and cross-ratios
        new_anchors.push_back(set_new_anchor(dart_ab, b, c, a, query));
        new_anchors.push_back(set_new_anchor(dart_bc, c, a, b, query));
        new_anchors.push_back(set_new_anchor(dart_cd, d, a, c, query));
        new_anchors.push_back(set_new_anchor(dart_da, a, c, d, query));
        this->_combinatorial_map.template set_attribute<1>(Base::ccw(dart_ab),
                this->_combinatorial_map.template create_attribute<1>(Base::cross_ratio(query, a, b, c)));
        this->_combinatorial_map.template set_attribute<1>(Base::ccw(dart_bc),
                this->_combinatorial_map.template create_attribute<1>(Base::cross_ratio(query, b, c, d)));
        this->_combinatorial_map.template set_attribute<1>(Base::ccw(dart_cd),
                this->_combinatorial_map.template create_attribute<1>(Base::cross_ratio(query, c, d, a)));
        this->_combinatorial_map.template set_attribute<1>(Base::ccw(dart_da),
                this->_combinatorial_map.template create_attribute<1>(Base::cross_ratio(query, d, a, b)));

        if (need_to_change_t_anchor) {
                this->_anchor = anchor(dart_ab);
        }

        return new_anchors;
}

// the darts of the new anchors correspond to the edges of the triangle the point was inserted in 
template<class Traits>
std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor> Anchored_hyperbolic_surface_triangulation_2<Traits>::insert(const Point& query, Anchor& anch)
{
        Anchor locate_anchor = std::get<0>(locate_visibility_walk(query, anch));
        // Anchor locate_anchor = std::get<0>(locate_straight_walk(query, anch));
        auto [lt, index] = lies_in_anchor(query, locate_anchor);
        CGAL_precondition(lt != OUTSIDE);

        std::vector<Anchor> new_anchors;
        if (lt == FACE) {
                new_anchors = insert_in_face(query, locate_anchor);
        } else if (lt == EDGE) {
                new_anchors = insert_on_edge(query, locate_anchor, index);
        }
        return new_anchors;
}

// the darts of the new anchors correspond to the edges of the triangle the point was inserted in
template<class Traits>
std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Anchor> Anchored_hyperbolic_surface_triangulation_2<Traits>::insert(const Point& query)
{
        return insert(query, this->_anchor);
}



//---------- Delaunay related methods


template<class Traits>
void Anchored_hyperbolic_surface_triangulation_2<Traits>::flip(Dart_handle dart)
{
        CGAL_expensive_precondition(is_valid());

        // first gather all the information needed
        Dart_handle a = Base::opposite(dart);  // get a fresh handle
        Dart_handle b = Base::ccw(a);
        Dart_handle c = Base::cw(a);
        Anchor& neighbor = anchor(a);

        Dart_handle d = Base::opposite(a);  // == dart
        Dart_handle e = Base::ccw(d);
        Dart_handle f = Base::cw(d);
        Anchor& anch = anchor(d);

        bool need_to_change_t_anchor = (is_in_anchor(a) || is_in_anchor(d));

        ComplexNumber cross_ratio_AB = Base::get_cross_ratio(e);
        ComplexNumber cross_ratio_BC = Base::get_cross_ratio(f);
        ComplexNumber cross_ratio_CD = Base::get_cross_ratio(b);
        ComplexNumber cross_ratio_DA = Base::get_cross_ratio(c);
        ComplexNumber cross_ratio_AC = Base::get_cross_ratio(a);

        int index = index_of_dart_in_anchor(dart);
        Point C = anch.vertices[index];
        Point A = anch.vertices[(index + 1) % NB_SIDES];
        Point B = anch.vertices[(index + 2) % NB_SIDES];
        Point D = Base::fourth_point_from_cross_ratio(A, B, C, cross_ratio_AC);

        ComplexNumber z_D (D.x(), D.y());
        CGAL_assertion(z_D.squared_modulus() < Number(1));

        // create the new anchors
        Anchor new_anchor = create_anchor(f, B, C, D);
        Anchor new_neighbor = create_anchor(c, D, A, B);

        // compute the new cross ratios
        ComplexNumber one (Number(1), Number(0));
        ComplexNumber cross_ratio_BD = (cross_ratio_AC) / ((cross_ratio_AC) - one) ;
        ComplexNumber cross_ratio_AB_2 = one - (one - (cross_ratio_AB)) * (cross_ratio_AC) ;
        ComplexNumber cross_ratio_BC_2 = one - (one - (cross_ratio_BC)) / (cross_ratio_BD) ;
        ComplexNumber cross_ratio_CD_2 = one - (one - (cross_ratio_CD)) * (cross_ratio_AC) ;
        ComplexNumber cross_ratio_DA_2 = one - (one - (cross_ratio_DA)) / (cross_ratio_BD) ;

        // make the topological flip
        this->_combinatorial_map.template unlink_beta<1>(a);
        this->_combinatorial_map.template unlink_beta<1>(b);
        this->_combinatorial_map.template unlink_beta<1>(c);

        this->_combinatorial_map.template unlink_beta<1>(d);
        this->_combinatorial_map.template unlink_beta<1>(e);
        this->_combinatorial_map.template unlink_beta<1>(f);

        this->_combinatorial_map.template link_beta<1>(b, a);
        this->_combinatorial_map.template link_beta<1>(a, f);
        this->_combinatorial_map.template link_beta<1>(f, b);

        this->_combinatorial_map.template link_beta<1>(e, d);
        this->_combinatorial_map.template link_beta<1>(d, c);
        this->_combinatorial_map.template link_beta<1>(c, e);

        CGAL_assertion(Base::opposite(a) == d);

        // and give the new cross ratios to the edges
        this->_combinatorial_map.template info<1>(a) = cross_ratio_BD;
        this->_combinatorial_map.template info<1>(e) = cross_ratio_AB_2;
        this->_combinatorial_map.template info<1>(f) = cross_ratio_BC_2;
        this->_combinatorial_map.template info<1>(b) = cross_ratio_CD_2;
        this->_combinatorial_map.template info<1>(c) = cross_ratio_DA_2;

        // take care of the particular cases where we need to "flip again"
        if (Base::opposite(e) == b) {
                this->_combinatorial_map.template info<1>(e) = one - (one - cross_ratio_AB_2) * (cross_ratio_AC) ;
        }
        if (Base::opposite(f) == c) {
                this->_combinatorial_map.template info<1>(f) = one - (one - cross_ratio_BC_2) / (cross_ratio_BD) ;
        }
        
        // give the new anchors to the faces
        this->_combinatorial_map.template info<2>(f) = new_anchor;
        this->_combinatorial_map.template info<2>(c) = new_neighbor;
        if (need_to_change_t_anchor) {
                this->_anchor = new_anchor;
        }

        CGAL_expensive_assertion(is_valid());
}

// Pushes dart in the list darts_to_flip if dart is flippable and is not already in the list
template<class Traits>
void Anchored_hyperbolic_surface_triangulation_2<Traits>::push_flippable_edge(const Dart_handle dart, std::list<Dart_handle>& darts_to_flip)
{       
        if (Base::is_delaunay_flippable(dart)) {
                bool already_there = false;
                for (const Dart_handle& dart_to_flip : darts_to_flip) {
                        if (dart_to_flip == dart || dart_to_flip == Base::opposite(dart)) {
                                already_there = true;
                                break;
                        }
                }
                if (not already_there) {
                        darts_to_flip.push_back(dart);
                }
        }
}

// Output: number of flips done to make the triangulation Delaunay given a list of darts to flip
template<class Traits>
std::tuple<int, std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Dart_handle>>
Anchored_hyperbolic_surface_triangulation_2<Traits>::restore_delaunay(std::list<Dart_handle>& darts_to_flip)
{       
        int number_of_flips_done = 0;
        std::vector<Dart_handle> flipped_darts;  // not useless: used in epsilon-net algo
        while (not darts_to_flip.empty()) {
                Dart_handle current_dart = darts_to_flip.front();
                if (Base::is_delaunay_flippable(current_dart)) {
                        flip(current_dart);
                        flipped_darts.push_back(current_dart);
                        number_of_flips_done++;
                        Dart_handle maybe_flippable[4] = {Base::ccw(current_dart), Base::cw(current_dart),
                                                          Base::ccw(Base::opposite(current_dart)), Base::cw(Base::opposite(current_dart))};
                        for (int i = 0; i < 4; i++) {
                                push_flippable_edge(maybe_flippable[i], darts_to_flip);
                        }
                }
                darts_to_flip.pop_front();
        }

        CGAL_expensive_assertion(is_valid());
        CGAL_expensive_assertion(Base::is_delaunay());
        return std::make_tuple(number_of_flips_done, flipped_darts);
}

template<class Traits>
int Anchored_hyperbolic_surface_triangulation_2<Traits>::make_delaunay() {
        int number_of_flips_done = 0;

        Dart_handle edge_to_flip = Base::pick_edge_to_flip();
        while (edge_to_flip != nullptr) {
                flip(edge_to_flip);
                edge_to_flip = Base::pick_edge_to_flip();
                number_of_flips_done++;
        }

        CGAL_expensive_assertion(is_valid());
        CGAL_expensive_assertion(Base::is_delaunay());
        return number_of_flips_done;
}

// Inserts query in the triangulation, with a search starting from the given anchor, and makes the triangulation Delaunay again
// Output: number of flips done to make the triangulation Delaunay after the insertion
template<class Traits>
std::tuple<int, std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Dart_handle>> Anchored_hyperbolic_surface_triangulation_2<Traits>::delaunay_insert(const Point& query, Anchor& anch)
{
        CGAL_expensive_precondition(Base::is_delaunay());

        std::vector<Anchor> new_anchors = insert(query, anch);
        std::list<Dart_handle> darts_to_flip;
        for (int i = 0; i < new_anchors.size(); i++) {
                push_flippable_edge(new_anchors[i].dart, darts_to_flip);
        }

        return restore_delaunay(darts_to_flip);
}

// Same but the search starts from main anchor
template<class Traits>
std::tuple<int, std::vector<typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Dart_handle>> Anchored_hyperbolic_surface_triangulation_2<Traits>::delaunay_insert(const Point& query)
{
        return delaunay_insert(query, this->_anchor);
}



//---------- eps-net methods


template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Number Anchored_hyperbolic_surface_triangulation_2<Traits>::delta(const Point& u, const Point& v) const
{
        Number num = (u.x() - v.x()) * (u.x() - v.x()) + (u.y() - v.y()) * (u.y() - v.y());
        Number den = (1 - (u.x() * u.x() + u.y() * u.y())) * (1 - (v.x() * v.x() + v.y() * v.y()));
        return 2 * num / den;
}

template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Point Anchored_hyperbolic_surface_triangulation_2<Traits>::approx_circumcenter(const Anchor& anch) const
{
        Traits gt;
        CGAL::internal::Construct_hyperbolic_circumcenter_CK_2<Traits> chc(gt);
        Voronoi_point exact_center = chc(anch.vertices[0], anch.vertices[1], anch.vertices[2]);
        Number x = to_double(exact_center.x());
        Number y = to_double(exact_center.y());
        return Point(x, y);
}

template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Number Anchored_hyperbolic_surface_triangulation_2<Traits>::delta_min(const Anchor& anch, const Point& approx_center) const
{
        Number res = Number(999);
        for (int i = 0; i < 3; ++i) {
                res = std::min(delta(anch.vertices[i], approx_center), res);
        }
        return res;
}

template<class Traits>
typename Anchored_hyperbolic_surface_triangulation_2<Traits>::Number Anchored_hyperbolic_surface_triangulation_2<Traits>::delta_max(const Anchor& anch, const Point& approx_center) const
{
        Number res = Number(0);
        for (int i = 0; i < 3; ++i) {
                res = std::max(delta(anch.vertices[i], approx_center), res);
        }
        return res;
}

template<class Traits>
void Anchored_hyperbolic_surface_triangulation_2<Traits>::push_unmarked_triangle(const Dart_handle dart, std::list<Dart_handle>& large_triangles, size_t& large_triangles_mark)
{
        Dart_handle comparison_dart = anchor(dart).dart;
        if (not this->_combinatorial_map.is_marked(comparison_dart, large_triangles_mark)) {
                large_triangles.push_back(comparison_dart);  // so that the darts in the list are exactly the marked ones
                this->_combinatorial_map.mark(comparison_dart, large_triangles_mark);
        }
}

template<class Traits>
void Anchored_hyperbolic_surface_triangulation_2<Traits>::epsilon_net(const double epsilon)
{       
        const double BOUND = std::cosh(epsilon)-1;
        CGAL_assertion(is_epsilon_packing(BOUND));
        size_t large_triangles_mark = this->_combinatorial_map.get_new_mark();
        
        std::list<Dart_handle> large_triangles;
        for (typename Face_range::iterator it = this->_combinatorial_map.template one_dart_per_cell<2>().begin();
                it != this->_combinatorial_map.template one_dart_per_cell<2>().end(); ++it) {
                Anchor& current_anchor = anchor(it);
                Point current_center = approx_circumcenter(current_anchor);
                if (delta_min(current_anchor, current_center) > BOUND) {
                        push_unmarked_triangle(it, large_triangles, large_triangles_mark);
                }
        }

        while (not large_triangles.empty()) {
                Dart_handle current_dart = large_triangles.front();
                Anchor& current_anchor = anchor(current_dart);
                this->_combinatorial_map.unmark(current_dart, large_triangles_mark); 
                large_triangles.pop_front();

                Point current_center = approx_circumcenter(current_anchor);
                if (delta_min(current_anchor, current_center) > BOUND) {
                        std::vector<Anchor> new_anchors = insert(current_center, current_anchor);
                        std::list<Dart_handle> darts_to_flip;
                        for (const Anchor& new_anchor : new_anchors) {
                                push_unmarked_triangle(new_anchor.dart, large_triangles, large_triangles_mark);
                                push_flippable_edge(new_anchor.dart, darts_to_flip);  //the darts of the new anchors correspond to the edges of the triangle the point was inserted in 
                        }

                        std::vector<Dart_handle> flipped_darts = std::get<1>(restore_delaunay(darts_to_flip));
                        for (const Dart_handle& dart : flipped_darts) {
                                push_unmarked_triangle(dart, large_triangles, large_triangles_mark);
                                push_unmarked_triangle(Base::opposite(dart), large_triangles, large_triangles_mark);
                        }
                }
        }
        this->_combinatorial_map.free_mark(large_triangles_mark);
        CGAL_assertion(is_epsilon_covering(BOUND));
        CGAL_assertion(is_epsilon_packing(BOUND));
        // std::cout << "covering:" << is_epsilon_covering(BOUND) << std::endl;
        // std::cout << "packing:" << is_epsilon_packing(0.99*BOUND) << std::endl;
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::is_epsilon_covering(const double BOUND)
{
        bool is_covering = true;
        for (typename Face_range::iterator it = this->_combinatorial_map.template one_dart_per_cell<2>().begin();
                it != this->_combinatorial_map.template one_dart_per_cell<2>().end(); ++it) {
                Anchor& current_anchor = anchor(it);
                Point approx_center = approx_circumcenter(current_anchor);
                if (delta_min(current_anchor, approx_center) > BOUND) {
                        is_covering = false;
                        break;
                }
        }
        return is_covering;
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::is_epsilon_packing(const double BOUND)
{
        bool is_packing = true;
        for (typename Edge_range::iterator it = this->_combinatorial_map.template one_dart_per_cell<1>().begin();
                                                                                it != this->_combinatorial_map.template one_dart_per_cell<1>().end(); ++it) {
                Anchor& current_anchor = anchor(it);
                int index = index_of_dart_in_anchor(it);
                Point a = current_anchor.vertices[index];
                Point b = current_anchor.vertices[(index + 1) % NB_SIDES];
                
                if (delta(a, b) < BOUND) {
                        is_packing = false;  // consider that it's not a packing
                        Dart_handle next = Base::ccw(it);
                        auto doc = this->_combinatorial_map.template darts_of_cell<0>(it);
                        // check if the edge is a loop:
                        // if next and it are darts of the same vertex, then a = b
                        for (auto dart = doc.begin(); dart != doc.end(); ++dart) {
                                if (next == dart) {
                                        is_packing = true;  // the edge is a loop, so actually it's ok
                                }
                        }
                }
        }
        return(is_packing);
}

template<class Traits>
bool Anchored_hyperbolic_surface_triangulation_2<Traits>::is_epsilon_net(const double BOUND)
{
        return is_epsilon_covering(BOUND)&&is_epsilon_packing(BOUND);
}

template<class Traits>
double Anchored_hyperbolic_surface_triangulation_2<Traits>::shortest_loop()
{
        Number min_delta_length = 999;
        for (typename Edge_range::iterator it = this->_combinatorial_map.template one_dart_per_cell<1>().begin();
                                                                                it != this->_combinatorial_map.template one_dart_per_cell<1>().end(); ++it) {
                Anchor& current_anchor = anchor(it);
                int index = index_of_dart_in_anchor(it);
                Point a = current_anchor.vertices[index];
                Point b = current_anchor.vertices[(index + 1) % NB_SIDES];
                
                Dart_handle next = Base::ccw(it);
                auto doc = this->_combinatorial_map.template darts_of_cell<0>(it);
                // check if the edge is a loop:
                // if next and it are darts of the same vertex, then a = b
                for (auto dart = doc.begin(); dart != doc.end(); ++dart) {
                        if (next == dart) {
                                Number delta_length = delta(a, b);  // the edge is a loop
                                min_delta_length = (min_delta_length > delta_length) ? delta_length : min_delta_length;
                        }
                }
        }
        return std::acosh(1+to_double(min_delta_length));
}

}  // namespace CGAL

#endif  //CGAL_ANCHORED_HYPERBOLIC_SURFACE_TRIANGULATION_2