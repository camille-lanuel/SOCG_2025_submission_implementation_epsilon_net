// #define CGAL_CHECK_EXPENSIVE

#include "window.h"

#include "../include/Anchored_hyperbolic_surface_triangulation_2.h"
#include "../include/Hyperbolic_Dirichlet_domain_2.h"

#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>

#include <time.h>

typedef CGAL::Exact_rational																					NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;

typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                        	Traits;
typedef typename Traits::Complex                                                                                ComplexNumber;
typedef typename Traits::Hyperbolic_point_2                                                                     Point;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                                                     Isometry;
typedef CGAL::Hyperbolic_surface_triangulation_2<Traits, CGAL::Anchored_Combinatorial_Map_Attributes<Traits>>   Triangulation;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                                   Factory;

typedef CGAL::Anchored_hyperbolic_surface_triangulation_2<Traits>                                               Anchored_Triangulation;
typedef typename Triangulation::Anchor                                                                          Anchor;
typedef CGAL::Combinatorial_map<2,CGAL::Anchored_Combinatorial_Map_Attributes<Traits>>                          CMap;

int main(int argc, char** argv){
	Domain domain;
	if (argc<=2){
		int seed = time(NULL);
		std::cout << "Generating surface with random seed " << seed << "..." << std::endl;
		Factory factory = Factory(seed);
		domain = factory.generate_domain_g2();
	} else{
		int seed = atoi(argv[2]);
		std::cout << "Generating surface with seed " << seed << "..." << std::endl;
		Factory factory = Factory(seed);
		domain = factory.generate_domain_g2();
	}
	Triangulation triangulation = Triangulation(domain);
	triangulation.make_delaunay();

	Anchor& anchor = triangulation.anchor(); 
	CMap& cmap = triangulation.combinatorial_map();
	Anchored_Triangulation my_triangulation = Anchored_Triangulation(cmap, anchor);

	CMap& my_cmap = my_triangulation.combinatorial_map();
	Anchor& my_anchor = my_triangulation.Triangulation::anchor();
	Point v0 = my_anchor.vertices[0];

	QApplication app(argc, argv);
	app.setApplicationName("Demo: Dirichlet domain and epsilon-net");
	DemoWindow window;
	window.item().draw_dirichlet(domain);

	double eps;
	if (argc == 1){
		eps = 0.1;
	} else{
		eps = std::stod(argv[1]);
	}
	std::cout << "Computing a " << eps << "-net..." << std::endl;
	my_triangulation.epsilon_net(eps);
	my_cmap.display_characteristics(std::cout) << std::endl;

	// change the triangulation's anchor such that it has v0 as its first vertex
	Anchor anch = std::get<0>(my_triangulation.locate_visibility_walk(v0));
	int index = 0;
	for (int i = 0; i < 3; i++) {
		if(v0 == anch.vertices[i]) {
			index = i;
		}
	}

	// set the triangulation's anchor s.t. its vertex #0 is v0 AND s.t. v0 is translated at 0
	my_anchor.dart = anch.dart;
	Isometry center_anchor = CGAL::hyperbolic_translation<Traits>(v0);
	for (int i = 0; i < 3; i++) {
		my_anchor.vertices[i] = center_anchor.evaluate(anch.vertices[(i + index) % 3]);
		if (i < index) {
			my_anchor.dart = my_triangulation.Triangulation::ccw(my_anchor.dart);
		}
	}
	
	// compute lifts by visiting them following a BFS algo
	std::queue<Anchor> anchors_queue;
	std::vector<Anchor> anchors;
	anchors_queue.push(my_anchor);
	size_t visited_mark = my_cmap.get_new_mark();

	while(my_cmap.number_of_unmarked_darts(visited_mark) > 0){
		Anchor current = anchors_queue.front();
		auto vdart = current.dart;
		if(!my_cmap.is_marked(vdart, visited_mark)){
			anchors.push_back(current);
			Anchor neighbor;
			for (int i=0; i<3; i++){
				ComplexNumber cross_ratio = my_triangulation.Triangulation::get_cross_ratio(vdart);
				Point c = current.vertices[i%3];
				Point a = current.vertices[(i+1)%3];
				Point b = current.vertices[(i+2)%3];
				Point d = my_triangulation.Triangulation::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
				neighbor = my_triangulation.create_anchor(my_triangulation.Triangulation::opposite(vdart), a, c, d);
				anchors_queue.push(neighbor);
				my_cmap.mark(vdart, visited_mark);
				vdart = my_triangulation.Triangulation::ccw(vdart);
			}
		}
		anchors_queue.pop();
	}

	my_cmap.free_mark(visited_mark);
	window.item().draw_triangles(anchors);
	// fin

	window.show();

	QStringList args = app.arguments();
	args.removeAt(0);
	return app.exec();
	}
