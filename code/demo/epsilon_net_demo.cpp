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
#include <CGAL/Timer.h>

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

int main(int argc, char* argv[]){

	// 1. GENERATE THE INPUT

	// std::ifstream data("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_triangulations/input_triangulation_1000.txt");
	// Triangulation triangulation;
	// triangulation.from_stream(data);
	
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
	// std::cout << triangulation.make_delaunay() << std::endl;


	// 2. GET A VERTEX OF THE ANCHOR
	// So that if you run the demo on a same surface but with different values of epsilon,
	// the drawing will be centered at the same vertex and it will look similar.

	Anchor& anchor = triangulation.anchor();
	CMap& cmap = triangulation.combinatorial_map();
	Anchored_Triangulation my_triangulation = Anchored_Triangulation(cmap, anchor);
	CMap& my_cmap = my_triangulation.combinatorial_map();
	Anchor& my_anchor = my_triangulation.Triangulation::anchor();
	Point v0 = my_anchor.vertices[0];

	// 3. COMPUTE EPSILON-NET and display useful info

	double eps;
	if (argc == 1){
		eps = 0.1;
	} else{
		eps = std::stod(argv[1]);
	}
	std::cout << "Computing a " << eps << "-net..." << std::endl;
	CGAL::Timer timer;
	timer.start();
	my_triangulation.epsilon_net(eps);
	timer.stop();
	std::cout << "Done in " << timer.time() << " seconds." << std::endl;
	my_cmap.display_characteristics(std::cout) << std::endl;

	// 4. CHANGE THE ANCHOR
	// To center the drawing at the vertex of step 2

	Anchor anch = std::get<0>(my_triangulation.locate_visibility_walk(v0));
	int index = 0;
	for (int i = 0; i < 3; i++) {
		if(v0 == anch.vertices[i]) {
			index = i;
		}
	}
	my_anchor.dart = anch.dart;
	for (int i = 0; i < 3; i++) {
		my_anchor.vertices[i] = anch.vertices[(i + index) % 3];
		if (i < index) {
			my_anchor.dart = my_triangulation.Triangulation::ccw(my_anchor.dart);
		}
	}
	
	// 5. DRAW

	QApplication app(argc, argv);
	app.setApplicationName("Demo: epsilon-net");
	DemoWindow window;
	window.item().draw_triangulation(my_triangulation);
	window.show();
	QStringList args = app.arguments();
	args.removeAt(0);
	return app.exec();
	
	
	return 0;
	}
