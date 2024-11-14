#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
// #include "../include/Anchored_hyperbolic_surface_triangulation_2.h"
#include <../include/AHST2_epsilon_net_combinatorics.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>

#include <time.h>
#include <string>
#include <functional>

typedef CGAL::Exact_rational																							NumberType; 
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>>	Kernel;

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


void test_epsilon_net(int surface_start, int surface_end, double epsilons[], int nb_eps)
{
	std::vector<std::string> parameters_names = {"insertions", "total locate", "max locate", "zero locate", "total flips", "max flips"};
	int n_param = parameters_names.size();
	double parameters[surface_end-surface_start][nb_eps][n_param];

	for (int n = surface_start; n < surface_end; n++) {
		std::cout << n << std::endl;
		std::ifstream data("../input_domains/input_domain_"+std::to_string(n)+".txt");
		Domain domain;
		domain.from_stream(data);
		Triangulation triangulation = Triangulation(domain);
		triangulation.make_delaunay();
		Anchor& anchor = triangulation.anchor();
		CMap& cmap = triangulation.combinatorial_map();

		for (int i = 0; i < nb_eps; i++) {
			Anchored_Triangulation at = Anchored_Triangulation(cmap, anchor);
			std::vector<int> res = at.epsilon_net(epsilons[i]);
			for (int j = 0; j < n_param; j++) {
				parameters[n - surface_start][i][j] = res[j];
			}
		}
	}
	
	for (int j = 0; j < n_param; j++) {
		std::cout << "---------- " << parameters_names[j] << " ----------" << std::endl;
		std::cout << "surface/epsilon ";
		for (int i = 0; i < nb_eps; i++) {
			std::cout << epsilons[i] << " ";
		}
		std::cout << '\n';

		for (int n = surface_start; n < surface_end; n++){
			std::cout << n << " ";
			for (int i = 0; i < nb_eps; i++) {
				std::cout << parameters[n - surface_start][i][j] << " ";
			}
			std::cout << '\n';
		}
		std::cout << "\n";
	}
}


// parameters are:
// - number of surfaces you want to test
// - start value for epsilon, end value (excluded), and step
// i.e. 0.5 0 0.1 will test for epsilon in {0.5, 0.4, 0.3, 0.2, 0.1} and 0.5 0.2 0.2 will test for epsilon in {0.5, 0.3}
int main(int argc, char* argv[]){

	if(argc == 1){
		std::cout << "Please provide the following parameters :" << std::endl;
		std::cout << "- range of surfaces indices you want to test (0-180), upper bound excluded" << std::endl;
		std::cout << "- range of epsilon values (start, end excluded, step)" << std::endl;
		std::cout << "Examples: " << std::endl;
		std::cout << "./epsilon_net_combinatorics 0 10 0.5 0 0.1 will run 'epsilon_net' on surfaces 0-9 for epsilon = 0.5, 0.4, 0.3, 0.2, 0.1" << std::endl;
		std::cout << "./epsilon_net_combinatorics 50 100 0.7 0.2 0.2 will run 'epsilon_net' on surfaces 50-99 for epsilon = 0.7, 0.5, 0.3" << std::endl;

		return 0;
	}

	int surface_start = atoi(argv[1]);
	int surface_end = atoi(argv[2]);;
	double eps_start = std::stod(argv[3]);
	double eps_end = std::stod(argv[4]);
	double eps_step = std::stod(argv[5]);

	int nb_eps = std::ceil((eps_start-eps_end)/eps_step);
	double epsilons[nb_eps];
	epsilons[0] = eps_start;
	for (int i=1; i<nb_eps; i++) {
		epsilons[i] = epsilons[i-1]-eps_step;
	}

	test_epsilon_net(surface_start, surface_end, epsilons, nb_eps);

	return 0;
}
