cmake_minimum_required(VERSION 3.1...3.15)

project( demos )

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_BINARY_DIR})

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)


# CGAL and its components
find_package(CGAL REQUIRED COMPONENTS Core Qt6)
find_package(Qt6 QUIET COMPONENTS Widgets)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

if ( NOT CGAL_Qt6_FOUND )

  message(STATUS "This project requires the Qt6 library, and will not be compiled.")
  return()

endif()

# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# ui files, created with Qt Designer
qt6_wrap_ui(UIS drawing_window_description.ui)

add_executable( epsilon_net_demo epsilon_net_demo.cpp window.cpp ../include/Hyperbolic_Dirichlet_domain_2.cpp ${UIS} )
add_executable( dirichlet_and_eps_net dirichlet_and_eps_net.cpp window.cpp ../include/Hyperbolic_Dirichlet_domain_2.cpp ${UIS} )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS epsilon_net_demo )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS dirichlet_and_eps_net )

target_link_libraries( epsilon_net_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets )
target_link_libraries( dirichlet_and_eps_net PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets )

set(CMAKE_BUILD_TYPE "Release")

