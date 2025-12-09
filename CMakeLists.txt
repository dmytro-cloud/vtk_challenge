
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(vtk_challenge)

find_package(VTK COMPONENTS
  CommonCore
  CommonDataModel
  FiltersSources
  FiltersCore
  FiltersHybrid
  FiltersGeneral
  FiltersStatistics
  CommonColor
  IOCore
  IOInfovis
  InteractionStyle
  RenderingCore
  RenderingFreeType
  RenderingOpenGL2
  RenderingAnnotation
  InteractionWidgets
  REQUIRED
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "vtk_challenge: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
add_executable(vtk_challenge MACOSX_BUNDLE main.cxx SeismicViewer.cxx utils.cxx )
  target_link_libraries(vtk_challenge PRIVATE ${VTK_LIBRARIES}
)
# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS vtk_challenge
  MODULES ${VTK_LIBRARIES}
)
