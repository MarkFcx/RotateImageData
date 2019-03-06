cmake_minimum_required(VERSION 3.3 FATAL_ERROR)


project (VolumeTransform)
find_package(VTK REQUIRED)
vtk_module_config(VTK
  vtkCommonCore
  vtkIOImage
  vtkIOXML
  vtkInteractionWidgets
  vtkRenderingVolume${VTK_RENDERING_BACKEND}
)
include(${VTK_USE_FILE})

file(GLOB H_FILES *.h)
file(GLOB CPP_FILES *.cpp)


add_executable(VolumeTransform MACOSX_BUNDLE ${CPP_FILES} ${H_FILES})
target_link_libraries(VolumeTransform ${VTK_LIBRARIES})




