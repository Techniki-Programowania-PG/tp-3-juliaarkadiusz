cmake_minimum_required(VERSION 3.15)
project(sygnaly LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

add_subdirectory(extern/pybind11)

add_subdirectory(extern/matplotplusplus)

set(SOURCES src/main.cpp)

pybind11_add_module(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PRIVATE matplot)
target_include_directories(${PROJECT_NAME} PRIVATE extern/matplotplusplus/source)

install(
  TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION ${PROJECT_NAME}
  RUNTIME DESTINATION ${PROJECT_NAME}
)
