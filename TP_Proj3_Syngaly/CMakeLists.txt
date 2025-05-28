# Wymagaj CMake 3.15+
cmake_minimum_required(VERSION 3.15...3.27)

# Ustaw nazwę projektu na "sygnaly"
project(sygnaly LANGUAGES CXX)

# Znajdź zależności
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Dodaj naszą bibliotekę C++ jako moduł Pythona o nazwie "sygnaly"
python_add_library(sygnaly MODULE src/main.cpp WITH_SOABI)

# Połącz moduł z nagłówkami pybind11
target_link_libraries(sygnaly PRIVATE pybind11::headers)

# Zainstaluj skompilowany moduł "sygnaly" w głównym folderze pakietu
install(
  TARGETS sygnaly
  DESTINATION ${PROJECT_NAME}
)