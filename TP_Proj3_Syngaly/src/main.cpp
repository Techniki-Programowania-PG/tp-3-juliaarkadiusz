#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Potrzebne do automatycznej konwersji std::vector
#include <vector>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Funkcja generujaca sygnal sinusoidalny
std::vector<double> generate_sine_wave(double frequency, double duration, int sample_rate) {
    std::vector<double> signal;
    int num_samples = static_cast<int>(duration * sample_rate);
    for (int i = 0; i < num_samples; ++i) {
        double time = static_cast<double>(i) / sample_rate;
        signal.push_back(sin(2.0 * M_PI * frequency * time));
    }
    return signal;
}


PYBIND11_MODULE(sygnaly, m) {
    m.doc() = "Biblioteka w C++ do przetwarzania sygnalow"; // Opis modulu

    m.def("generate_sine", &generate_sine_wave, "Generuje sygnal sinusoidalny", pybind11::arg("frequency"), pybind11::arg("duration"), pybind11::arg("sample_rate"));
}