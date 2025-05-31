#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <vector>
#include <string>
#include <cmath>
#include <thread>
#include <complex>
#include <matplot/matplot.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace py = pybind11;



//==================== PODSTAWOWE SYGNALY ====================
std::vector<double> generate_sine_wave(double frequency, double duration, int sample_rate) {
    std::vector<double> signal;
    int num_samples = static_cast<int>(duration * sample_rate);
    signal.reserve(num_samples); // Prealokacja pamięci
    for (int i = 0; i < num_samples; ++i) {
        double time = static_cast<double>(i) / sample_rate;
        signal.push_back(sin(2.0 * M_PI * frequency * time));
    }
    return signal;
}

std::vector<double> generate_cosine_wave(double frequency, double duration, int sample_rate) {
    std::vector<double> signal;
    int num_samples = static_cast<int>(duration * sample_rate);
    signal.reserve(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        double time = static_cast<double>(i) / sample_rate;
        signal.push_back(cos(2.0 * M_PI * frequency * time));
    }
    return signal;
}

std::vector<double> generate_square_wave(double frequency, double duration, int sample_rate) {
    std::vector<double> signal;
    int num_samples = static_cast<int>(duration * sample_rate);
    signal.reserve(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        double time = static_cast<double>(i) / sample_rate;
        if (sin(2.0 * M_PI * frequency * time) > 0) {
            signal.push_back(1.0);
        } else {
            signal.push_back(-1.0);
        }
    }
    return signal;
}

std::vector<double> generate_sawtooth_wave(double frequency, double duration, int sample_rate) {
    std::vector<double> signal;
    int num_samples = static_cast<int>(duration * sample_rate);
    signal.reserve(num_samples);
    for (int i = 0; i < num_samples; ++i) {
        double time = static_cast<double>(i) / sample_rate;
        double value = 2.0 * fmod(time * frequency, 1.0) - 1.0;
        signal.push_back(value);
    }
    return signal;
}



//==================== RYSOWANIE ====================
void plot_signal_cpp(const std::vector<double>& y_data, const std::string& title = "Sygnał z C++") {
    std::vector<double> y_data_copy = y_data;
    std::string title_copy = title;
    std::thread t([y_data_copy, title_copy]() {
        matplot::figure();
        matplot::plot(y_data_copy);
        matplot::title(title_copy);
        matplot::xlabel("Próbki");
        matplot::ylabel("Amplituda");
        matplot::grid(true);
        matplot::show(); 
    });
    t.detach();
}



//==================== DFT ====================
std::vector<std::complex<double>> calculate_dft(const std::vector<double>& input_signal) {
    int N = input_signal.size();
    std::vector<std::complex<double>> output_spectrum(N);
    const std::complex<double> j(0, 1);

    for (int k = 0; k < N; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n) {
            double angle = 2.0 * M_PI * k * n / N;
            sum += input_signal[n] * (std::cos(angle) - j * std::sin(angle));
        }
        output_spectrum[k] = sum;
    }
    return output_spectrum;
}



//==================== IDFT ====================
std::vector<double> calculate_idft(const std::vector<std::complex<double>>& input_spectrum) {
    int N = input_spectrum.size();
    std::vector<double> output_signal(N);
    const std::complex<double> j(0, 1);

    for (int n = 0; n < N; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k) {
            double angle = 2.0 * M_PI * k * n / N;
            sum += input_spectrum[k] * (std::cos(angle) + j * std::sin(angle));
        }
        output_signal[n] = sum.real() / N;
    }
    return output_signal;
}



//==================== PORÓWNANIE DFT vs IDFT ====================
void plot_idft_comparison_cpp(
    const std::vector<double>& time_vector,
    const std::vector<double>& original_signal,
    const std::vector<double>& reconstructed_signal,
    const std::string& title = "Porównanie IDFT z C++") {

    std::vector<double> time_copy = time_vector;
    std::vector<double> original_copy = original_signal;
    std::vector<double> reconstructed_copy = reconstructed_signal;
    std::string title_copy = title;

    std::thread t([time_copy, original_copy, reconstructed_copy, title_copy]() {
        matplot::figure();
        auto p1 = matplot::plot(time_copy, original_copy);
        p1->line_style("-");
        p1->display_name("Oryginalny");

        matplot::hold(matplot::on);
        auto p2 = matplot::plot(time_copy, reconstructed_copy);
        p2->line_style("--");
        p2->display_name("Odtworzony (IDFT)");
        matplot::hold(matplot::off);
        
        matplot::title(title_copy);
        matplot::xlabel("Czas [s]");
        matplot::ylabel("Amplituda");
        matplot::legend();
        matplot::grid(true);
        matplot::show();
    });
    t.detach();
}



//==================== PORÓWNANIE PRZED vs PO FILTRACJI ====================
void plot_filter_comparison_cpp(
    const std::vector<double>& time_vector,
    const std::vector<double>& original_signal,
    const std::vector<double>& filtered_signal,
    const std::string& title = "Porównanie Filtracji 1D z C++") {

    std::vector<double> time_copy = time_vector;
    std::vector<double> original_copy = original_signal;
    std::vector<double> filtered_copy = filtered_signal;
    std::string title_copy = title;

    std::thread t([time_copy, original_copy, filtered_copy, title_copy]() {
        matplot::figure();
        auto p1 = matplot::plot(time_copy, original_copy);
        p1->line_style("-");
        p1->color({0.7f, 0.7f, 0.7f}); // Jasnoszary dla oryginału
        p1->display_name("Oryginalny");
        
        matplot::hold(matplot::on);
        auto p2 = matplot::plot(time_copy, filtered_copy);
        p2->line_style("-");
        p2->color("b");
        p2->display_name("Przefiltrowany");
        matplot::hold(matplot::off);
        
        matplot::title(title_copy);
        matplot::xlabel("Czas [s]");
        matplot::ylabel("Amplituda");
        matplot::legend();
        matplot::grid(true);
        matplot::show();
    });
    t.detach();
}



//==================== FILTRACJA 1D ====================
std::vector<double> convolve_1d(const std::vector<double>& input_signal, const std::vector<double>& kernel) {
    int signal_len = input_signal.size();
    int kernel_len = kernel.size();
    if (signal_len == 0 || kernel_len == 0) {
        return {};
    }
    int output_len = signal_len + kernel_len - 1;

    std::vector<double> output_signal(output_len, 0.0);

    for (int n = 0; n < output_len; ++n) {
        double sum = 0.0;
        for (int k = 0; k < kernel_len; ++k) {
            int signal_index = n - k;
            if (signal_index >= 0 && signal_index < signal_len) {
                sum += input_signal[signal_index] * kernel[k];
            }
        }
        output_signal[n] = sum;
    }
    return output_signal;
}



//==================== FILTRACJA 2D ====================
std::vector<std::vector<double>> convolve_2d(const std::vector<std::vector<double>>& input_image, const std::vector<std::vector<double>>& kernel) {
    if (input_image.empty() || input_image[0].empty() || kernel.empty() || kernel[0].empty()) {
        return {};
    }

    int image_rows = input_image.size();
    int image_cols = input_image[0].size();
    int kernel_rows = kernel.size();
    int kernel_cols = kernel[0].size();

    std::vector<std::vector<double>> output_image(image_rows, std::vector<double>(image_cols, 0.0));
    int kernel_center_row = kernel_rows / 2;
    int kernel_center_col = kernel_cols / 2;

    for (int r = 0; r < image_rows; ++r) {
        for (int c = 0; c < image_cols; ++c) {
            double sum = 0.0;
            for (int kr = 0; kr < kernel_rows; ++kr) {
                for (int kc = 0; kc < kernel_cols; ++kc) {
                    int image_pixel_row = r + (kr - kernel_center_row);
                    int image_pixel_col = c + (kc - kernel_center_col);
                    if (image_pixel_row >= 0 && image_pixel_row < image_rows &&
                        image_pixel_col >= 0 && image_pixel_col < image_cols) {
                        sum += input_image[image_pixel_row][image_pixel_col] * kernel[kr][kc];
                    }
                }
            }
            output_image[r][c] = sum;
        }
    }
    return output_image;
}


//==================== RYSOWANIE 2D ====================
void plot_image_cpp(
    const std::vector<std::vector<double>>& image_data,
    const std::string& title = "Obraz") {

    std::vector<std::vector<double>> image_copy = image_data;
    std::string title_copy = title;

    std::thread t([image_copy, title_copy]() {
        if (image_copy.empty() || image_copy[0].empty()) {
            return;
        }
        matplot::figure();
        matplot::imagesc(image_copy);

        matplot::title(title_copy);
        matplot::xlabel("Kolumny");
        matplot::ylabel("Wiersze");
        matplot::show();
    });
    t.detach();
}

//==================== PYBIND11 ====================
PYBIND11_MODULE(sygnaly, m) {
    m.def("generate_sine", &generate_sine_wave, "Generuje sygnał sinusoidalny", 
          py::arg("frequency"), py::arg("duration"), py::arg("sample_rate"));

    m.def("generate_cosine", &generate_cosine_wave, "Generuje sygnał cosinusoidalny", 
          py::arg("frequency"), py::arg("duration"), py::arg("sample_rate"));

    m.def("generate_square", &generate_square_wave, "Generuje sygnał prostokątny", 
          py::arg("frequency"), py::arg("duration"), py::arg("sample_rate"));

    m.def("generate_sawtooth", &generate_sawtooth_wave, "Generuje sygnał piłokształtny", 
          py::arg("frequency"), py::arg("duration"), py::arg("sample_rate"));

    m.def("plot_signal_cpp", &plot_signal_cpp, "Rysuje wykres pojedynczego sygnału z C++", 
          py::arg("data"), py::arg("title") = "Sygnał z C++");

    m.def("calculate_dft", &calculate_dft, "Oblicza DFT", 
          py::arg("input_signal"));

    m.def("calculate_idft", &calculate_idft, "Oblicza IDFT", 
          py::arg("input_spectrum"));

    m.def("plot_idft_comparison_cpp", &plot_idft_comparison_cpp, "Rysuje porównanie sygnału oryginalnego i odtworzonego po IDFT z C++",
          py::arg("time_vector"), py::arg("original_signal"), py::arg("reconstructed_signal"), py::arg("title") = "Porównanie IDFT z C++");

    m.def("convolve_1d", &convolve_1d, "Oblicza splot 1D dwóch sygnałów",
          py::arg("input_signal"), py::arg("kernel"));
        
    m.def("plot_filter_comparison_cpp", &plot_filter_comparison_cpp, "Rysuje porównanie sygnału przed i po filtracji 1D z C++",
          py::arg("time_vector"), py::arg("original_signal"), py::arg("filtered_signal"), py::arg("title") = "Porównanie Filtracji 1D z C++");

    m.def("convolve_2d", &convolve_2d, "Oblicza splot 2D (filtrację) obrazu z jądrem",
          py::arg("input_image"), py::arg("kernel"));

    m.def("plot_image_cpp", &plot_image_cpp, "Rysuje obraz 2D (macierz) z C++ przy użyciu matplotplusplus",
          py::arg("image_data"), py::arg("title") = "Obraz 2D z C++");
}