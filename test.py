import time
from sygnaly import sygnaly as sgl
import numpy as np

frequency_1d = 2.0
duration_1d = 2.0 
sample_rate_1d = 100
num_samples_1d = int(duration_1d * sample_rate_1d)
time_vector_1d_list = np.linspace(0, duration_1d, num_samples_1d, endpoint=False).tolist() 

def wait_for_plot(seconds=4, message="Następny wykres za"):
    print(f"\n{message} {seconds}s...")
    time.sleep(seconds)

sine_wave = sgl.generate_sine(frequency_1d, duration_1d, sample_rate_1d)
sgl.plot_signal_cpp(sine_wave, f"Sygnał Sinusoidalny ({frequency_1d} Hz)")
wait_for_plot()

cosine_wave = sgl.generate_cosine(frequency_1d, duration_1d, sample_rate_1d)
sgl.plot_signal_cpp(cosine_wave, f"Sygnał Cosinusoidalny ({frequency_1d} Hz)")
wait_for_plot()

square_wave = sgl.generate_square(frequency_1d, duration_1d, sample_rate_1d)
sgl.plot_signal_cpp(square_wave, f"Sygnał Prostokątny ({frequency_1d} Hz)")
wait_for_plot()

sawtooth_wave = sgl.generate_sawtooth(frequency_1d, duration_1d, sample_rate_1d)
sgl.plot_signal_cpp(sawtooth_wave, f"Sygnał Piłokształtny ({frequency_1d} Hz)")
wait_for_plot()

signal_for_dft_idft = sine_wave
spectrum_complex = sgl.calculate_dft(signal_for_dft_idft)
reconstructed_signal = sgl.calculate_idft(spectrum_complex)

sgl.plot_idft_comparison_cpp(
    time_vector_1d_list,
    signal_for_dft_idft,    
    reconstructed_signal,        
    f"DFT vs IDFT (sin {frequency_1d}Hz)"
)
wait_for_plot()

noise_amplitude = 0.4
noise = noise_amplitude * np.random.normal(0, 1, len(sine_wave))
noisy_sine_wave_list = [s + n for s, n in zip(sine_wave, noise)] 

kernel_size_1d = 9
kernel_moving_avg = [1.0/kernel_size_1d] * kernel_size_1d 

full_filtered_signal_1d = sgl.convolve_1d(noisy_sine_wave_list, kernel_moving_avg)

offset_1d = (kernel_size_1d - 1) // 2
filtered_signal_trimmed_1d = full_filtered_signal_1d[offset_1d : offset_1d + num_samples_1d]

sgl.plot_filter_comparison_cpp(
    time_vector_1d_list, 
    noisy_sine_wave_list,
    filtered_signal_trimmed_1d,
    f"Filtracja 1D"
)
wait_for_plot()

size_2d = 60 
duration_2d = 1.0 
sample_rate_2d = int(size_2d / duration_2d)
freq_x_2d = 2.5 
freq_y_2d = 1.5 

signal_x_comp = sgl.generate_sine(freq_x_2d, duration_2d, sample_rate_2d)
signal_y_comp = sgl.generate_cosine(freq_y_2d, duration_2d, sample_rate_2d)

signal_x_comp = (signal_x_comp + [0]*size_2d)[:size_2d]
signal_y_comp = (signal_y_comp + [0]*size_2d)[:size_2d]

sample_image_2d = [[(signal_x_comp[i] + signal_y_comp[j]) / 2.0 for j in range(size_2d)] for i in range(size_2d)]

sgl.plot_image_cpp(sample_image_2d, "Oryginalny (sin+cos)")
wait_for_plot()

kernel_blur_size_2d = 5 
kernel_blur_value = 1.0 / (kernel_blur_size_2d * kernel_blur_size_2d)
blur_kernel_2d = [[kernel_blur_value for _ in range(kernel_blur_size_2d)] for _ in range(kernel_blur_size_2d)]

filtered_image_blur = sgl.convolve_2d(sample_image_2d, blur_kernel_2d)
sgl.plot_image_cpp(filtered_image_blur, f"Filtracja 2D")
wait_for_plot(15)