# Zamiast "import sygnaly", importujemy moduł "sygnaly" z pakietu "sygnaly"
# Dajemy mu alias "sgl" dla przejrzystości
from sygnaly import sygnaly as sgl
import matplotlib.pyplot as plt

print("Moduł C++ zaimportowany poprawnie!")
print(f"Zawartość modułu C++: {dir(sgl)}")


# Wygeneruj sygnał za pomocą naszej nowej funkcji z C++
# Używamy teraz aliasu "sgl"
sine_wave = sgl.generate_sine(frequency=1.0, duration=2.0, sample_rate=100)
print("Sygnał został wygenerowany pomyślnie!")

# Narysuj wykres, żeby zobaczyć wynik
plt.figure(figsize=(10, 4))
plt.plot(sine_wave)
plt.title("Sygnał sinusoidalny z C++ (DZIAŁA!)")
plt.xlabel("Próbka")
plt.ylabel("Amplituda")
plt.grid(True)
plt.show()

print("Test zakończony sukcesem!")