import numpy as np
import matplotlib.pyplot as plt

# Параметры сетки (должны совпадать с программой на C)
T = 100.0  # Конечное время
X = 100.0  # Конечное пространство
Nx = 4000  # Число узлов по x
Nt = 5000  # Число узлов по t

# Шаги сетки
tau = T / Nt
h = X / Nx

# Функция аналитического решения
def analytical_solution(x, t):
    if x >= t:
        return np.sin(np.pi * (x - t) / X)  # u(x, t) = sin(pi * (x - t) / 100) для x >= t
    else:
        return np.exp(-(t - x) / T)  # u(x, t) = e^{-(t - x)/100} для x < t

# Чтение численного решения из файла
data = np.loadtxt("output.txt")  # Формат: t, x, u
t_vals = data[:, 0]  # Время
x_vals = data[:, 1]  # Координата x
u_vals = data[:, 2]  # Значения u

# Преобразование численного решения в двумерный массив u(t, x)
u_numerical = np.zeros((Nt + 1, Nx + 1))
for i in range(len(t_vals)):
    k = int(round(t_vals[i] * Nt / T))  # Индекс времени
    m = int(round(x_vals[i] * Nx / X))  # Индекс пространства
    u_numerical[k, m] = u_vals[i]

# Вычисление аналитического решения
u_analytical = np.zeros((Nt + 1, Nx + 1))
for k in range(Nt + 1):
    t = k * tau
    for m in range(Nx + 1):
        x = m * h
        u_analytical[k, m] = analytical_solution(x, t)

# Отладочный вывод для проверки аналитического решения
#print("Проверка аналитического решения:")
#print(f"u_analytical[0][0] (t=0, x=0): {u_analytical[0][0]}")  # Должно быть 0
#print(f"u_analytical[0][2000] (t=0, x=50): {u_analytical[0][2000]}")  # Должно быть 1
#print(f"u_analytical[0][4000] (t=0, x=100): {u_analytical[0][4000]}")  # Должно быть 0
#print(f"u_analytical[100][0] (t=2, x=0): {u_analytical[100][0]}")  # Должно быть e^{-0.02} \approx 0.951229
#print(f"u_analytical[100][2000] (t=2, x=50): {u_analytical[100][2000]}")  # Должно быть sin(0.48pi) \approx 0.649448

# Сравнение решений
error = np.abs(u_numerical - u_analytical)
max_error = np.max(error)
#print(f"\nМаксимальная ошибка между численным и аналитическим решением: {max_error}")
#print(f"Мин u (численное): {np.min(u_numerical)}, Макс u (численное): {np.max(u_numerical)}")
#print(f"Мин u (аналитическое): {np.min(u_analytical)}, Макс u (аналитическое): {np.max(u_analytical)}")

# Построение тепловых карт рядом
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

# Тепловая карта численного решения
im1 = ax1.imshow(u_numerical, extent=[0, X, 0, T], origin='lower', aspect='auto', cmap='viridis')
ax1.set_title('Численное решение')
ax1.set_xlabel('x')
ax1.set_ylabel('t')
plt.colorbar(im1, ax=ax1, label='u(x, t)')

# Тепловая карта аналитического решения
im2 = ax2.imshow(u_analytical, extent=[0, X, 0, T], origin='lower', aspect='auto', cmap='viridis')
ax2.set_title('Аналитическое решение')
ax2.set_xlabel('x')
plt.colorbar(im2, ax=ax2, label='u(x, t)')

plt.tight_layout()
plt.savefig('comp_heatmap.png')

print("Успех! Карты сохранены как 'comp_heatmap.png'.")
