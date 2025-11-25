import numpy as np
import matplotlib.pyplot as plt
from kalman_models import KalmanFilter

def moving_average(data, window_size):
    """Простое скользящее среднее."""
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def exponential_moving_average(data, alpha):
    """Экспоненциальное скользящее среднее (EMA)."""
    ema = [data[0]]
    for x in data[1:]:
        ema.append(alpha * x + (1 - alpha) * ema[-1])
    return np.array(ema)

def run_comparison():
    dt = 0.01
    t_end = 30.0 # Увеличено время
    times = np.arange(0, t_end, dt)
    
    Omega_true = 2.0
    Initial_Amplitude = 5.0
    sigma_meas = 2.5 # Очень сильный шум, чтобы показать преимущество Калмана
    
    np.random.seed(42)
    # Истинный сигнал (скорость, так как мы измеряем скорость)
    # Для наглядности будем сравнивать фильтрацию именно измеряемой величины (скорости)
    true_signal = -Initial_Amplitude * Omega_true * np.sin(Omega_true * times)
    measurements = true_signal + np.random.normal(0, sigma_meas, size=len(times))
    
    # 2. Метод 1: Скользящее среднее
    window = 20 # окно 0.2 сек
    ma_result = moving_average(measurements, window)
    # Корректируем время для MA
    ma_times = times[window-1:]
    
    # 3. Метод 2: EMA
    alpha = 0.1
    ema_result = exponential_moving_average(measurements, alpha)
    
    # 4. Метод 3: Фильтр Калмана (Гармонический)
    # Настраиваем на правильную частоту
    x0 = np.array([[0.0], [0.0]]) # Не знаем старт
    P0 = np.eye(2) * 10.0
    H = np.array([[0.0, 1.0]]) # Измеряем скорость
    R = np.array([[sigma_meas**2]])
    
    c = np.cos(Omega_true * dt)
    s = np.sin(Omega_true * dt)
    F_osc = np.array([
        [c, s/Omega_true],
        [-Omega_true*s, c]
    ])
    Q_osc = np.eye(2) * 1e-4 # Немного даем свободы
    
    kf = KalmanFilter(F_osc, H, Q_osc, R, x0, P0)
    
    kalman_result = []
    for z in measurements:
        z_vec = np.array([[z]])
        kf.predict()
        kf.update(z_vec)
        kalman_result.append(kf.x[1,0]) # Берем СКОРОСТЬ (2-й элемент), чтобы сравнивать с другими методами
        
    kalman_result = np.array(kalman_result)
    
    # --- ВЕСОМОЕ ДОКАЗАТЕЛЬСТВО: RMSE ---
    # Для MA нужно обрезать true_signal
    rmse_ma = np.sqrt(np.mean((true_signal[window-1:] - ma_result)**2))
    rmse_ema = np.sqrt(np.mean((true_signal - ema_result)**2))
    rmse_kalman = np.sqrt(np.mean((true_signal - kalman_result)**2))
    
    print("\n[VERIFICATION] Method Comparison (RMSE):")
    print(f"Moving Average (w={window}): {rmse_ma:.4f}")
    print(f"EMA (alpha={alpha}):        {rmse_ema:.4f}")
    print(f"Kalman Filter:             {rmse_kalman:.4f}")
    print(f"Kalman Improvement over MA:  {rmse_ma - rmse_kalman:.4f}")
    
    # 5. Построение графика
    plt.figure(figsize=(12, 7))
    
    # Истина
    plt.plot(times, true_signal, 'k--', linewidth=2, label='Истинный сигнал')
    
    # Измерения (шум)
    # Используем точки ('.') вместо линии, чтобы показать "облако" измерений
    # alpha=0.3 делает их полупрозрачными, чтобы не перекрывать графики
    plt.plot(times, measurements, '.', color='gray', alpha=0.3, label='Зашумленные измерения')
    
    # Скользящее среднее
    # MA (Moving Average) - берет среднее арифметическое за последние N точек.
    # Плюс: Простота. Минус: Все точки имеют равный вес, поэтому реакция на изменения запаздывает.
    plt.plot(ma_times, ma_result, color='orange', linewidth=2, label=f'Скользящее Среднее (RMSE={rmse_ma:.2f})')
    
    # EMA
    # EMA (Exponential Moving Average) - новые точки имеют больший вес, чем старые.
    # Вес убывает экспоненциально. Это позволяет быстрее реагировать на изменения, чем обычное MA.
    plt.plot(times, ema_result, color='purple', linewidth=2, label=f'Экспоненциальное Сглаживание (RMSE={rmse_ema:.2f})')
    
    # Калман
    # Фильтр Калмана учитывает физику процесса (модель) и статистику шума.
    # Он "предсказывает", где должна быть точка, и корректирует это предсказание измерением.
    plt.plot(times, kalman_result, color='green', linewidth=2, label=f'Фильтр Калмана (RMSE={rmse_kalman:.2f})')
    
    plt.title("Сравнение методов фильтрации: Калман vs Простые методы")
    plt.xlabel("Время (с)")
    plt.ylabel("Скорость (усл. ед.)")
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig("method_comparison.png")
    print("График сравнения методов сохранен в 'method_comparison.png'")
    plt.show() # Отобразить график на экране

if __name__ == "__main__":
    run_comparison()
