import numpy as np
import plotly.graph_objects as go
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
    # 1. Параметры и Генерация данных
    dt = 0.01
    t_end = 20.0
    times = np.arange(0, t_end, dt)
    
    Omega_true = 2.0
    Initial_Amplitude = 5.0
    sigma_meas = 2.0 # Сильный шум, чтобы усложнить задачу
    
    np.random.seed(42)
    # Истинный сигнал (скорость, так как мы измеряем скорость)
    # Для наглядности будем сравнивать фильтрацию именно измеряемой величины (скорости)
    true_signal = -Initial_Amplitude * Omega_true * np.sin(Omega_true * times)
    measurements = true_signal + np.random.normal(0, sigma_meas, size=len(times))
    
    # 2. Метод 1: Скользящее среднее (Moving Average)
    window = 20 # окно 0.2 сек
    ma_result = moving_average(measurements, window)
    # Корректируем время для MA (оно укорачивает массив)
    ma_times = times[window-1:]
    
    # 3. Метод 2: EMA (Low-pass filter)
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
    
    # 5. Построение графика
    fig = go.Figure()
    
    # Истина
    fig.add_trace(go.Scatter(x=times, y=true_signal, name='Истинный сигнал',
                             line=dict(color='black', width=3, dash='dash')))
    
    # Измерения (шум) - делаем полупрозрачными
    fig.add_trace(go.Scatter(x=times, y=measurements, name='Зашумленные измерения',
                             line=dict(color='lightgray', width=1), opacity=0.5))
    
    # Скользящее среднее
    fig.add_trace(go.Scatter(x=ma_times, y=ma_result, name=f'Скользящее среднее (окно {window})',
                             line=dict(color='orange', width=2)))
    
    # EMA
    fig.add_trace(go.Scatter(x=times, y=ema_result, name=f'EMA (alpha={alpha})',
                             line=dict(color='purple', width=2)))
    
    # Калман
    fig.add_trace(go.Scatter(x=times, y=kalman_result, name='Фильтр Калмана',
                             line=dict(color='green', width=3)))
    
    fig.update_layout(
        title="Сравнение методов фильтрации: Калман vs Простые методы",
        xaxis_title="Время (с)",
        yaxis_title="Скорость (усл. ед.)",
        height=700,
        hovermode="x unified"
    )
    
    fig.write_html("method_comparison.html")
    print("График сравнения методов сохранен в 'method_comparison.html'")

if __name__ == "__main__":
    run_comparison()
