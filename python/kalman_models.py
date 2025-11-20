import numpy as np
import matplotlib.pyplot as plt

class KalmanFilter:
    """
    Стандартный дискретный фильтр Калмана.
    
    Состояние: x_k
    Измерения: z_k
    
    Модель:
    x_k = F * x_{k-1} + w_k,  w_k ~ N(0, Q)
    z_k = H * x_k + v_k,      v_k ~ N(0, R)
    """
    def __init__(self, F, H, Q, R, x0, P0):
        self.F = F  # Матрица перехода состояния (State Transition Matrix)
        self.H = H  # Матрица измерений (Measurement Matrix)
        self.Q = Q  # Ковариация шума процесса (Process Noise Covariance)
        self.R = R  # Ковариация шума измерений (Measurement Noise Covariance)
        self.x = x0 # Начальная оценка состояния (Initial State Estimate)
        self.P = P0 # Начальная ковариация ошибки (Initial Error Covariance)
        
    def predict(self):
        # 1. Прогноз состояния: x_{k|k-1} = F * x_{k-1|k-1}
        self.x = self.F @ self.x
        
        # 2. Прогноз ковариации: P_{k|k-1} = F * P_{k-1|k-1} * F^T + Q
        self.P = self.F @ self.P @ self.F.T + self.Q
        
    def update(self, z):
        # 3. Инновация (невязка): y = z - H * x_{k|k-1}
        y = z - self.H @ self.x
        
        # 4. Ковариация инновации: S = H * P_{k|k-1} * H^T + R
        S = self.H @ self.P @ self.H.T + self.R
        
        # 5. Коэффициент усиления Калмана: K = P_{k|k-1} * H^T * S^{-1}
        K = self.P @ self.H.T @ np.linalg.inv(S)
        
        # 6. Обновление состояния: x_{k|k} = x_{k|k-1} + K * y
        self.x = self.x + K @ y
        
        # 7. Обновление ковариации: P_{k|k} = (I - K * H) * P_{k|k-1}
        I = np.eye(self.P.shape[0])
        self.P = (I - K @ self.H) @ self.P
        
        return y # Возвращаем инновацию для адаптивного фильтра

class AdaptiveKalmanFilter(KalmanFilter):
    """
    Адаптивный фильтр Калмана с уточнением ковариации (Таблица 3 в PDF).
    Оценивает матрицу Q на основе последовательности инноваций.
    """
    def __init__(self, F, H, Q, R, x0, P0, window_size=100):
        super().__init__(F, H, Q, R, x0, P0)
        self.window_size = window_size
        self.innovations = [] # Хранение инноваций (y)
        self.K_history = []   # Хранение коэффициентов Калмана (K)
        
    def update(self, z):
        # Стандартное обновление
        y = z - self.H @ self.x
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P @ self.H.T @ np.linalg.inv(S)
        self.x = self.x + K @ y
        I = np.eye(self.P.shape[0])
        self.P = (I - K @ self.H) @ self.P
        
        # Адаптивная часть
        self.innovations.append(y)
        self.K_history.append(K)
        
        if len(self.innovations) > self.window_size:
            self.innovations.pop(0)
            self.K_history.pop(0)
            
        # Оценка Q, если накоплено достаточно данных
        # Формула из PDF для оценки Q (приближенная логика из текста):
        # Q_new = K * Cov(y) * K^T
        if len(self.innovations) == self.window_size:
            # Вычисление выборочной ковариации инноваций
            y_matrix = np.hstack(self.innovations) # (m, N)
            cov_y = np.cov(y_matrix)
            if np.ndim(cov_y) == 0: cov_y = np.array([[cov_y]]) # Обработка скалярного случая
            
            # Используем последний K для обновления (упрощение)
            K_curr = self.K_history[-1]
            
            # Обновление Q
            # Q_est = K * Cov(y) * K.T
            Q_est = K_curr @ cov_y @ K_curr.T
            
            # Смешивание со старым Q для стабильности
            alpha = 0.05
            self.Q = (1 - alpha) * self.Q + alpha * Q_est

def run_simulation():
    # ==========================================
    # 1. ПАРАМЕТРЫ СИСТЕМЫ
    # ==========================================
    dt = 0.01
    t_end = 50.0           # Увеличили время симуляции для наглядности
    Omega_true = 2.0       # Немного уменьшили частоту, чтобы график был читаемее
    Initial_Amplitude = 5.0 # Увеличили амплитуду
    sigma_meas = 1.0       # Увеличили шум, чтобы работа фильтра была заметнее
    
    times = np.arange(0, t_end, dt)
    
    # ==========================================
    # 2. ГЕНЕРАЦИЯ ДАННЫХ (Истинные значения)
    # ==========================================
    # Истинная физика: Гармонический осциллятор
    true_angles = Initial_Amplitude * np.cos(Omega_true * times)
    true_velocities = -Initial_Amplitude * Omega_true * np.sin(Omega_true * times)
    
    # Измерения: Скорость + Шум
    np.random.seed(42)
    measurements = true_velocities + np.random.normal(0, sigma_meas, size=len(times))
    
    # ==========================================
    # 3. НАСТРОЙКА ФИЛЬТРОВ
    # ==========================================
    # Начальное состояние [Угол, Скорость]
    # Допустим, мы не знаем точного начального состояния
    x0 = np.array([[0.0], [0.0]]) 
    P0 = np.eye(2) * 10.0 # Большая начальная неопределенность
    
    # Матрица измерений (измеряем скорость, т.е. 2-ю переменную состояния)
    H = np.array([[0.0, 1.0]])
    
    # Ковариация шума измерений
    R = np.array([[sigma_meas**2]])
    
    # --- МОДЕЛЬ 1: Перманентное вращение (Постоянная скорость) ---
    # x_k = x_{k-1} + dt * v_{k-1}
    # v_k = v_{k-1}
    F_perm = np.array([
        [1.0, dt],
        [0.0, 1.0]
    ])
    # Высокий шум процесса, так как модель неверна
    Q_perm = np.eye(2)
    Q_perm[0,0] = 1e-3
    Q_perm[1,1] = 1.0
    
    kf_perm = KalmanFilter(F_perm, H, Q_perm, R, x0, P0)
    
    # --- МОДЕЛЬ 2: Крутильные колебания (Гармонический осциллятор) ---
    # Точная дискретизация x'' + w^2*x = 0
    c = np.cos(Omega_true * dt)
    s = np.sin(Omega_true * dt)
    F_osc = np.array([
        [c, s/Omega_true],
        [-Omega_true*s, c]
    ])
    # Низкий шум процесса, так как модель верна
    Q_osc = np.eye(2) * 1e-6
    
    kf_osc = KalmanFilter(F_osc, H, Q_osc, R, x0, P0)
    
    # --- МОДЕЛЬ 3: Адаптивный фильтр (изначально как Модель 2) ---
    # Начинаем с немного неверного Q, чтобы увидеть адаптацию
    Q_adapt_init = np.eye(2) * 1e-5 
    kf_adapt = AdaptiveKalmanFilter(F_osc, H, Q_adapt_init, R, x0, P0)
    
    # ==========================================
    # 4. ЦИКЛ СИМУЛЯЦИИ
    # ==========================================
    est_perm = []
    est_osc = []
    est_adapt = []
    
    for z in measurements:
        z_vec = np.array([[z]])
        
        # Прогноз (Predict)
        kf_perm.predict()
        kf_osc.predict()
        kf_adapt.predict()
        
        # Обновление (Update)
        kf_perm.update(z_vec)
        kf_osc.update(z_vec)
        kf_adapt.update(z_vec)
        
        # Сохраняем оценку угла (1-я переменная состояния)
        est_perm.append(kf_perm.x[0,0])
        est_osc.append(kf_osc.x[0,0])
        est_adapt.append(kf_adapt.x[0,0])
        
    # ==========================================
    # 5. АНАЛИЗ РЕЗУЛЬТАТОВ
    # ==========================================
    est_perm = np.array(est_perm)
    est_osc = np.array(est_osc)
    est_adapt = np.array(est_adapt)
    
    mse_perm = np.mean((true_angles - est_perm)**2)
    mse_osc = np.mean((true_angles - est_osc)**2)
    mse_adapt = np.mean((true_angles - est_adapt)**2)
    
    print(f"MSE Модель 1 (Пост. скорость): {mse_perm:.6f}")
    print(f"MSE Модель 2 (Гармоническая):  {mse_osc:.6f}")
    print(f"MSE Модель 3 (Адаптивная):     {mse_adapt:.6f}")
    
    # ==========================================
    # 6. ПОСТРОЕНИЕ ИНТЕРАКТИВНОГО ГРАФИКА (PLOTLY)
    # ==========================================
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, 
                            vertical_spacing=0.1,
                            subplot_titles=("Сравнение оценки угла", "Ошибка оценки"))

        # График 1: Углы
        fig.add_trace(go.Scatter(x=times, y=true_angles, name='Истинный угол',
                                 line=dict(color='black', width=2, dash='dash')), row=1, col=1)
        fig.add_trace(go.Scatter(x=times, y=est_perm, name='Модель 1 (Пост. скорость)',
                                 line=dict(color='red', width=1)), row=1, col=1)
        fig.add_trace(go.Scatter(x=times, y=est_osc, name='Модель 2 (Гармоническая)',
                                 line=dict(color='green', width=1)), row=1, col=1)
        fig.add_trace(go.Scatter(x=times, y=est_adapt, name='Модель 3 (Адаптивная)',
                                 line=dict(color='blue', width=2, dash='dot')), row=1, col=1)

        # График 2: Ошибки
        fig.add_trace(go.Scatter(x=times, y=true_angles - est_perm, name='Ошибка Модель 1',
                                 line=dict(color='red', width=1), showlegend=False), row=2, col=1)
        fig.add_trace(go.Scatter(x=times, y=true_angles - est_osc, name='Ошибка Модель 2',
                                 line=dict(color='green', width=1), showlegend=False), row=2, col=1)
        fig.add_trace(go.Scatter(x=times, y=true_angles - est_adapt, name='Ошибка Модель 3',
                                 line=dict(color='blue', width=2, dash='dot'), showlegend=False), row=2, col=1)

        fig.update_layout(height=800, title_text="Результаты Фильтра Калмана (Интерактивный график)",
                          hovermode="x unified")
        fig.update_xaxes(title_text="Время (с)", row=2, col=1)
        fig.update_yaxes(title_text="Угол (рад)", row=1, col=1)
        fig.update_yaxes(title_text="Ошибка (рад)", row=2, col=1)

        fig.write_html("kalman_results.html")
        print("Интерактивный график сохранен в 'kalman_results.html'")
        
    except ImportError:
        print("Plotly не установлен. Используется Matplotlib.")
        # Fallback to Matplotlib if Plotly fails (though we installed it)
        plt.figure(figsize=(12, 6))
        plt.plot(times, true_angles, 'k--', label='Истинный угол')
        plt.plot(times, est_perm, 'r-', label='Модель 1')
        plt.plot(times, est_osc, 'g-', label='Модель 2')
        plt.plot(times, est_adapt, 'b:', label='Модель 3')
        plt.legend()
        plt.grid(True)
        plt.savefig('kalman_results.png')
        print("График сохранен в 'kalman_results.png'")

if __name__ == "__main__":
    run_simulation()
