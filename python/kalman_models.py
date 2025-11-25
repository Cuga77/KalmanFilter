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
    t_end = 100.0            # Время симуляции (увеличено для наглядности)
    Omega_true = 2.0        
    Initial_Amplitude = 10.0 
    sigma_meas = 0.8        # Шум измерений (увеличен)
    
    times = np.arange(0, t_end, dt)
    
    # ==========================================
    # 2. ГЕНЕРАЦИЯ ДАННЫХ
    # ==========================================
    true_angles = Initial_Amplitude * np.cos(Omega_true * times)
    true_velocities = -Initial_Amplitude * Omega_true * np.sin(Omega_true * times)
    
    np.random.seed(42)
    measurements = true_velocities + np.random.normal(0, sigma_meas, size=len(times))
    
    # ==========================================
    # 3. НАСТРОЙКА ФИЛЬТРОВ
    # ==========================================
    x0 = np.array([[0.0], [0.0]]) 
    P0 = np.eye(2) * 10.0 
    H = np.array([[0.0, 1.0]])
    R = np.array([[sigma_meas**2]])
    
    # --- МОДЕЛЬ 2 (Гармоническая - "Эталон") ---
    c = np.cos(Omega_true * dt)
    s = np.sin(Omega_true * dt)
    F_osc = np.array([
        [c, s/Omega_true],
        [-Omega_true*s, c]
    ])
    # Маленький шум процесса (мы доверяем модели)
    Q_osc = np.eye(2) * 1e-5
    
    kf_osc = KalmanFilter(F_osc, H, Q_osc, R, x0, P0)
    
    # ==========================================
    # 4. ЦИКЛ СИМУЛЯЦИИ
    # ==========================================
    est_osc = []
    P_history_osc = [] # <--- НОВОЕ: Храним диагональ ковариации
    
    for z in measurements:
        z_vec = np.array([[z]])
        
        kf_osc.predict()
        kf_osc.update(z_vec)
        
        est_osc.append(kf_osc.x[0,0])
        
        # Сохраняем дисперсию ошибки угла (элемент P[0,0])
        # P показывает теоретический квадрат ошибки оценки
        P_history_osc.append(kf_osc.P[0,0])
        
    est_osc = np.array(est_osc)
    
    # ==========================================
    # 5. ПОСТРОЕНИЕ ГРАФИКА С 3-SIGMA (PLOTLY)
    # ==========================================
    # ==========================================
    # 5. ПОСТРОЕНИЕ ГРАФИКА С 3-SIGMA (MATPLOTLIB)
    # ==========================================
    
    # Расчет 3-сигма коридора
    # ПРАВИЛО 3-Х СИГМ (3-Sigma Rule):
    # В нормальном распределении (гауссиане):
    # - 68.2% значений лежат в пределах ±1σ (одна сигма)
    # - 95.4% значений лежат в пределах ±2σ
    # - 99.7% значений лежат в пределах ±3σ
    #
    # Если наш фильтр работает правильно, то истинная ошибка (красная линия)
    # должна находиться внутри коридора ±3σ (серая зона) в 99.7% случаев.
    # Если ошибка часто вылетает за этот коридор — значит, фильтр "слишком самоуверен"
    # (думает, что он точнее, чем есть на самом деле) или модель неверна.
    
    sigma_angle = np.sqrt(np.array(P_history_osc))
    upper_bound = 3 * sigma_angle
    lower_bound = -3 * sigma_angle
    
    error = true_angles - est_osc
    
    # --- ВЕСОМОЕ ДОКАЗАТЕЛЬСТВО: RMSE ---
    rmse = np.sqrt(np.mean(error**2))
    print(f"\n[VERIFICATION] Model 2 (Harmonic) RMSE: {rmse:.4f}")
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    # График 1: Углы
    ax1.plot(times, true_angles, 'k--', label='Истинный угол')
    ax1.plot(times, est_osc, 'g-', label='Оценка Калмана')
    ax1.set_title("Оценка угла (Model 2)")
    ax1.set_ylabel("Угол (рад)")
    ax1.legend()
    ax1.grid(True)
    
    # График 2: Ошибки + Коридор
    ax2.plot(times, upper_bound, 'gray', alpha=0.3)
    ax2.plot(times, lower_bound, 'gray', alpha=0.3)
    ax2.fill_between(times, lower_bound, upper_bound, color='gray', alpha=0.2, label='Коридор 3σ')
    ax2.plot(times, error, 'r-', label='Ошибка (True - Est)')
    
    ax2.set_title(f"Ошибка оценки и 3σ коридор (RMSE = {rmse:.4f})")
    ax2.set_xlabel("Время (с)")
    ax2.set_ylabel("Ошибка (рад)")
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig("kalman_3sigma.png")
    print("График сохранен в 'kalman_3sigma.png'")
    plt.show() # Отобразить график на экране

if __name__ == "__main__":
    run_simulation()
