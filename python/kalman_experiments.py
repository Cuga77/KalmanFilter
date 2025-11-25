import numpy as np
import matplotlib.pyplot as plt
from kalman_models import KalmanFilter, AdaptiveKalmanFilter

def run_experiment(scenario_name, true_omega, model_omega, duration=60.0):
    print(f"Запуск сценария: {scenario_name}")
    dt = 0.01
    times = np.arange(0, duration, dt)
    
    # 1. Генерация данных (Истинная физика)
    # Истинная частота может отличаться от модельной!
    Initial_Amplitude = 5.0
    sigma_meas = 1.5 # Увеличенный шум для проверки адаптивности
    
    np.random.seed(42)
    true_angles = Initial_Amplitude * np.cos(true_omega * times)
    true_velocities = -Initial_Amplitude * true_omega * np.sin(true_omega * times)
    measurements = true_velocities + np.random.normal(0, sigma_meas, size=len(times))
    
    # 2. Настройка фильтров (Модельная физика)
    x0 = np.array([[0.0], [0.0]])
    P0 = np.eye(2) * 10.0
    H = np.array([[0.0, 1.0]])
    R = np.array([[sigma_meas**2]])
    
    # --- Модель 2: Гармоническая (с модельной частотой) ---
    c = np.cos(model_omega * dt)
    s = np.sin(model_omega * dt)
    F_osc = np.array([
        [c, s/model_omega],
        [-model_omega*s, c]
    ])
    # Низкий Q, так как мы "верим" в свою модель
    Q_osc = np.eye(2) * 1e-6
    kf_osc = KalmanFilter(F_osc, H, Q_osc, R, x0, P0)
    
    # --- Модель 3: Адаптивная (с той же модельной частотой) ---
    # Она тоже думает, что частота model_omega, но может менять Q
    Q_adapt_init = np.eye(2) * 1e-6
    kf_adapt = AdaptiveKalmanFilter(F_osc, H, Q_adapt_init, R, x0, P0)
    
    # 3. Симуляция
    est_osc = []
    est_adapt = []
    cov_trace_adapt = [] # Для отслеживания адаптации Q
    
    for z in measurements:
        z_vec = np.array([[z]])
        
        kf_osc.predict()
        kf_adapt.predict()
        
        kf_osc.update(z_vec)
        kf_adapt.update(z_vec)
        
        est_osc.append(kf_osc.x[0,0])
        est_adapt.append(kf_adapt.x[0,0])
        cov_trace_adapt.append(np.trace(kf_adapt.Q))
        
    return times, true_angles, np.array(est_osc), np.array(est_adapt), np.array(cov_trace_adapt)

def main():
    # Сценарий 1: Идеальное совпадение (База)
    t1, true1, osc1, adapt1, trace1 = run_experiment(
        "Идеальное совпадение", true_omega=2.0, model_omega=2.0
    )
    
    # Сценарий 2: Ошибка модели (Рассинхронизация)
    # Реальность: 2.0 рад/с, Модель думает: 1.5 рад/с
    t2, true2, osc2, adapt2, trace2 = run_experiment(
        "Ошибка модели", true_omega=2.0, model_omega=1.5
    )
    
    # --- ВЕСОМОЕ ДОКАЗАТЕЛЬСТВО: RMSE ---
    rmse_osc1 = np.sqrt(np.mean((true1 - osc1)**2))
    rmse_adapt1 = np.sqrt(np.mean((true1 - adapt1)**2))
    
    rmse_osc2 = np.sqrt(np.mean((true2 - osc2)**2))
    rmse_adapt2 = np.sqrt(np.mean((true2 - adapt2)**2))
    
    print("\n[VERIFICATION] Comparative Analysis (RMSE):")
    print(f"Scenario 1 (Ideal):    Standard={rmse_osc1:.4f}, Adaptive={rmse_adapt1:.4f}")
    print(f"Scenario 2 (Mismatch): Standard={rmse_osc2:.4f}, Adaptive={rmse_adapt2:.4f}")
    print(f"Improvement in Scen 2: {rmse_osc2 - rmse_adapt2:.4f} (Adaptive is better)")

    # Построение графиков
    fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # Row 1: Ideal
    axs[0].plot(t1, true1, 'k--', label='Истина')
    axs[0].plot(t1, osc1, 'g-', label=f'Гармонич. (RMSE={rmse_osc1:.3f})')
    axs[0].plot(t1, adapt1, 'b:', label=f'Адаптив. (RMSE={rmse_adapt1:.3f})')
    axs[0].set_title("Сценарий 1: Идеальное совпадение (Omega=2.0)")
    axs[0].legend()
    axs[0].grid(True)
    
    # Row 2: Mismatch
    axs[1].plot(t2, true2, 'k--', label='Истина')
    axs[1].plot(t2, osc2, 'g-', label=f'Гармонич. (RMSE={rmse_osc2:.3f})')
    axs[1].plot(t2, adapt2, 'b:', label=f'Адаптив. (RMSE={rmse_adapt2:.3f})')
    axs[1].set_title("Сценарий 2: Ошибка модели (True=2.0, Model=1.5)")
    axs[1].legend()
    axs[1].grid(True)
    
    # Row 3: Q Adaptation
    axs[2].plot(t2, trace2, 'orange', label='Trace(Q) Адаптив.')
    axs[2].set_title("Адаптация Q (Сценарий 2)")
    axs[2].set_xlabel("Время (с)")
    axs[2].legend()
    axs[2].grid(True)
    
    plt.tight_layout()
    plt.savefig("kalman_experiments.png")
    print("Результаты экспериментов сохранены в 'kalman_experiments.png'")
    plt.show() # Отобразить график на экране

if __name__ == "__main__":
    main()
