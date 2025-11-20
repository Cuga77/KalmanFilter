import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from kalman_models import KalmanFilter, AdaptiveKalmanFilter

def run_experiment(scenario_name, true_omega, model_omega, duration=50.0):
    print(f"Запуск сценария: {scenario_name}")
    dt = 0.01
    times = np.arange(0, duration, dt)
    
    # 1. Генерация данных (Истинная физика)
    # Истинная частота может отличаться от модельной!
    Initial_Amplitude = 5.0
    sigma_meas = 1.0
    
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
    
    # Построение графиков
    fig = make_subplots(
        rows=3, cols=1, 
        shared_xaxes=True,
        vertical_spacing=0.08,
        subplot_titles=(
            "Сценарий 1: Идеальное совпадение (Omega=2.0)", 
            "Сценарий 2: Ошибка модели (True=2.0, Model=1.5)",
            "Адаптация Q (Сценарий 2)"
        )
    )
    
    # Row 1: Ideal
    fig.add_trace(go.Scatter(x=t1, y=true1, name='Истина (Сц.1)', line=dict(color='black', dash='dash')), row=1, col=1)
    fig.add_trace(go.Scatter(x=t1, y=osc1, name='Гармонич. (Сц.1)', line=dict(color='green')), row=1, col=1)
    fig.add_trace(go.Scatter(x=t1, y=adapt1, name='Адаптив. (Сц.1)', line=dict(color='blue', dash='dot')), row=1, col=1)
    
    # Row 2: Mismatch
    fig.add_trace(go.Scatter(x=t2, y=true2, name='Истина (Сц.2)', line=dict(color='black', dash='dash'), showlegend=False), row=2, col=1)
    fig.add_trace(go.Scatter(x=t2, y=osc2, name='Гармонич. (Сц.2)', line=dict(color='green'), showlegend=False), row=2, col=1)
    fig.add_trace(go.Scatter(x=t2, y=adapt2, name='Адаптив. (Сц.2)', line=dict(color='blue', dash='dot'), showlegend=False), row=2, col=1)
    
    # Row 3: Q Adaptation
    fig.add_trace(go.Scatter(x=t2, y=trace2, name='Trace(Q) Адаптив.', line=dict(color='orange')), row=3, col=1)
    
    fig.update_layout(height=1000, title_text="Эксперимент: Устойчивость к ошибкам модели")
    fig.write_html("kalman_experiments.html")
    print("Результаты экспериментов сохранены в 'kalman_experiments.html'")

if __name__ == "__main__":
    main()
