from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from air import required_airflow


def make_plots(config, props, heat_w, cfm_req, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    powers = np.linspace(0, max(6500, heat_w*1.15), 200)
    req = [required_airflow(p, config.AIR_RISE_K, props)[2] for p in powers]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(powers/1000, req, label='Required inlet airflow')
    ax.axhline(config.AIR_COOLED_CFM[-1], ls='--', label='Reported air-cooled max (~240 CFM)')
    ax.scatter([heat_w/1000], [cfm_req], label='Selected heat load', zorder=5)
    ax.set(xlabel='Heat transferred to air (kW)', ylabel='Inlet airflow (CFM)',
           title='Required airflow vs. reported flow at 3000 m')
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(output_dir/'required_vs_available_airflow.png', dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(config.RPM, config.AIR_COOLED_CFM, 'o-', label='Air-cooled')
    ax.plot(config.RPM, config.WATER_COOLED_CFM, 'o-', label='Water-cooled')
    ax.set(xlabel='Fan speed (RPM)', ylabel='Reported whole-server airflow (CFM)',
           title='Approximate readings: airflow vs. fan speed')
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(output_dir/'airflow_vs_fan_speed.png', dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    for name, values in config.LOCAL_VELOCITY_M_S.items():
        ax.plot(config.PWM, values, marker='o', lw=1.3, label=name)
    ax.set(xlabel='Main fan PWM (%)', ylabel='Local velocity (m/s)',
           title='Approximate local velocity readings')
    ax.grid(alpha=.3); ax.legend(ncol=3, fontsize=8); fig.tight_layout()
    fig.savefig(output_dir/'local_air_velocity_vs_pwm.png', dpi=180); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(config.PWM, config.PWM_TOTAL_AIRFLOW_CFM, 'o-', label='Right-axis whole-server CFM')
    ax.set(xlabel='Main fan PWM (%)', ylabel='Total airflow (CFM)',
           title='Approximate total airflow vs. PWM')
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(output_dir/'total_airflow_vs_pwm.png', dpi=180); plt.close(fig)

    if config.FAN_CURVE is not None and (config.SERIAL_LOSSES or config.PARALLEL_BRANCHES):
        from hydraulics import fan_bank_pressure_pa, network_pressure_drop_pa, operating_point
        curve = config.FAN_CURVE
        lo, hi = config.N_FANS*curve['cfm'][0], config.N_FANS*curve['cfm'][-1]
        q = np.linspace(lo, hi, 80)
        pf = [fan_bank_pressure_pa(x, curve, props, config.N_FANS) for x in q]
        ps = [network_pressure_drop_pa(x, config.SERIAL_LOSSES, config.PARALLEL_BRANCHES, props)[0] for x in q]
        op = operating_point(curve, config.SERIAL_LOSSES, config.PARALLEL_BRANCHES, props, config.N_FANS)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(q, pf, label='Fan bank'); ax.plot(q, ps, label='System')
        ax.scatter([op['cfm']], [op['pressure_pa']], label='Operating point', zorder=5)
        ax.set(xlabel='Total airflow (CFM)', ylabel='Static pressure (Pa)',
               title='Fan bank and system resistance (specified curve)')
        ax.grid(alpha=.3); ax.legend(); fig.tight_layout()
        fig.savefig(output_dir/'fan_system_operating_point.png', dpi=180); plt.close(fig)
