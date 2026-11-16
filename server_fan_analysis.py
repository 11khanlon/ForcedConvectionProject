#Preliminary 2U server fan sizing based on user-provided airflow graphs
#4 fans side-by-side --> parallel flow


import numpy as np
import matplotlib.pyplot as plt

# -------- USER INPUTS --------
T_IN_C = 35.0
ALTITUDE_M = 3000.0
AIR_RISE_K = 15.0
CP_AIR = 1005.0
N_FANS = 4
CPU_COUNT, CPU_W = 2, 600.0  # ASSUMPTION: earlier 600 W per CPU
GPU_COUNT, GPU_W = 4, 1000.0  # ASSUMPTION: truncated statement means each GPU
EXPANSION_W = 75.0
OTHER_W = 0.0  # Add memory, drives, fans, board, PSU losses as known
# Alternative package case: set CPU_W=190 or 250 if applicable to same server.

# Graph 2: APPROXIMATE total system airflow, not per-fan airflow.
RPM = np.array([3400., 6800., 10400.])
AIR_COOLED_CFM = np.array([85., 160., 240.])
WATER_COOLED_CFM = np.array([60., 120., 185.])

# Graph 1: APPROXIMATE local velocities (m/s), read visually.
PWM = np.array([30., 45., 60., 75., 100.])
LOCAL_VELOCITY_M_S = {
    'CPU-A': [1.8, 3.1, 4.4, 5.8, 6.6],
    'CPU-B': [2.7, 4.1, 5.5, 6.6, 7.4],
    'GPU-A1': [2.4, 3.8, 5.0, 6.2, 6.8],
    'GPU-A2': [1.9, 2.8, 3.8, 4.8, 5.5],
    'GPU-B1': [0.9, 1.8, 2.7, 3.5, 4.2],
    'PS': [0.9, 1.2, 1.4, 1.8, 1.6],
}

# Optional hydraulic estimate: these are NOT supplied by the assignment.
# Enter estimated flow area and loss coefficient for each SERIAL segment.
# For parallel branches, model each branch separately; don't sum them as serial.
SERIAL_LOSSES = {
    # 'inlet': {'area_m2': 0.02, 'K': 2.0},
    # 'exhaust': {'area_m2': 0.02, 'K': 1.5},
}

CFM_PER_M3_S = 2118.880003

def isa_pressure_pa(altitude_m):
    """ISA troposphere approximation, valid for 0-11 km altitude."""
    return 101325.0 * (1.0 - 0.0065 * altitude_m / 288.15) ** 5.25588

def air_density_kg_m3(altitude_m, temp_c):
    return isa_pressure_pa(altitude_m) / (287.05 * (temp_c + 273.15))

def required_airflow(heat_w, rise_k, rho):
    mdot = heat_w / (CP_AIR * rise_k)
    return mdot, mdot / rho * CFM_PER_M3_S

def estimated_serial_pressure_drop_pa(cfm, rho):
    if not SERIAL_LOSSES:
        return None
    vdot = cfm / CFM_PER_M3_S
    return sum(0.5 * rho * (vdot / s['area_m2']) ** 2 * s['K']
               for s in SERIAL_LOSSES.values())

def main():
    rho = air_density_kg_m3(ALTITUDE_M, T_IN_C)
    heat_w = CPU_COUNT * CPU_W + GPU_COUNT * GPU_W + EXPANSION_W + OTHER_W
    mdot, required_cfm = required_airflow(heat_w, AIR_RISE_K, rho)
    max_cfm = AIR_COOLED_CFM[-1]
    max_heat_removed_w = (max_cfm / CFM_PER_M3_S) * rho * CP_AIR * AIR_RISE_K
    # Simple fan failure: 3/4 flow only if each fan maintains same operating flow.
    failure_cfm = max_cfm * (N_FANS - 1) / N_FANS
    failure_heat_removed_w = (failure_cfm / CFM_PER_M3_S) * rho * CP_AIR * AIR_RISE_K
    dp = estimated_serial_pressure_drop_pa(max_cfm, rho)

    print('--- 2U SERVER AIRFLOW ANALYSIS ---')
    print(f'Inlet temperature: {T_IN_C:.1f} C')
    print(f'Altitude: {ALTITUDE_M:.0f} m')
    print(f'Atmospheric pressure (ISA): {isa_pressure_pa(ALTITUDE_M)/1000:.2f} kPa')
    print(f'Inlet air density: {rho:.3f} kg/m^3')
    print(f'Total assumed air-cooled heat load: {heat_w:.0f} W')
    print(f'Required mass flow: {mdot:.3f} kg/s')
    print(f'Required total airflow: {required_cfm:.1f} CFM')
    print(f'Max reported air-cooled system airflow: {max_cfm:.0f} CFM')
    print(f'Air-side heat removal at max flow, {AIR_RISE_K:g} K rise: {max_heat_removed_w:.0f} W')
    print(f'Available / required airflow: {max_cfm/required_cfm:.1%}')
    print(f'Approx. one-fan-failure flow: {failure_cfm:.0f} CFM (simplifying assumption)')
    print(f'Approx. failure heat removal at {AIR_RISE_K:g} K rise: {failure_heat_removed_w:.0f} W')
    print('Estimated serial pressure drop:',
          f'{dp:.1f} Pa' if dp is not None else 'Not available: K and area inputs missing')
    print('CAUTION: fan pressure-flow curve required to verify operating point.')
    print('CAUTION: check if graph 2 is total-system flow at applicable altitude.')

    # Plot 1: Required airflow versus heat load, compared with measured system flow.
    powers = np.linspace(0, 6500, 250)
    required = np.array([required_airflow(p, AIR_RISE_K, rho)[1] for p in powers])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(powers / 1000, required, label='Required flow (15 C air rise)')
    ax.axhline(max_cfm, linestyle='--', label='Reported max air-cooled flow (~240 CFM)')
    ax.scatter([heat_w / 1000], [required_cfm], zorder=4, label='Selected load')
    ax.set(xlabel='Total heat transferred to air (kW)', ylabel='Total airflow (CFM)',
           title='Required vs. reported server airflow at 3,000 m')
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig('required_vs_available_airflow.png', dpi=180)

    # Plot 2: approximate graph 2 data.
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(RPM, AIR_COOLED_CFM, 'o-', label='Air-cooled configuration')
    ax.plot(RPM, WATER_COOLED_CFM, 'o-', label='Water-cooled configuration')
    ax.set(xlabel='Fan speed (rpm)', ylabel='Total system airflow (CFM)',
           title='Approximate digitization of provided fan-speed graph')
    ax.grid(alpha=.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig('airflow_vs_fan_speed.png', dpi=180)

    # Plot 3: local velocities from graph 1.
    fig, ax = plt.subplots(figsize=(8, 5))
    for name, velocities in LOCAL_VELOCITY_M_S.items():
        ax.plot(PWM, velocities, 'o-', label=name)
    ax.set(xlabel='Main fan PWM (%)', ylabel='Local air velocity (m/s)',
           title='Approximate digitization of local velocity graph')
    ax.grid(alpha=.3)
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig('local_air_velocity_vs_pwm.png', dpi=180)
    plt.show()

if __name__ == '__main__':
    main()
