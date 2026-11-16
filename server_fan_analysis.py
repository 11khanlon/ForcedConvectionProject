"""Preliminary 2U server fan/thermal-fluid analysis.
Install: pip install numpy matplotlib CoolProp
All graph values are approximate visual readings, NOT manufacturer specifications.
No fan static-pressure curve supplied: cannot solve actual fan operating point.
"""
import numpy as np
import matplotlib.pyplot as plt

try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

# ---------------- INPUTS ----------------
T_IN_C = 35.0
ALTITUDE_M = 3000.0
AIR_RISE_K = 15.0
N_FANS = 4
CPU_COUNT, CPU_W = 2, 600.0  # verify 600 W is PER CPU
GPU_COUNT, GPU_W = 4, 1000.0  # verify original statement
EXPANSION_W = 75.0
OTHER_W = 0.0

# Graph 2: approximate WHOLE-SERVER flow readings, not flow per fan.
RPM = np.array([3400., 6800., 10400.])
AIR_COOLED_CFM = np.array([85., 160., 240.])
WATER_COOLED_CFM = np.array([60., 120., 185.])

# Graph 1: approximate local air velocities from LEFT axis, m/s.
PWM = np.array([30., 45., 60., 75., 100.])
LOCAL_VELOCITY_M_S = {
    'CPU-A':       [1.8, 3.1, 4.4, 5.8, 6.6],
    'CPU-B':       [2.7, 4.1, 5.5, 6.6, 7.4],
    'GPU-A1':      [2.4, 3.8, 5.0, 6.2, 6.8],
    'GPU-A2':      [1.9, 2.8, 3.8, 4.8, 5.5],
    'GPU-B1':      [0.9, 1.8, 2.7, 3.5, 4.2],
    'GPU-B2':      [1.9, 2.7, 3.7, 4.7, 5.5],
    'PCI-1':       [2.0, 3.2, 4.6, 5.2, 6.1],
    'PCI-2':       [0.5, 0.6, 0.8, 1.1, 1.4],
    'PCI-4-Exit':  [2.5, 4.0, 5.4, 6.6, 7.5],
    'BXE':         [3.2, 4.6, 5.9, 7.0, 7.8],
    'PS':          [0.9, 1.2, 1.4, 1.8, 1.6],
}
# Graph 1 PURPLE dotted curve, read against RIGHT axis in CFM.
# These values are especially approximate; replace with digitized originals.
PWM_TOTAL_AIRFLOW_CFM = np.array([45., 135., 180., 230., 255.])

# OPTIONAL flow-loss geometry. Enter only real or explicitly assumed values.
# Each serial section uses its own cross-sectional area and K.
SERIAL_LOSSES = {
    # 'inlet': {'area_m2': 0.020, 'K': 0.5},  # example only
    # 'outlet': {'area_m2': 0.020, 'K': 1.0}, # example only
}
# Optional fin-channel geometry, e.g. s=fin gap, H=channel height, L=flow length.
# Leave None until actual heatsink dimensions are supplied.
HEATSINK_CHANNEL = None
# Example (NOT real geometry):
# HEATSINK_CHANNEL = {'gap_m': 0.002, 'height_m': 0.025, 'length_m': 0.10,
#                     'K_entrance': 0.5, 'K_exit': 1.0}
# If using measured local velocities, only evaluate Nu and pressure drop for
# a channel if velocity represents the CHANNEL MEAN velocity.

CFM_PER_M3_S = 2118.880003
R_AIR = 287.05


def isa_pressure_pa(altitude_m):
    """ISA troposphere pressure, 0 <= altitude <= 11000 m."""
    if not 0 <= altitude_m <= 11000:
        raise ValueError('ISA troposphere model supports 0 to 11 km')
    return 101325.0 * (1.0 - 0.0065 * altitude_m / 288.15) ** 5.25588


def air_properties(temp_c, pressure_pa):
    T = temp_c + 273.15
    if PropsSI is not None:
        rho = PropsSI('Dmass', 'T', T, 'P', pressure_pa, 'Air')
        cp = PropsSI('Cpmass', 'T', T, 'P', pressure_pa, 'Air')
        mu = PropsSI('VISCOSITY', 'T', T, 'P', pressure_pa, 'Air')
        k = PropsSI('CONDUCTIVITY', 'T', T, 'P', pressure_pa, 'Air')
        source = 'CoolProp'
    else:
        rho = pressure_pa / (R_AIR * T)
        cp = 1005.0
        mu = 1.9e-5  # approximate near 35-50 C
        k = 0.0275   # approximate near 35-50 C
        source = 'ideal gas / approximate constant properties'
    return {'rho': rho, 'cp': cp, 'mu': mu, 'k': k,
            'Pr': cp * mu / k, 'source': source}


def required_airflow(heat_w, rise_k, props):
    mdot = heat_w / (props['cp'] * rise_k)
    vdot = mdot / props['rho']
    return mdot, vdot, vdot * CFM_PER_M3_S

#for steady, incompressible flow between two points in a duct, the pressure drop is given by the Darcy-Weisbach equation:
def serial_pressure_drop_pa(cfm, rho):
    if not SERIAL_LOSSES:
        return None
    vdot = cfm / CFM_PER_M3_S
    return sum(0.5 * rho * (vdot / item['area_m2']) ** 2 * item['K']
               for item in SERIAL_LOSSES.values())


def heatsink_channel_analysis(velocity_m_s, props, geometry):
    """Rectangular channel estimate; assumes supplied velocity is channel mean.
    Gnielinski applied only for Re>=3000, subject to smooth fully-developed
    channel assumptions; laminar rectangular-duct Nu requires another model.
    """
    gap = geometry['gap_m']
    height = geometry['height_m']
    L = geometry['length_m']
    dh = 2.0 * gap * height / (gap + height)
    re = props['rho'] * velocity_m_s * dh / props['mu']
    if re < 3000:
        return {'Re': re, 'Dh_m': dh, 'Nu': None, 'h_W_m2K': None,
                'dP_Pa': None, 'note': 'Gnielinski not valid; use laminar/developing rectangular-duct correlation'}
    f = (0.79 * np.log(re) - 1.64) ** -2
    pr = props['Pr']
    nu = (f / 8) * (re - 1000) * pr / (1 + 12.7 * np.sqrt(f / 8) * (pr ** (2/3) - 1))
    h = nu * props['k'] / dh
    k_ent = geometry.get('K_entrance', 0)
    k_exit = geometry.get('K_exit', 0)
    dp = (f * L / dh + k_ent + k_exit) * 0.5 * props['rho'] * velocity_m_s ** 2
    return {'Re': re, 'Dh_m': dh, 'Nu': nu, 'h_W_m2K': h,
            'dP_Pa': dp, 'note': 'Preliminary smooth-channel approximation'}


def main():
    p = isa_pressure_pa(ALTITUDE_M)
    props = air_properties(T_IN_C, p)
    heat_w = CPU_COUNT * CPU_W + GPU_COUNT * GPU_W + EXPANSION_W + OTHER_W
    mdot, vdot, cfm_req = required_airflow(heat_w, AIR_RISE_K, props)
    cfm_available = AIR_COOLED_CFM[-1]
    capacity = (cfm_available / CFM_PER_M3_S) * props['rho'] * props['cp'] * AIR_RISE_K
    failure_cfm = cfm_available * (N_FANS - 1) / N_FANS
    dp = serial_pressure_drop_pa(cfm_available, props['rho'])

    print('--- PRELIMINARY SERVER FAN ANALYSIS ---')
    print(f'Property model: {props["source"]}')
    print(f'Altitude: {ALTITUDE_M:.0f} m, inlet T: {T_IN_C:.1f} C, outlet target: {T_IN_C+AIR_RISE_K:.1f} C')
    print(f'Pressure: {p/1000:.2f} kPa, inlet density: {props["rho"]:.4f} kg/m3')
    print(f'cp: {props["cp"]:.1f} J/kg-K, mu: {props["mu"]:.3e} Pa-s, k: {props["k"]:.4f} W/m-K, Pr: {props["Pr"]:.3f}')
    print(f'Assumed heat load: {heat_w:.1f} W')
    print(f'Required mass flow: {mdot:.4f} kg/s')
    print(f'Required volumetric flow: {vdot:.4f} m3/s = {cfm_req:.1f} CFM')
    print(f'Graph 2 air-cooled max: ~{cfm_available:.0f} CFM (graph conditions unverified)')
    print(f'Heat removal at that flow and 15 K rise: {capacity:.0f} W')
    print(f'One-fan-out rough estimate: {failure_cfm:.0f} CFM (NOT validated)')
    print(f'Optional serial pressure drop: {dp:.1f} Pa' if dp is not None
          else 'Pressure drop NOT calculable without loss/geometry inputs.')
    if HEATSINK_CHANNEL is not None:
        result = heatsink_channel_analysis(float(LOCAL_VELOCITY_M_S['CPU-A'][-1]), props, HEATSINK_CHANNEL)
        print('Illustrative CPU-A channel estimate:', result)
    print('NOTE: graph 1 and graph 2 are distinct datasets; do not equate their CFM without verification.')

    powers = np.linspace(0, 6500, 200)
    req = np.array([required_airflow(power, AIR_RISE_K, props)[2] for power in powers])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(powers / 1000, req, label='Required inlet airflow')
    ax.axhline(cfm_available, linestyle='--', label='Graph 2 maximum (~240 CFM)')
    ax.scatter([heat_w / 1000], [cfm_req], label='Selected heat load', zorder=5)
    ax.set(xlabel='Heat transferred to air (kW)', ylabel='Inlet airflow (CFM)',
           title='Required airflow vs. reported flow at 3,000 m')
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig('required_vs_available_airflow.png', dpi=180)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(RPM, AIR_COOLED_CFM, 'o-', label='Air-cooled configuration')
    ax.plot(RPM, WATER_COOLED_CFM, 'o-', label='Water-cooled configuration')
    ax.set(xlabel='Fan speed (rpm)', ylabel='Reported airflow (CFM)',
           title='Graph 2: approximate visual readings')
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig('airflow_vs_fan_speed.png', dpi=180)

    fig, ax = plt.subplots(figsize=(10, 6))
    for name, values in LOCAL_VELOCITY_M_S.items():
        ax.plot(PWM, values, marker='o', linewidth=1.3, label=name)
    ax.set(xlabel='Main fan PWM (%)', ylabel='Local velocity (m/s)',
           title='Graph 1: all local velocity curves (approximate)')
    ax.grid(alpha=0.3)
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig('local_air_velocity_vs_pwm.png', dpi=180)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(PWM, PWM_TOTAL_AIRFLOW_CFM, 'o-', label='Graph 1: right-hand CFM axis')
    ax.set(xlabel='Main fan PWM (%)', ylabel='Reported total airflow (CFM)',
           title='Graph 1: total airflow vs. PWM (rough visual readings)')
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig('total_airflow_vs_pwm.png', dpi=180)
    plt.show()


if __name__ == '__main__':
    main()
