from pathlib import Path
import config as cfg
from air import isa_pressure_pa, air_properties, heat_load_w, required_airflow, heat_capacity_w
from hydraulics import serial_pressure_drop_pa, network_pressure_drop_pa, operating_point
from convection import analyze_fin_channel
from plots import make_plots


def main():
    pressure = isa_pressure_pa(cfg.ALTITUDE_M)
    props = air_properties(cfg.T_IN_C, pressure, cfg.USE_COOLPROP)
    heat = heat_load_w(cfg.POWER_CASES[cfg.ACTIVE_POWER_CASE])
    mdot, vdot, required_cfm = required_airflow(heat, cfg.AIR_RISE_K, props)
    reported_cfm = cfg.AIR_COOLED_CFM[-1]

    print('--- PRELIMINARY SERVER FAN / THERMAL-FLUID ANALYSIS ---')
    print(f'Case: {cfg.ACTIVE_POWER_CASE}; heat: {heat:.0f} W')
    print(f'Air properties: {props["source"]}')
    print(f'Altitude: {cfg.ALTITUDE_M:.0f} m; inlet: {cfg.T_IN_C:.1f} °C; rise: {cfg.AIR_RISE_K:.1f} K')
    print(f'Atmospheric pressure: {pressure/1000:.2f} kPa; density: {props["rho"]:.4f} kg/m³')
    print(f'cp: {props["cp"]:.1f} J/(kg K); viscosity: {props["mu"]:.3e} Pa s')
    print(f'Conductivity: {props["k"]:.4f} W/(m K); Pr: {props["Pr"]:.3f}')
    print(f'Required mass flow: {mdot:.4f} kg/s; volumetric: {vdot:.4f} m³/s = {required_cfm:.1f} CFM')
    print(f'Paper air-cooled graph max: ~{reported_cfm:.0f} CFM (different conditions possible)')
    print(f'At this inlet density and 15 K rise, that flow transports: {heat_capacity_w(reported_cfm,cfg.AIR_RISE_K,props):.0f} W')
    print('Caution: paper airflow does not establish available airflow at the assignment operating conditions.')

    if cfg.SERIAL_LOSSES or cfg.PARALLEL_BRANCHES:
        dp, flows = network_pressure_drop_pa(reported_cfm, cfg.SERIAL_LOSSES, cfg.PARALLEL_BRANCHES, props)
        print(f'System pressure loss at {reported_cfm:.0f} CFM (ASSUMED network): {dp:.2f} Pa')
        for name, flow in flows.items():
            print(f'  Branch {name}: {flow:.2f} CFM')
    else:
        print('Pressure drop NOT calculated: specify actual/assumed loss geometry in config.py.')

    if cfg.FAN_CURVE is not None and (cfg.SERIAL_LOSSES or cfg.PARALLEL_BRANCHES):
        for active in ([cfg.N_FANS, cfg.N_FANS-1] if cfg.EVALUATE_ONE_FAN_FAILURE else [cfg.N_FANS]):
            try:
                op = operating_point(cfg.FAN_CURVE, cfg.SERIAL_LOSSES, cfg.PARALLEL_BRANCHES, props, active)
                print(f'{active} active fans: {op["cfm"]:.1f} CFM at {op["pressure_pa"]:.1f} Pa')
            except ValueError as exc:
                print(f'{active} active fans: operating point unavailable ({exc})')
    else:
        print('Fan operating point and fan-failure CFM NOT calculated: need measured fan curve + network.')

    if cfg.HEATSINK_CHANNEL is not None:
        print('Optional heatsink channel:', analyze_fin_channel(cfg.HEATSINK_CHANNEL, props))
    else:
        print('Heatsink Nu/h NOT calculated: channel geometry and channel-mean velocity/flow needed.')

    out = Path(__file__).resolve().parent / 'outputs'
    make_plots(cfg, props, heat, required_cfm, out)
    print(f'Plots saved to: {out}')


if __name__ == '__main__':
    main()
