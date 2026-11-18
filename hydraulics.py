#Quadratic loss network, rectangular-duct friction, and fan intersection
#Determine actual air branch airflow --> vdot i and pressure loss Pi
import numpy as np
from air import CFM


def _positive(value, name):
    if value <= 0:
        raise ValueError(f'{name} must be > 0')
    return value

'''Darcy–Weisbach was originally developed using the diameter of a circular pipe. 
To account for this --> Calculate hydraulic diamter for rectangular duct 
'''
def hydraulic_diameter_rect(width_m, height_m):
    w, h = _positive(width_m, 'width'), _positive(height_m, 'height')
    return 2 * w * h / (w + h)


"""
This function calculates the Darcy friction factor, fd
Laminar 64/Re for circular ducts; turbulent Haaland approximation
Laminar rectangular channels require aspect-ratio correction elsewhere
Transitional Re 2300-4000 is not predicted reliably
"""
def friction_factor_darcy(re, relative_roughness=0.0):
    _positive(re, 'Re')
    if re < 2300:
        return 64 / re       #for laminar flow in circular pipe
    if re < 4000:
        return None
    return (-1.8 * np.log10((relative_roughness / 3.7) ** 1.11 + 6.9 / re)) ** -2   #Use Haaland approx for turbulent flow 


'''
Add correction factor for laminar region in pipe 
For rectangular passages, friction factor depends on aspect ratio
Fully-developed rectangular duct Poiseuille number f_D*Re
Shah-London polynomial; aspect ratio=min(w,h)/max(w,h)
A very narrow 
'''
def rectangular_laminar_friction_factor(re, width_m, height_m):
    a = min(width_m, height_m) / max(width_m, height_m)
    po = 96 * (1 - 1.3553*a + 1.9467*a**2 - 1.7012*a**3 + 0.9564*a**4 - 0.2537*a**5)
    return po / re


def branch_geometry(branch):
    typ = branch['type']
    if typ == 'equivalent':
        return _positive(branch['area_m2'], 'area'), None, 0.0, branch['K']
    if typ == 'rectangular':
        w, h = branch['width_m'], branch['height_m']
        return w*h, hydraulic_diameter_rect(w, h), branch['length_m'], branch.get('K_minor', 0.)
    if typ == 'fin_array':
        s, h, n = branch['gap_m'], branch['height_m'], branch['n_channels']
        if not isinstance(n, int) or n < 1:
            raise ValueError('n_channels must be positive integer')
        return n*s*h, hydraulic_diameter_rect(s, h), branch['length_m'], branch.get('K_minor', 0.)
    raise ValueError(f'Unsupported branch type: {typ}')


#Calculate pressure drop through one path
def branch_pressure_drop_pa(cfm, branch, props):
    area, dh, length, k_minor = branch_geometry(branch)
    _positive(area, 'flow area')
    if cfm < 0 or k_minor < 0 or length < 0:
        raise ValueError('Flow, K and length must be nonnegative')
    if cfm == 0:
        return 0.0
    velocity = (cfm / CFM) / area      #convert CFM to a usable velocity
    k_total = k_minor
    if dh is not None and length:
        re = props['rho'] * velocity * dh / props['mu']

        #Select either laminar or turbulent region 
        if re < 2300:
            if branch['type'] == 'rectangular':
                f = rectangular_laminar_friction_factor(re, branch['width_m'], branch['height_m'])
            else:
                f = rectangular_laminar_friction_factor(re, branch['gap_m'], branch['height_m'])
        else:
            f = friction_factor_darcy(re)
            if f is None:
                raise ValueError(f'Re={re:.0f} transitional: friction correlation not defined')
        k_total += f * length / dh
    return 0.5 * props['rho'] * velocity**2 * k_total


#add pressure losses in series 
def serial_pressure_drop_pa(cfm, losses, props):
    return sum(branch_pressure_drop_pa(cfm, {'type':'equivalent', 'area_m2':item['area_m2'], 'K':item['K']}, props)
               for item in losses.values())

#know pressure difference and need to find the airflow 
def _bisect(fn, lo, hi, iterations=75):
    flo, fhi = fn(lo), fn(hi)
    if flo == 0: return lo
    if fhi == 0: return hi
    if flo * fhi > 0:
        raise ValueError('Root not bracketed; check fan and resistance curves')
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        fm = fn(mid)
        if fm * flo > 0:
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)

#calculate flow splitting
def parallel_network_at_pressure(dp_pa, branches, props):
    if dp_pa < 0:
        raise ValueError('Pressure must be nonnegative')
    flows = {}
    for name, branch in branches.items():
        if dp_pa == 0:
            flows[name] = 0.0
            continue
        hi = 1.0
        while branch_pressure_drop_pa(hi, branch, props) < dp_pa:
            hi *= 2
            if hi > 1e7:
                raise ValueError(f'Could not bracket branch flow: {name}')
        flows[name] = _bisect(lambda q: branch_pressure_drop_pa(q, branch, props)-dp_pa, 0, hi)
    return flows


#Serial losses + ONE parallel group. If no branches, only serial losses
def network_pressure_drop_pa(total_cfm, serial_losses, branches, props):
    dp_serial = serial_pressure_drop_pa(total_cfm, serial_losses, props)
    if not branches:
        return dp_serial, {}
    if total_cfm == 0:
        return 0., {key:0. for key in branches}
    hi = 1.0
    while sum(parallel_network_at_pressure(hi, branches, props).values()) < total_cfm:
        hi *= 2
        if hi > 1e9:
            raise ValueError('Could not bracket parallel network pressure')
    dp_parallel = _bisect(lambda dp: sum(parallel_network_at_pressure(dp, branches, props).values())-total_cfm, 0, hi)
    return dp_serial + dp_parallel, parallel_network_at_pressure(dp_parallel, branches, props)

#Calculate pressure generated by fans 
#Density-corrected one-module fan curve, with identical modules in parallel
def fan_bank_pressure_pa(total_cfm, curve, props, n_active):
    q = np.asarray(curve['cfm'], float)
    dp = np.asarray(curve['pressure_pa'], float)
    if len(q) < 2 or len(q) != len(dp) or np.any(np.diff(q) <= 0):
        raise ValueError('Fan curve needs >=2 increasing CFM points')
    per_fan = total_cfm / n_active
    if per_fan < q[0] or per_fan > q[-1]:
        raise ValueError('Flow outside supplied fan curve; extrapolation disabled')
    return float(np.interp(per_fan, q, dp) * props['rho'] / curve['rho_ref'])


#Solve fan-bank pressure = system pressure without extrapolating fan curve
def operating_point(curve, serial_losses, branches, props, n_active):
    if not (serial_losses or branches):
        raise ValueError('System resistance is undefined: add geometry/loss coefficients')
    lo = n_active * curve['cfm'][0]
    hi = n_active * curve['cfm'][-1]
    def residual(q):
        return fan_bank_pressure_pa(q, curve, props, n_active) - network_pressure_drop_pa(q, serial_losses, branches, props)[0]
    q_op = _bisect(residual, lo, hi)
    dp, branch_flows = network_pressure_drop_pa(q_op, serial_losses, branches, props)
    return {'cfm':q_op, 'pressure_pa':dp, 'branch_cfm':branch_flows}
