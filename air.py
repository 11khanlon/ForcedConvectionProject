#Air properties and thermal energy balance 
R_AIR = 287.05            #specific gas constant for dry
CFM = 2118.880003      #It comes from converting cubic meters per second into cubic feet per minute
#m^3/s. 2118.88 converts your volumetric flow rate (vdot) into CFM, which is the unit commonly used to describe fan airflow.

try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

#calculate atmospheric pressure at a given altitude
def isa_pressure_pa(altitude_m):
    if not 0 <= altitude_m <= 11000:
        raise ValueError('ISA troposphere model supports 0 to 11 km')
    return 101325.0 * (1 - 0.0065 * altitude_m / 288.15) ** 5.25588

#Use coolprop to find air properties 
def air_properties(temp_c, pressure_pa, use_coolprop=True):
    t = temp_c + 273.15
    if use_coolprop and PropsSI is not None:
        rho = PropsSI('Dmass', 'T', t, 'P', pressure_pa, 'Air')
        cp = PropsSI('Cpmass', 'T', t, 'P', pressure_pa, 'Air')
        mu = PropsSI('VISCOSITY', 'T', t, 'P', pressure_pa, 'Air')
        k = PropsSI('CONDUCTIVITY', 'T', t, 'P', pressure_pa, 'Air')
        source = 'CoolProp'
    else:
        rho, cp, mu, k = pressure_pa / (R_AIR * t), 1005., 1.9e-5, 0.0275
        source = 'ideal gas / approximate transport properties'
    return dict(rho=rho, cp=cp, mu=mu, k=k, Pr=cp * mu / k, source=source)

#calculate total heat load 
def heat_load_w(case):
    return (case['cpu_count'] * case['cpu_w'] + case['gpu_count'] * case['gpu_w']
            + case['expansion_w'] + case['other_w'])

#using Qdot = mdotcp(T_out - Iin), you can get mdot 
#Vdot = volumetric flow rate/ m^3/s
def required_airflow(heat_w, rise_k, props):
    if rise_k <= 0:
        raise ValueError('Air temperature rise must be positive')
    mdot = heat_w / (props['cp'] * rise_k)
    vdot = mdot / props['rho']
    return mdot, vdot, vdot * CFM

#Output: Heat carried by air
def heat_capacity_w(cfm, rise_k, props):
    return cfm / CFM * props['rho'] * props['cp'] * rise_k
