#calculate heat-transfer performance of each channel: Re, Pr, Nu, hc 

import numpy as np
from air import CFM
from hydraulics import hydraulic_diameter_rect, rectangular_laminar_friction_factor


"""Straight rectangular channels; no heat transfer through the base included
    Uses fully-developed constant-wall-temperature rectangular laminar Nu, or Gnielinski for turbulent flow. 
    Transitional and developing heat transfer are not estimated. K_minor is referenced to channel velocity.
    """
def analyze_fin_channel(geometry, props):
    s, h, length = geometry['gap_m'], geometry['height_m'], geometry['length_m']
    dh = hydraulic_diameter_rect(s, h)
    n = geometry.get('n_channels', 1)
    if n <= 0 or length <= 0:
        raise ValueError('n_channels and channel length must be positive')
    area = n * s * h
    if 'velocity_m_s' in geometry:
        v = geometry['velocity_m_s']
    elif 'flow_cfm' in geometry:
        v = geometry['flow_cfm'] / CFM / area
    else:
        raise ValueError('Provide channel mean velocity_m_s or total flow_cfm')
    if v <= 0:
        raise ValueError('Channel velocity must be positive')
    re = props['rho'] * v * dh / props['mu']
    pr = props['Pr']
    g = length / (re * pr * dh)  # inverse thermal entrance length parameter
    note = ''
    if re < 2300:
        a = min(s,h)/max(s,h)
        nu = 7.541 * (1 - 2.610*a + 4.970*a*a - 5.119*a**3 + 2.702*a**4 - 0.548*a**5)
        f = rectangular_laminar_friction_factor(re, s, h)
        if g < 0.05:
            note = 'Thermal development may matter: fully developed Nu may be inaccurate.'
        else:
            note = 'Laminar, fully developed, constant wall temperature approximation.'
    elif re < 4000:
        return {'Re':re, 'Dh_m':dh, 'velocity_m_s':v, 'Nu':None, 'h_W_m2K':None,
                'dP_Pa':None, 'note':'Transitional flow: choose validated correlation.'}
    else:
        f = (0.79 * np.log(re) - 1.64) ** -2
        nu = (f/8) * (re-1000) * pr / (1 + 12.7*np.sqrt(f/8)*(pr**(2/3)-1))
        note = 'Gnielinski smooth-duct approximation; validate geometry/development.'
    htc = nu * props['k'] / dh
    dp = (f * length / dh + geometry.get('K_minor', 0)) * props['rho'] * v*v / 2
    return {'Re':re, 'Dh_m':dh, 'velocity_m_s':v, 'Nu':nu,
            'h_W_m2K':htc, 'dP_Pa':dp, 'note':note}
