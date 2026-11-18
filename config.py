from dataclasses import dataclass

T_IN_C = 35.0
ALTITUDE_M = 3000.0
AIR_RISE_K = 15.0
N_FANS = 4
USE_COOLPROP = True

# Assignment scenario is distinct from the actual IBM POWER AC922 paper
POWER_CASES = {
    'assignment': {'cpu_count': 2, 'cpu_w': 600.0, 'gpu_count': 4,
                   'gpu_w': 1000.0, 'expansion_w': 75.0, 'other_w': 0.0},
    'ibm_paper_cpu_gpu_only': {'cpu_count': 2, 'cpu_w': 190.0, 'gpu_count': 4,
                               'gpu_w': 300.0, 'expansion_w': 0.0, 'other_w': 0.0}
}

ACTIVE_POWER_CASE = 'assignment'

# Approximate visual readings from two distinct graphs; NOT digitized/validated
RPM = [3400., 6800., 10400.]
AIR_COOLED_CFM = [85., 160., 240.]
WATER_COOLED_CFM = [60., 120., 185.]
PWM = [30., 45., 60., 75., 100.]
LOCAL_VELOCITY_M_S = {
    'CPU-A': [1.8, 3.1, 4.4, 5.8, 6.6],
    'CPU-B': [2.7, 4.1, 5.5, 6.6, 7.4],
    'GPU-A1': [2.4, 3.8, 5.0, 6.2, 6.8],
    'GPU-A2': [1.9, 2.8, 3.8, 4.8, 5.5],
    'GPU-B1': [0.9, 1.8, 2.7, 3.5, 4.2],
    'GPU-B2': [1.9, 2.7, 3.7, 4.7, 5.5],
    'PCI-1': [2.0, 3.2, 4.6, 5.2, 6.1],
    'PCI-2': [0.5, 0.6, 0.8, 1.1, 1.4],
    'PCI-4-Exit': [2.5, 4.0, 5.4, 6.6, 7.5],
    'BXE': [3.2, 4.6, 5.9, 7.0, 7.8],
    'PS': [0.9, 1.2, 1.4, 1.8, 1.6]
}

PWM_TOTAL_AIRFLOW_CFM = [45., 135., 180., 230., 255.]

# All dimensions below are nominal or illustrative; do not claim IBM measurements.
#Double check CAD Data 
GEOMETRY_M = {
    'fan_frame_width': 0.080, 'fan_frame_height': 0.080,
    'fan_module_depth_assumed': 0.080,
    'pcie_card_length_max': 0.16765, 'pcie_card_height_max': 0.06890,
    'pcie_pcb_thickness_nominal': 0.00157,
    'pcie_slot_pitch_assumed': 0.02032}




# Optional serial elements, each with its OWN reference flow area and loss K.
# Example: 'inlet': {'area_m2': 0.020, 'K': 0.5}
SERIAL_LOSSES = {}

'''Parallel branches must connect the same upstream/downstream pressure nodes.
 Supported type 'equivalent': {'type':'equivalent','area_m2':...,'K':...}
Supported type 'rectangular': {'type':'rectangular','width_m':...,
   'height_m':...,'length_m':...,'K_minor':...}
For fin arrays use 'fin_array': {'type':'fin_array','gap_m':...,
 'height_m':...,'length_m':...,'n_channels':...,'K_minor':...}
K_minor refers to channel mean velocity for rectangular/fin_array types'''

PARALLEL_BRANCHES = {}
''''
Example: 
"pcie": {
        "type": "fin_array",
        "gap_m": 0.01875,
        "height_m": 0.0689,
        "n_channels": 3,
        "length_m": 0.16765,
        "K_minor": 1.5
    }

'''

'''Optional heat-transfer channel: requires known channel MEAN velocity or flow.
{'gap_m':0.002,'height_m':0.025,'length_m':0.10,
  'n_channels':20,'K_minor':1.5,'velocity_m_s':5.0}'''
HEATSINK_CHANNEL = None

# Manufacturer FAN CURVE: one fan module, pressure Pa vs volumetric CFM.
# IMPORTANT: Do not use max free-air CFM and shutoff pressure as a fabricated curve.
# Fill with digitized points at a SINGLE RPM/PWM and specified air density.
# Example structure only: {'cfm':[...], 'pressure_pa':[...], 'rho_ref':1.2}
FAN_CURVE = None

# If enabled and FAN_CURVE provided, also evaluate 3 identical active modules.
EVALUATE_ONE_FAN_FAILURE = True
