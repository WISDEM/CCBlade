from ccblade import CCAirfoil, CCBlade as CCBlade
from openmdao.api import Component
import numpy as np

cosd = lambda x: np.cos(np.deg2rad(x))
sind = lambda x: np.sin(np.deg2rad(x))

# ---------------------
# Base Components
# ---------------------


class GeometrySetupBase(Component):
    """a base component that computes the rotor radius from the geometry"""
    def __init__(self):
        super(GeometrySetupBase, self).__init__()
        self.add_output('R', val=0.0, units='m', desc='rotor radius')

# ---------------------
# Default Implementations of Base Classes
# ---------------------
class CCBladeGeometry(GeometrySetupBase):
    def __init__(self):
        super(CCBladeGeometry, self).__init__()
        self.add_param('Rtip', val=0.0, units='m', desc='tip radius')
        self.add_param('precurveTip', val=0.0, units='m', desc='tip radius')
        self.add_param('precone', val=0.0, desc='precone angle', units='deg')
        self.add_output('diameter', val=0.0, units='m')
        
    def solve_nonlinear(self, params, unknowns, resids):

        self.Rtip = params['Rtip']
        self.precurveTip = params['precurveTip']
        self.precone = params['precone']

        self.R = self.Rtip*cosd(self.precone) + self.precurveTip*sind(self.precone)

        unknowns['R'] = self.R
        unknowns['diameter'] = self.R*2

    def list_deriv_vars(self):

        inputs = ('Rtip', 'precurveTip', 'precone')
        outputs = ('R',)

        return inputs, outputs

    def linearize(self, params, unknowns, resids):

        J_sub = np.array([[cosd(self.precone), sind(self.precone),
            (-self.Rtip*sind(self.precone) + self.precurveTip*sind(self.precone))*np.pi/180.0]])
        J = {}
        J['R', 'Rtip'] = J_sub[0][0]
        J['R', 'precurveTip'] = J_sub[0][1]
        J['R', 'precone'] = J_sub[0][2]
        J['diameter', 'Rtip'] = 2.0*J_sub[0][0]
        J['diameter', 'precurveTip'] = 2.0*J_sub[0][1]
        J['diameter', 'precone'] = 2.0*J_sub[0][2]
        J['diameter', 'R'] = 2.0

        return J


class CCBladePower(Component):
    def __init__(self, naero, npower):
        super(CCBladePower, self).__init__()
        """blade element momentum code"""

        # inputs
        self.add_param('Uhub', val=np.zeros(npower), units='m/s', desc='hub height wind speed')
        self.add_param('Omega', val=np.zeros(npower), units='rpm', desc='rotor rotation speed')
        self.add_param('pitch', val=np.zeros(npower), units='deg', desc='blade pitch setting')

        # outputs
        self.add_output('T', val=np.zeros(npower), units='N', desc='rotor aerodynamic thrust')
        self.add_output('Q', val=np.zeros(npower), units='N*m', desc='rotor aerodynamic torque')
        self.add_output('P', val=np.zeros(npower), units='W', desc='rotor aerodynamic power')

        
        # (potential) variables
        self.add_param('r', val=np.zeros(naero), units='m', desc='radial locations where blade is defined (should be increasing and not go all the way to hub or tip)')
        self.add_param('chord', val=np.zeros(naero), units='m', desc='chord length at each section')
        self.add_param('theta', val=np.zeros(naero),  units='deg', desc='twist angle at each section (positive decreases angle of attack)')
        self.add_param('Rhub', val=0.0, units='m', desc='hub radius')
        self.add_param('Rtip', val=0.0, units='m', desc='tip radius')
        self.add_param('hubHt', val=0.0, units='m', desc='hub height')
        self.add_param('precone', val=0.0, desc='precone angle', units='deg')
        self.add_param('tilt', val=0.0, desc='shaft tilt', units='deg')
        self.add_param('yaw', val=0.0, desc='yaw error', units='deg')

        # TODO: I've not hooked up the gradients for these ones yet.
        self.add_param('precurve', val=np.zeros(naero), units='m', desc='precurve at each section')
        self.add_param('precurveTip', val=0.0, units='m', desc='precurve at tip')

        # parameters
        self.add_param('airfoil_files', val=[0]*naero, desc='names of airfoil file', pass_by_obj=True)
        self.add_param('B', val=0, desc='number of blades', pass_by_obj=True)
        self.add_param('rho', val=0.0, units='kg/m**3', desc='density of air')
        self.add_param('mu', val=0.0, units='kg/(m*s)', desc='dynamic viscosity of air')
        self.add_param('shearExp', val=0.0, desc='shear exponent')
        self.add_param('nSector', val=4, desc='number of sectors to divide rotor face into in computing thrust and power', pass_by_obj=True)
        self.add_param('tiploss', val=True, desc='include Prandtl tip loss model', pass_by_obj=True)
        self.add_param('hubloss', val=True, desc='include Prandtl hub loss model', pass_by_obj=True)
        self.add_param('wakerotation', val=True, desc='include effect of wake rotation (i.e., tangential induction factor is nonzero)', pass_by_obj=True)
        self.add_param('usecd', val=True, desc='use drag coefficient in computing induction factors', pass_by_obj=True)

	self.naero = naero
	self.deriv_options['form'] = 'central'
	self.deriv_options['step_calc'] = 'relative'
        
    def solve_nonlinear(self, params, unknowns, resids):

        self.r = params['r']
        self.chord = params['chord']
        self.theta = params['theta']
        self.Rhub = params['Rhub']
        self.Rtip = params['Rtip']
        self.hubHt = params['hubHt']
        self.precone = params['precone']
        self.tilt = params['tilt']
        self.yaw = params['yaw']
        self.precurve = params['precurve']
        self.precurveTip = params['precurveTip']
        self.airfoil_files = params['airfoil_files']
        self.B = params['B']
        self.rho = params['rho']
        self.mu = params['mu']
        self.shearExp = params['shearExp']
        self.nSector = params['nSector']
        self.tiploss = params['tiploss']
        self.hubloss = params['hubloss']
        self.wakerotation = params['wakerotation']
        self.usecd = params['usecd']
        self.Uhub = params['Uhub']
        self.Omega = params['Omega']
        self.pitch = params['pitch']
        
        if len(self.precurve) == 0:
            self.precurve = np.zeros_like(self.r)

        # airfoil files
        n = len(self.airfoil_files)
        af = [0]*n
        for i in range(n):
            af[i] = CCAirfoil.initFromAerodynFile(self.airfoil_files[i])

        self.ccblade = CCBlade(self.r, self.chord, self.theta, af, self.Rhub, self.Rtip, self.B,
            self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp, self.hubHt,
            self.nSector, self.precurve, self.precurveTip, tiploss=self.tiploss, hubloss=self.hubloss,
            wakerotation=self.wakerotation, usecd=self.usecd, derivatives=True)

        # power, thrust, torque
        self.P, self.T, self.Q, self.dP, self.dT, self.dQ \
            = self.ccblade.evaluate(self.Uhub, self.Omega, self.pitch, coefficient=False)
        unknowns['T'] = self.T
        unknowns['Q'] = self.Q
        unknowns['P'] = self.P
        

    def list_deriv_vars(self):

        inputs = ('precone', 'tilt', 'hubHt', 'Rhub', 'Rtip', 'yaw',
                  'Uhub', 'Omega', 'pitch', 'r', 'chord', 'theta', 'precurve', 'precurveTip')
        outputs = ('P', 'T', 'Q')
        return inputs, outputs


    def linearize(self, params, unknowns, resids):

        dP = self.dP
        dT = self.dT
        dQ = self.dQ

        J = {}
        J['P', 'precone'] = dP['dprecone']
        J['P', 'tilt'] = dP['dtilt']
        J['P', 'hubHt'] = dP['dhubHt']
        J['P', 'Rhub'] = dP['dRhub']
        J['P', 'Rtip'] = dP['dRtip']
        J['P', 'yaw'] = dP['dyaw']
        J['P', 'Uhub'] = dP['dUinf']
        J['P', 'Omega'] = dP['dOmega']
        J['P', 'pitch'] =  dP['dpitch']
        J['P', 'r'] = dP['dr']
        J['P', 'chord'] = dP['dchord']
        J['P', 'theta'] = dP['dtheta']
        J['P', 'precurve'] = dP['dprecurve']
        J['P', 'precurveTip'] = dP['dprecurveTip']

        J['T', 'precone'] = dT['dprecone']
        J['T', 'tilt'] = dT['dtilt']
        J['T', 'hubHt'] = dT['dhubHt']
        J['T', 'Rhub'] = dT['dRhub']
        J['T', 'Rtip'] = dT['dRtip']
        J['T', 'yaw'] = dT['dyaw']
        J['T', 'Uhub'] = dT['dUinf']
        J['T', 'Omega'] = dT['dOmega']
        J['T', 'pitch'] =  dT['dpitch']
        J['T', 'r'] = dT['dr']
        J['T', 'chord'] = dT['dchord']
        J['T', 'theta'] = dT['dtheta']
        J['T', 'precurve'] = dT['dprecurve']
        J['T', 'precurveTip'] = dT['dprecurveTip']

        J['Q', 'precone'] = dQ['dprecone']
        J['Q', 'tilt'] = dQ['dtilt']
        J['Q', 'hubHt'] = dQ['dhubHt']
        J['Q', 'Rhub'] = dQ['dRhub']
        J['Q', 'Rtip'] = dQ['dRtip']
        J['Q', 'yaw'] = dQ['dyaw']
        J['Q', 'Uhub'] = dQ['dUinf']
        J['Q', 'Omega'] = dQ['dOmega']
        J['Q', 'pitch'] =  dQ['dpitch']
        J['Q', 'r'] = dQ['dr']
        J['Q', 'chord'] = dQ['dchord']
        J['Q', 'theta'] = dQ['dtheta']
        J['Q', 'precurve'] = dQ['dprecurve']
        J['Q', 'precurveTip'] = dQ['dprecurveTip']

        return J


    
class CCBladeLoads(Component):
    def __init__(self, naero, npower):
        super(CCBladeLoads, self).__init__()
        """blade element momentum code"""

        # inputs
        self.add_param('V_load', val=0.0, units='m/s', desc='hub height wind speed')
        self.add_param('Omega_load', val=0.0, units='rpm', desc='rotor rotation speed')
        self.add_param('pitch_load', val=0.0, units='deg', desc='blade pitch setting')
        self.add_param('azimuth_load', val=0.0, units='deg', desc='blade azimuthal location')

        # outputs
        self.add_output('loads:r', val=np.zeros(naero+2), units='m', desc='radial positions along blade going toward tip')
        self.add_output('loads:Px', val=np.zeros(naero+2), units='N/m', desc='distributed loads in blade-aligned x-direction')
        self.add_output('loads:Py', val=np.zeros(naero+2), units='N/m', desc='distributed loads in blade-aligned y-direction')
        self.add_output('loads:Pz', val=np.zeros(naero+2), units='N/m', desc='distributed loads in blade-aligned z-direction')

        # corresponding setting for loads
        self.add_output('loads:V', val=0.0, units='m/s', desc='hub height wind speed')
        self.add_output('loads:Omega', val=0.0, units='rpm', desc='rotor rotation speed')
        self.add_output('loads:pitch', val=0.0, units='deg', desc='pitch angle')
        self.add_output('loads:azimuth', val=0.0, units='deg', desc='azimuthal angle')


        
        # (potential) variables
        self.add_param('r', val=np.zeros(naero), units='m', desc='radial locations where blade is defined (should be increasing and not go all the way to hub or tip)')
        self.add_param('chord', val=np.zeros(naero), units='m', desc='chord length at each section')
        self.add_param('theta', val=np.zeros(naero),  units='deg', desc='twist angle at each section (positive decreases angle of attack)')
        self.add_param('Rhub', val=0.0, units='m', desc='hub radius')
        self.add_param('Rtip', val=0.0, units='m', desc='tip radius')
        self.add_param('hubHt', val=0.0, units='m', desc='hub height')
        self.add_param('precone', val=0.0, desc='precone angle', units='deg')
        self.add_param('tilt', val=0.0, desc='shaft tilt', units='deg')
        self.add_param('yaw', val=0.0, desc='yaw error', units='deg')

        # TODO: I've not hooked up the gradients for these ones yet.
        self.add_param('precurve', val=np.zeros(naero), units='m', desc='precurve at each section')
        self.add_param('precurveTip', val=0.0, units='m', desc='precurve at tip')

        # parameters
        self.add_param('airfoil_files', val=[0]*naero, desc='names of airfoil file', pass_by_obj=True)
        self.add_param('B', val=0, desc='number of blades', pass_by_obj=True)
        self.add_param('rho', val=0.0, units='kg/m**3', desc='density of air')
        self.add_param('mu', val=0.0, units='kg/(m*s)', desc='dynamic viscosity of air')
        self.add_param('shearExp', val=0.0, desc='shear exponent')
        self.add_param('nSector', val=4, desc='number of sectors to divide rotor face into in computing thrust and power', pass_by_obj=True)
        self.add_param('tiploss', val=True, desc='include Prandtl tip loss model', pass_by_obj=True)
        self.add_param('hubloss', val=True, desc='include Prandtl hub loss model', pass_by_obj=True)
        self.add_param('wakerotation', val=True, desc='include effect of wake rotation (i.e., tangential induction factor is nonzero)', pass_by_obj=True)
        self.add_param('usecd', val=True, desc='use drag coefficient in computing induction factors', pass_by_obj=True)

	self.naero = naero
	self.deriv_options['form'] = 'central'
	self.deriv_options['step_calc'] = 'relative'
        
    def solve_nonlinear(self, params, unknowns, resids):

        self.r = params['r']
        self.chord = params['chord']
        self.theta = params['theta']
        self.Rhub = params['Rhub']
        self.Rtip = params['Rtip']
        self.hubHt = params['hubHt']
        self.precone = params['precone']
        self.tilt = params['tilt']
        self.yaw = params['yaw']
        self.precurve = params['precurve']
        self.precurveTip = params['precurveTip']
        self.airfoil_files = params['airfoil_files']
        self.B = params['B']
        self.rho = params['rho']
        self.mu = params['mu']
        self.shearExp = params['shearExp']
        self.nSector = params['nSector']
        self.tiploss = params['tiploss']
        self.hubloss = params['hubloss']
        self.wakerotation = params['wakerotation']
        self.usecd = params['usecd']
        self.V_load = params['V_load']
        self.Omega_load = params['Omega_load']
        self.pitch_load = params['pitch_load']
        self.azimuth_load = params['azimuth_load']


        if len(self.precurve) == 0:
            self.precurve = np.zeros_like(self.r)

        # airfoil files
        n = len(self.airfoil_files)
        af = [0]*n
        for i in range(n):
            af[i] = CCAirfoil.initFromAerodynFile(self.airfoil_files[i])

        self.ccblade = CCBlade(self.r, self.chord, self.theta, af, self.Rhub, self.Rtip, self.B,
            self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp, self.hubHt,
            self.nSector, self.precurve, self.precurveTip, tiploss=self.tiploss, hubloss=self.hubloss,
            wakerotation=self.wakerotation, usecd=self.usecd, derivatives=True)

        # distributed loads
        Np, Tp, self.dNp, self.dTp \
            = self.ccblade.distributedAeroLoads(self.V_load, self.Omega_load, self.pitch_load, self.azimuth_load)

        # concatenate loads at root/tip
        unknowns['loads:r'] = np.concatenate([[self.Rhub], self.r, [self.Rtip]])
        Np = np.concatenate([[0.0], Np, [0.0]])
        Tp = np.concatenate([[0.0], Tp, [0.0]])

        # conform to blade-aligned coordinate system
        unknowns['loads:Px'] = Np
        unknowns['loads:Py'] = -Tp
        unknowns['loads:Pz'] = 0*Np

        # return other outputs needed
        unknowns['loads:V'] = self.V_load
        unknowns['loads:Omega'] = self.Omega_load
        unknowns['loads:pitch'] = self.pitch_load
        unknowns['loads:azimuth'] = self.azimuth_load


    def list_deriv_vars(self):
        inputs = ('r', 'chord', 'theta', 'Rhub', 'Rtip', 'hubHt', 'precone',
                  'tilt', 'yaw', 'V_load', 'Omega_load', 'pitch_load', 'azimuth_load', 'precurve')
        outputs = ('loads:r', 'loads:Px', 'loads:Py', 'loads:Pz', 'loads:V',
                   'loads:Omega', 'loads:pitch', 'loads:azimuth')

        return inputs, outputs


    def linearize(self, params, unknowns, resids):

        dNp = self.dNp
        dTp = self.dTp
        n = len(self.r)

        dr_dr = np.vstack([np.zeros(n), np.eye(n), np.zeros(n)])
        dr_dRhub = np.zeros(n+2)
        dr_dRtip = np.zeros(n+2)
        dr_dRhub[0] = 1.0
        dr_dRtip[-1] = 1.0

        dV = np.zeros(4*n+10)
        dV[3*n+6] = 1.0
        dOmega = np.zeros(4*n+10)
        dOmega[3*n+7] = 1.0
        dpitch = np.zeros(4*n+10)
        dpitch[3*n+8] = 1.0
        dazimuth = np.zeros(4*n+10)
        dazimuth[3*n+9] = 1.0

        J = {}
        zero = np.zeros(self.naero)
        J['loads:r', 'r'] = dr_dr
        J['loads:r', 'Rhub'] = dr_dRhub
        J['loads:r', 'Rtip'] = dr_dRtip
        J['loads:Px', 'r'] = np.vstack([zero, dNp['dr'], zero])
        J['loads:Px', 'chord'] = np.vstack([zero, dNp['dchord'], zero])
        J['loads:Px', 'theta'] = np.vstack([zero, dNp['dtheta'], zero])
        J['loads:Px', 'Rhub'] = np.concatenate([[0.0], np.squeeze(dNp['dRhub']), [0.0]])
        J['loads:Px', 'Rtip'] = np.concatenate([[0.0], np.squeeze(dNp['dRtip']), [0.0]])
        J['loads:Px', 'hubHt'] = np.concatenate([[0.0], np.squeeze(dNp['dhubHt']), [0.0]])
        J['loads:Px', 'precone'] = np.concatenate([[0.0], np.squeeze(dNp['dprecone']), [0.0]])
        J['loads:Px', 'tilt'] = np.concatenate([[0.0], np.squeeze(dNp['dtilt']), [0.0]])
        J['loads:Px', 'yaw'] = np.concatenate([[0.0], np.squeeze(dNp['dyaw']), [0.0]])
        J['loads:Px', 'V_load'] = np.concatenate([[0.0], np.squeeze(dNp['dUinf']), [0.0]])
        J['loads:Px', 'Omega_load'] = np.concatenate([[0.0], np.squeeze(dNp['dOmega']), [0.0]])
        J['loads:Px', 'pitch_load'] = np.concatenate([[0.0], np.squeeze(dNp['dpitch']), [0.0]])
        J['loads:Px', 'azimuth_load'] = np.concatenate([[0.0], np.squeeze(dNp['dazimuth']), [0.0]])
        J['loads:Px', 'precurve'] = np.vstack([zero, dNp['dprecurve'], zero])
        J['loads:Py', 'r'] = np.vstack([zero, -dTp['dr'], zero])
        J['loads:Py', 'chord'] = np.vstack([zero, -dTp['dchord'], zero])
        J['loads:Py', 'theta'] = np.vstack([zero, -dTp['dtheta'], zero])
        J['loads:Py', 'Rhub'] = np.concatenate([[0.0], -np.squeeze(dTp['dRhub']), [0.0]])
        J['loads:Py', 'Rtip'] = np.concatenate([[0.0], -np.squeeze(dTp['dRtip']), [0.0]])
        J['loads:Py', 'hubHt'] = np.concatenate([[0.0], -np.squeeze(dTp['dhubHt']), [0.0]])
        J['loads:Py', 'precone'] = np.concatenate([[0.0], -np.squeeze(dTp['dprecone']), [0.0]])
        J['loads:Py', 'tilt'] = np.concatenate([[0.0], -np.squeeze(dTp['dtilt']), [0.0]])
        J['loads:Py', 'yaw'] = np.concatenate([[0.0], -np.squeeze(dTp['dyaw']), [0.0]])
        J['loads:Py', 'V_load'] = np.concatenate([[0.0], -np.squeeze(dTp['dUinf']), [0.0]])
        J['loads:Py', 'Omega_load'] = np.concatenate([[0.0], -np.squeeze(dTp['dOmega']), [0.0]])
        J['loads:Py', 'pitch_load'] = np.concatenate([[0.0], -np.squeeze(dTp['dpitch']), [0.0]])
        J['loads:Py', 'azimuth_load'] = np.concatenate([[0.0], -np.squeeze(dTp['dazimuth']), [0.0]])
        J['loads:Py', 'precurve'] = np.vstack([zero, -dTp['dprecurve'], zero])
        J['loads:V', 'V_load'] = 1.0
        J['loads:Omega', 'Omega_load'] = 1.0
        J['loads:pitch', 'pitch_load'] = 1.0
        J['loads:azimuth', 'azimuth_load'] = 1.0

        return J
    
'''
def common_io(group, varspeed, varpitch):

    regulated = varspeed or varpitch

    # add inputs
    group.add_param('npts_coarse_power_curve', val=20, desc='number of points to evaluate aero analysis at')
    group.add_param('npts_spline_power_curve', val=200, desc='number of points to use in fitting spline to power curve')
    group.add_param('AEP_loss_factor', val=1.0, desc='availability and other losses (soiling, array, etc.)')
    if varspeed:
        group.add_param('control:Vin', units='m/s', desc='cut-in wind speed')
        group.add_param('control:Vout', units='m/s', desc='cut-out wind speed')
        group.add_param('control:ratedPower', units='W', desc='rated power')
        group.add_param('control:minOmega', units='rpm', desc='minimum allowed rotor rotation speed')
        group.add_param('control:maxOmega', units='rpm', desc='maximum allowed rotor rotation speed')
        group.add_param('control:tsr', desc='tip-speed ratio in Region 2 (should be optimized externally)')
        group.add_param('control:pitch', units='deg', desc='pitch angle in region 2 (and region 3 for fixed pitch machines)')
    else:
        group.add_param('control:Vin', units='m/s', desc='cut-in wind speed')
        group.add_param('control:Vout', units='m/s', desc='cut-out wind speed')
        group.add_param('control:ratedPower', units='W', desc='rated power')
        group.add_param('control:Omega', units='rpm', desc='fixed rotor rotation speed')
        group.add_param('control:pitch', units='deg', desc='pitch angle in region 2 (and region 3 for fixed pitch machines)')
        group.add_param('control:npts', val=20, desc='number of points to evalute aero code to generate power curve')


    # # add slots (must replace)
    # group.add('geom', Slot(GeomtrySetupBase))
    # group.add('analysis', Slot(AeroBase))
    # group.add('dt', Slot(DrivetrainLossesBase))
    # group.add('cdf', Slot(CDFBase))


    # add outputs
    group.add_output('AEP', units='kW*h', desc='annual energy production')
    group.add_output('V', units='m/s', desc='wind speeds (power curve)')
    group.add_output('P', units='W', desc='power (power curve)')
    group.add_output('diameter', units='m', desc='rotor diameter')
    if regulated:
        group.add_output('ratedConditions:V', units='m/s', desc='rated wind speed')
        group.add_output('ratedConditions:Omega', units='rpm', desc='rotor rotation speed at rated')
        group.add_output('ratedConditions:pitch', units='deg', desc='pitch setting at rated')
        group.add_output('ratedConditions:T', units='N', desc='rotor aerodynamic thrust at rated')
        group.add_output('ratedConditions:Q', units='N*m', desc='rotor aerodynamic torque at rated')


def common_configure(group, varspeed, varpitch):

    regulated = varspeed or varpitch

    # add components
    group.add('geom', GeomtrySetupBase())

    if varspeed:
        group.add('setup', SetupRunVarSpeed(20))
    else:
        group.add('setup', SetupRunFixedSpeed())

    group.add('analysis', AeroBase())
    group.add('dt', DrivetrainLossesBase())

    if varspeed or varpitch:
        group.add('powercurve', RegulatedPowerCurve(20))
        group.add('brent', Brent())
        group.brent.workflow.add(['powercurve'])
    else:
        group.add('powercurve', UnregulatedPowerCurve())

    group.add('cdf', CDFBase())
    group.add('aep', AEP(200))

    if regulated:
        group.driver.workflow.add(['geom', 'setup', 'analysis', 'dt', 'brent', 'cdf', 'aep'])
    else:
        group.driver.workflow.add(['geom', 'setup', 'analysis', 'dt', 'powercurve', 'cdf', 'aep'])


    # connections to setup
    group.connect('control', 'setup.control')
    group.connect('npts_coarse_power_curve', 'setup.npts')
    if varspeed:
        group.connect('geom.R', 'setup.R')


    # connections to analysis
    group.connect('setup.Uhub', 'analysis.Uhub')
    group.connect('setup.Omega', 'analysis.Omega')
    group.connect('setup.pitch', 'analysis.pitch')
    group.analysis.run_case = 'power'


    # connections to drivetrain
    group.connect('analysis.P', 'dt.aeroPower')
    group.connect('analysis.Q', 'dt.aeroTorque')
    group.connect('analysis.T', 'dt.aeroThrust')
    group.connect('control:ratedPower', 'dt.ratedPower')


    # connections to powercurve
    group.connect('control', 'powercurve.control')
    group.connect('setup.Uhub', 'powercurve.Vcoarse')
    group.connect('dt.power', 'powercurve.Pcoarse')
    group.connect('analysis.T', 'powercurve.Tcoarse')
    group.connect('npts_spline_power_curve', 'powercurve.npts')

    if regulated:
        group.connect('geom.R', 'powercurve.R')

        # setup Brent method to find rated speed
        group.connect('control:Vin', 'brent.lower_bound')
        group.connect('control:Vout', 'brent.upper_bound')
        group.brent.add_parameter('powercurve.Vrated', low=-1e-15, high=1e15)
        group.brent.add_constraint('powercurve.residual = 0')
        group.brent.invalid_bracket_return = 1.0


    # connections to cdf
    group.connect('powercurve.V', 'cdf.x')


    # connections to aep
    group.connect('cdf.F', 'aep.CDF_V')
    group.connect('powercurve.P', 'aep.P')
    group.connect('AEP_loss_factor', 'aep.lossFactor')


    # connections to outputs
    group.connect('powercurve.V', 'V')
    group.connect('powercurve.P', 'P')
    group.connect('aep.AEP', 'AEP')
    group.connect('2*geom.R', 'diameter')
    if regulated:
        group.connect('powercurve.ratedConditions', 'ratedConditions')





def common_io_with_ccblade(group, varspeed, varpitch, cdf_type):

    regulated = varspeed or varpitch

    # add inputs
    group.add_param('r_af', units='m', desc='locations where airfoils are defined on unit radius')
    group.add_param('r_max_chord')
    group.add_param('chord_sub', units='m', desc='chord at control points')
    group.add_param('theta_sub', units='deg', desc='twist at control points')
    group.add_param('Rhub', units='m', desc='hub radius')
    group.add_param('Rtip', units='m', desc='tip radius')
    group.add_param('hubHt', units='m')
    group.add_param('precone', desc='precone angle', units='deg')
    group.add_param('tilt', val=0.0, desc='shaft tilt', units='deg')
    group.add_param('yaw', val=0.0, desc='yaw error', units='deg')
    group.add_param('airfoil_files', desc='names of airfoil file')
    group.add_param('idx_cylinder', desc='location where cylinder section ends on unit radius')
    group.add_param('B', val=3, desc='number of blades')
    group.add_param('rho', val=1.225, units='kg/m**3', desc='density of air')
    group.add_param('mu', val=1.81206e-5, units='kg/m/s', desc='dynamic viscosity of air')
    group.add_param('shearExp', val=0.2, desc='shear exponent')
    group.add_param('nSector', val=4, desc='number of sectors to divide rotor face into in computing thrust and power')
    group.add_param('tiploss', val=True, desc='include Prandtl tip loss model')
    group.add_param('hubloss', val=True, desc='include Prandtl hub loss model')
    group.add_param('wakerotation', val=True, desc='include effect of wake rotation (i.e., tangential induction factor is nonzero)')
    group.add_param('usecd', val=True, desc='use drag coefficient in computing induction factors')
    group.add_param('npts_coarse_power_curve', val=20, desc='number of points to evaluate aero analysis at')
    group.add_param('npts_spline_power_curve', val=200, desc='number of points to use in fitting spline to power curve')
    group.add_param('AEP_loss_factor', val=1.0, desc='availability and other losses (soiling, array, etc.)')

    if varspeed:
        group.add_param('control:Vin', units='m/s', desc='cut-in wind speed')
        group.add_param('control:Vout', units='m/s', desc='cut-out wind speed')
        group.add_param('control:ratedPower', units='W', desc='rated power')
        group.add_param('control:minOmega', units='rpm', desc='minimum allowed rotor rotation speed')
        group.add_param('control:maxOmega', units='rpm', desc='maximum allowed rotor rotation speed')
        group.add_param('control:tsr', desc='tip-speed ratio in Region 2 (should be optimized externally)')
        group.add_param('control:pitch', units='deg', desc='pitch angle in region 2 (and region 3 for fixed pitch machines)')
    else:
        group.add_param('control:Vin', units='m/s', desc='cut-in wind speed')
        group.add_param('control:Vout', units='m/s', desc='cut-out wind speed')
        group.add_param('control:ratedPower', units='W', desc='rated power')
        group.add_param('control:Omega', units='rpm', desc='fixed rotor rotation speed')
        group.add_param('control:pitch', units='deg', desc='pitch angle in region 2 (and region 3 for fixed pitch machines)')
        group.add_param('control:npts', val=20, desc='number of points to evalute aero code to generate power curve')

    group.add_param('drivetrainType', val=DRIVETRAIN_TYPE['GEARED'])
    group.add_param('cdf_mean_wind_speed', units='m/s', desc='mean wind speed of site cumulative distribution function')

    if cdf_type == 'weibull':
        group.add_param('weibull_shape_factor', desc='(shape factor of weibull distribution)')

    # outputs
    group.add_output('AEP', units='kW*h', desc='annual energy production')
    group.add_output('V', units='m/s', desc='wind speeds (power curve)')
    group.add_output('P', units='W', desc='power (power curve)')
    group.add_output('diameter', units='m')
    if regulated:
        group.add_output('ratedConditions:V', units='m/s', desc='rated wind speed')
        group.add_output('ratedConditions:Omega', units='rpm', desc='rotor rotation speed at rated')
        group.add_output('ratedConditions:pitch', units='deg', desc='pitch setting at rated')
        group.add_output('ratedConditions:T', units='N', desc='rotor aerodynamic thrust at rated')
        group.add_output('ratedConditions:Q', units='N*m', desc='rotor aerodynamic torque at rated')



def common_configure_with_ccblade(group, varspeed, varpitch, cdf_type):
    common_configure(group, varspeed, varpitch)

    # put in parameterization for CCBlade
    group.add('spline', GeometrySpline())
    group.replace('geom', CCBladeGeometry())
    group.replace('analysis', CCBlade())
    group.replace('dt', CSMDrivetrain())
    if cdf_type == 'rayleigh':
        group.replace('cdf', RayleighCDF())
    elif cdf_type == 'weibull':
        group.replace('cdf', WeibullWithMeanCDF())


    # add spline to workflow
    group.driver.workflow.add('spline')

    # connections to spline
    group.connect('r_af', 'spline.r_af')
    group.connect('r_max_chord', 'spline.r_max_chord')
    group.connect('chord_sub', 'spline.chord_sub')
    group.connect('theta_sub', 'spline.theta_sub')
    group.connect('idx_cylinder', 'spline.idx_cylinder')
    group.connect('Rhub', 'spline.Rhub')
    group.connect('Rtip', 'spline.Rtip')

    # connections to geom
    group.connect('Rtip', 'geom.Rtip')
    group.connect('precone', 'geom.precone')

    # connections to analysis
    group.connect('spline.r', 'analysis.r')
    group.connect('spline.chord', 'analysis.chord')
    group.connect('spline.theta', 'analysis.theta')
    group.connect('spline.precurve', 'analysis.precurve')
    group.connect('Rhub', 'analysis.Rhub')
    group.connect('Rtip', 'analysis.Rtip')
    group.connect('hubHt', 'analysis.hubHt')
    group.connect('precone', 'analysis.precone')
    group.connect('tilt', 'analysis.tilt')
    group.connect('yaw', 'analysis.yaw')
    group.connect('airfoil_files', 'analysis.airfoil_files')
    group.connect('B', 'analysis.B')
    group.connect('rho', 'analysis.rho')
    group.connect('mu', 'analysis.mu')
    group.connect('shearExp', 'analysis.shearExp')
    group.connect('nSector', 'analysis.nSector')
    group.connect('tiploss', 'analysis.tiploss')
    group.connect('hubloss', 'analysis.hubloss')
    group.connect('wakerotation', 'analysis.wakerotation')
    group.connect('usecd', 'analysis.usecd')

    # connections to dt
    group.connect('drivetrainType', 'dt.drivetrainType')
    group.dt.missing_deriv_policy = 'assume_zero'  # TODO: openmdao bug remove later

    # connnections to cdf
    group.connect('cdf_mean_wind_speed', 'cdf.xbar')
    if cdf_type == 'weibull':
        group.connect('weibull_shape_factor', 'cdf.k')




class RotorAeroVSVP(Group):

    def configure(self):
        varspeed = True
        varpitch = True
        common_io(self, varspeed, varpitch)
        common_configure(self, varspeed, varpitch)


class RotorAeroVSFP(Group):

    def configure(self):
        varspeed = True
        varpitch = False
        common_io(self, varspeed, varpitch)
        common_configure(self, varspeed, varpitch)


class RotorAeroFSVP(Group):

    def configure(self):
        varspeed = False
        varpitch = True
        common_io(self, varspeed, varpitch)
        common_configure(self, varspeed, varpitch)


class RotorAeroFSFP(Group):

    def configure(self):
        varspeed = False
        varpitch = False
        common_io(self, varspeed, varpitch)
        common_configure(self, varspeed, varpitch)

        
class RotorAeroVSVPWithCCBlade(Group):
    def __init__(self, cdf_type='weibull'):
        super(RotorAeroVSVPWithCCBlade, self).__init__()
        self.cdf_type = cdf_type

    def configure(self):
        varspeed = True
        varpitch = True
        common_io_with_ccblade(self, varspeed, varpitch, self.cdf_type)
        common_configure_with_ccblade(self, varspeed, varpitch, self.cdf_type)


class RotorAeroVSFPWithCCBlade(Group):

    def __init__(self, cdf_type='weibull'):
        self.cdf_type = cdf_type
        super(RotorAeroVSFPWithCCBlade, self).__init__()

    def configure(self):
        varspeed = True
        varpitch = False
        common_io_with_ccblade(self, varspeed, varpitch, self.cdf_type)
        common_configure_with_ccblade(self, varspeed, varpitch, self.cdf_type)



class RotorAeroFSVPWithCCBlade(Group):

    def __init__(self, cdf_type='weibull'):
        self.cdf_type = cdf_type
        super(RotorAeroFSVPWithCCBlade, self).__init__()

    def configure(self):
        varspeed = False
        varpitch = True
        common_io_with_ccblade(self, varspeed, varpitch, self.cdf_type)
        common_configure_with_ccblade(self, varspeed, varpitch, self.cdf_type)



class RotorAeroFSFPWithCCBlade(Group):

    def __init__(self, cdf_type='weibull'):
        self.cdf_type = cdf_type
        super(RotorAeroFSFPWithCCBlade, self).__init__()

    def configure(self):
        varspeed = False
        varpitch = False
        common_io_with_ccblade(self, varspeed, varpitch, self.cdf_type)
        common_configure_with_ccblade(self, varspeed, varpitch, self.cdf_type)


'''
