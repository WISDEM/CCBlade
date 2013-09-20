import unittest
import numpy as np
from math import pi

from ccblade import CCAirfoil, CCBlade


class TestGradients(unittest.TestCase):

    def setUp(self):

        # geometry
        self.Rhub = 1.5
        self.Rtip = 63.0

        self.r = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
                      28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
                      56.1667, 58.9000, 61.6333])
        self.chord = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
                          3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
        self.theta = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
                          6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])
        self.B = 3  # number of blades

        # atmosphere
        self.rho = 1.225
        self.mu = 1.81206e-5

        import os
        afinit = CCAirfoil.initFromAerodynFile  # just for shorthand
        basepath = '5MW_AFFiles' + os.path.sep

        # load all airfoils
        airfoil_types = [0]*8
        airfoil_types[0] = afinit(basepath + 'Cylinder1.dat')
        airfoil_types[1] = afinit(basepath + 'Cylinder2.dat')
        airfoil_types[2] = afinit(basepath + 'DU40_A17.dat')
        airfoil_types[3] = afinit(basepath + 'DU35_A17.dat')
        airfoil_types[4] = afinit(basepath + 'DU30_A17.dat')
        airfoil_types[5] = afinit(basepath + 'DU25_A17.dat')
        airfoil_types[6] = afinit(basepath + 'DU21_A17.dat')
        airfoil_types[7] = afinit(basepath + 'NACA64_A17.dat')

        # place at appropriate radial stations
        af_idx = [0, 0, 1, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7]

        self.af = [0]*len(self.r)
        for i in range(len(self.r)):
            self.af[i] = airfoil_types[af_idx[i]]


        self.tilt = -5.0
        self.precone = 2.5
        self.yaw = 0.0
        self.shearExp = 0.2
        self.hubHt = 80.0
        self.nSector = 8

        # create CCBlade object
        self.rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True)

        # set conditions
        self.Uinf = 10.0
        tsr = 7.55
        self.pitch = 0.0
        self.Omega = self.Uinf*tsr/self.Rtip * 30.0/pi  # convert to RPM
        self.azimuth = 90


        self.Np, self.Tp, self.dNp_dX, self.dTp_dX, self.dNp_dprecurve, self.dTp_dprecurve = \
            self.rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        self.P, self.T, self.Q, self.dP_ds, self.dT_ds, self.dQ_ds, self.dP_dv, self.dT_dv, \
            self.dQ_dv = self.rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        self.rotor.derivatives = False
        self.n = len(self.r)

        # X = [r, chord, theta, Rhub, Rtip, presweep, precone, tilt, hubHt]
        # scalars = [precone, tilt, hubHt, Rhub, Rtip, precurvetip, presweeptip]
        # vectors = [r, chord, theta, precurve, presweep]

    def test_dr(self):

        dNp_dr = self.dNp_dX[0, :]
        dTp_dr = self.dTp_dX[0, :]
        dT_dr = self.dT_dv[0, 0, :]
        dQ_dr = self.dQ_dv[0, 0, :]
        dP_dr = self.dP_dv[0, 0, :]
        dNp_dr_fd = np.zeros(self.n)
        dTp_dr_fd = np.zeros(self.n)
        dT_dr_fd = np.zeros(self.n)
        dQ_dr_fd = np.zeros(self.n)
        dP_dr_fd = np.zeros(self.n)

        for i in range(self.n):
            r = np.array(self.r)
            delta = 1e-6*r[i]
            r[i] += delta

            rotor = CCBlade(r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dNp_dr_fd[i] = (Npd[i] - self.Np[i]) / delta
            dTp_dr_fd[i] = (Tpd[i] - self.Tp[i]) / delta
            dT_dr_fd[i] = (Td - self.T) / delta
            dQ_dr_fd[i] = (Qd - self.Q) / delta
            dP_dr_fd[i] = (Pd - self.P) / delta


        np.testing.assert_allclose(dNp_dr_fd, dNp_dr, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dr_fd, dTp_dr, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dT_dr_fd, dT_dr, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dr_fd, dQ_dr, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dr_fd, dP_dr, rtol=3e-4, atol=1e-8)




    def test_dchord(self):

        dNp_dchord = self.dNp_dX[1, :]
        dTp_dchord = self.dTp_dX[1, :]
        dT_dchord = self.dT_dv[0, 1, :]
        dQ_dchord = self.dQ_dv[0, 1, :]
        dP_dchord = self.dP_dv[0, 1, :]
        dNp_dchord_fd = np.zeros(self.n)
        dTp_dchord_fd = np.zeros(self.n)
        dT_dchord_fd = np.zeros(self.n)
        dQ_dchord_fd = np.zeros(self.n)
        dP_dchord_fd = np.zeros(self.n)

        for i in range(self.n):
            chord = np.array(self.chord)
            delta = 1e-6*chord[i]
            chord[i] += delta

            rotor = CCBlade(self.r, chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dNp_dchord_fd[i] = (Npd[i] - self.Np[i]) / delta
            dTp_dchord_fd[i] = (Tpd[i] - self.Tp[i]) / delta
            dT_dchord_fd[i] = (Td - self.T) / delta
            dQ_dchord_fd[i] = (Qd - self.Q) / delta
            dP_dchord_fd[i] = (Pd - self.P) / delta


        np.testing.assert_allclose(dNp_dchord_fd, dNp_dchord, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dchord_fd, dTp_dchord, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dT_dchord_fd, dT_dchord, rtol=5e-6, atol=1e-8)
        np.testing.assert_allclose(dQ_dchord_fd, dQ_dchord, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dchord_fd, dP_dchord, rtol=7e-5, atol=1e-8)







    def test_dtheta(self):

        dNp_dtheta = self.dNp_dX[2, :]
        dTp_dtheta = self.dTp_dX[2, :]
        dT_dtheta = self.dT_dv[0, 2, :]
        dQ_dtheta = self.dQ_dv[0, 2, :]
        dP_dtheta = self.dP_dv[0, 2, :]
        dNp_dtheta_fd = np.zeros(self.n)
        dTp_dtheta_fd = np.zeros(self.n)
        dT_dtheta_fd = np.zeros(self.n)
        dQ_dtheta_fd = np.zeros(self.n)
        dP_dtheta_fd = np.zeros(self.n)

        for i in range(self.n):
            theta = np.array(self.theta)
            delta = 1e-6*theta[i]
            theta[i] += delta

            rotor = CCBlade(self.r, self.chord, theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dNp_dtheta_fd[i] = (Npd[i] - self.Np[i]) / delta
            dTp_dtheta_fd[i] = (Tpd[i] - self.Tp[i]) / delta
            dT_dtheta_fd[i] = (Td - self.T) / delta
            dQ_dtheta_fd[i] = (Qd - self.Q) / delta
            dP_dtheta_fd[i] = (Pd - self.P) / delta


        np.testing.assert_allclose(dNp_dtheta_fd, dNp_dtheta, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dtheta_fd, dTp_dtheta, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dT_dtheta_fd, dT_dtheta, rtol=5e-6, atol=1e-8)
        np.testing.assert_allclose(dQ_dtheta_fd, dQ_dtheta, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dtheta_fd, dP_dtheta, rtol=7e-5, atol=1e-8)


    def test_dRhub(self):

        dNp_dRhub = self.dNp_dX[3, :]
        dTp_dRhub = self.dTp_dX[3, :]
        dT_dRhub = self.dT_ds[0, 3]
        dQ_dRhub = self.dQ_ds[0, 3]
        dP_dRhub = self.dP_ds[0, 3]


        Rhub = float(self.Rhub)
        delta = 1e-6*Rhub
        Rhub += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dNp_dRhub_fd = (Npd - self.Np) / delta
        dTp_dRhub_fd = (Tpd - self.Tp) / delta
        dT_dRhub_fd = (Td - self.T) / delta
        dQ_dRhub_fd = (Qd - self.Q) / delta
        dP_dRhub_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dNp_dRhub_fd, dNp_dRhub, rtol=1e-5, atol=1e-7)
        np.testing.assert_allclose(dTp_dRhub_fd, dTp_dRhub, rtol=1e-4, atol=1e-7)
        np.testing.assert_allclose(dT_dRhub_fd, dT_dRhub, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dRhub_fd, dQ_dRhub, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dRhub_fd, dP_dRhub, rtol=5e-5, atol=1e-8)


    def test_dRtip(self):

        dNp_dRtip = self.dNp_dX[4, :]
        dTp_dRtip = self.dTp_dX[4, :]
        dT_dRtip = self.dT_ds[0, 4]
        dQ_dRtip = self.dQ_ds[0, 4]
        dP_dRtip = self.dP_ds[0, 4]


        Rtip = float(self.Rtip)
        delta = 1e-6*Rtip
        Rtip += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dNp_dRtip_fd = (Npd - self.Np) / delta
        dTp_dRtip_fd = (Tpd - self.Tp) / delta
        dT_dRtip_fd = (Td - self.T) / delta
        dQ_dRtip_fd = (Qd - self.Q) / delta
        dP_dRtip_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dNp_dRtip_fd, dNp_dRtip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dRtip_fd, dTp_dRtip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dT_dRtip_fd, dT_dRtip, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dRtip_fd, dQ_dRtip, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dRtip_fd, dP_dRtip, rtol=5e-5, atol=1e-8)


    def test_dprecone(self):

        dNp_dprecone = self.dNp_dX[6, :]
        dTp_dprecone = self.dTp_dX[6, :]
        dT_dprecone = self.dT_ds[0, 0]
        dQ_dprecone = self.dQ_ds[0, 0]
        dP_dprecone = self.dP_ds[0, 0]


        precone = float(self.precone)
        delta = 1e-6*precone
        precone += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dNp_dprecone_fd = (Npd - self.Np) / delta
        dTp_dprecone_fd = (Tpd - self.Tp) / delta
        dT_dprecone_fd = (Td - self.T) / delta
        dQ_dprecone_fd = (Qd - self.Q) / delta
        dP_dprecone_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dNp_dprecone_fd, dNp_dprecone, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dprecone_fd, dTp_dprecone, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dT_dprecone_fd, dT_dprecone, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dprecone_fd, dQ_dprecone, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dprecone_fd, dP_dprecone, rtol=5e-5, atol=1e-8)


    def test_dtilt(self):

        dNp_dtilt = self.dNp_dX[7, :]
        dTp_dtilt = self.dTp_dX[7, :]
        dT_dtilt = self.dT_ds[0, 1]
        dQ_dtilt = self.dQ_ds[0, 1]
        dP_dtilt = self.dP_ds[0, 1]


        tilt = float(self.tilt)
        delta = 1e-6*tilt
        tilt += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dNp_dtilt_fd = (Npd - self.Np) / delta
        dTp_dtilt_fd = (Tpd - self.Tp) / delta
        dT_dtilt_fd = (Td - self.T) / delta
        dQ_dtilt_fd = (Qd - self.Q) / delta
        dP_dtilt_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dNp_dtilt_fd, dNp_dtilt, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dtilt_fd, dTp_dtilt, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dT_dtilt_fd, dT_dtilt, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dtilt_fd, dQ_dtilt, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dtilt_fd, dP_dtilt, rtol=5e-5, atol=1e-8)

    def test_dhubht(self):

        dNp_dhubht = self.dNp_dX[8, :]
        dTp_dhubht = self.dTp_dX[8, :]
        dT_dhubht = self.dT_ds[0, 2]
        dQ_dhubht = self.dQ_ds[0, 2]
        dP_dhubht = self.dP_ds[0, 2]

        hubht = float(self.hubHt)
        delta = 1e-6*hubht
        hubht += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            hubht, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dNp_dhubht_fd = (Npd - self.Np) / delta
        dTp_dhubht_fd = (Tpd - self.Tp) / delta
        dT_dhubht_fd = (Td - self.T) / delta
        dQ_dhubht_fd = (Qd - self.Q) / delta
        dP_dhubht_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dNp_dhubht_fd, dNp_dhubht, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dTp_dhubht_fd, dTp_dhubht, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dT_dhubht_fd, dT_dhubht, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dhubht_fd, dQ_dhubht, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dhubht_fd, dP_dhubht, rtol=5e-5, atol=1e-8)

    def test_dprecurve(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        Np, Tp, dNp_dX, dTp_dX, dNp_dprecurve, dTp_dprecurve = \
            rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        P, T, Q, dP_ds, dT_ds, dQ_ds, dP_dv, dT_dv, dQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dprecurve = dT_dv[0, 3, :]
        dQ_dprecurve = dQ_dv[0, 3, :]
        dP_dprecurve = dP_dv[0, 3, :]


        dNp_dprecurve_fd = np.zeros((self.n, self.n))
        dTp_dprecurve_fd = np.zeros((self.n, self.n))
        dT_dprecurve_fd = np.zeros(self.n)
        dQ_dprecurve_fd = np.zeros(self.n)
        dP_dprecurve_fd = np.zeros(self.n)
        for i in range(self.n):
            pc = np.array(precurve)
            delta = 1e-6*pc[i]
            pc[i] += delta

            rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False, precurve=pc, precurveTip=precurveTip)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dNp_dprecurve_fd[i, :] = (Npd - Np) / delta
            dTp_dprecurve_fd[i, :] = (Tpd - Tp) / delta
            dT_dprecurve_fd[i] = (Td - T) / delta
            dQ_dprecurve_fd[i] = (Qd - Q) / delta
            dP_dprecurve_fd[i] = (Pd - P) / delta

        np.testing.assert_allclose(dNp_dprecurve_fd, dNp_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dprecurve_fd, dTp_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dT_dprecurve_fd, dT_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dprecurve_fd, dQ_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dprecurve_fd, dP_dprecurve, rtol=3e-4, atol=1e-8)


    def test_dpresweep(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        Np, Tp, dNp_dX, dTp_dX, dNp_dprecurve, dTp_dprecurve = \
            rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        P, T, Q, dP_ds, dT_ds, dQ_ds, dP_dv, dT_dv, dQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dNp_dpresweep = dNp_dX[5, :]
        dTp_dpresweep = dTp_dX[5, :]
        dT_dpresweep = dT_dv[0, 4, :]
        dQ_dpresweep = dQ_dv[0, 4, :]
        dP_dpresweep = dP_dv[0, 4, :]


        dNp_dpresweep_fd = np.zeros(self.n)
        dTp_dpresweep_fd = np.zeros(self.n)
        dT_dpresweep_fd = np.zeros(self.n)
        dQ_dpresweep_fd = np.zeros(self.n)
        dP_dpresweep_fd = np.zeros(self.n)
        for i in range(self.n):
            ps = np.array(presweep)
            delta = 1e-6*ps[i]
            ps[i] += delta

            rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False, presweep=ps, presweepTip=presweepTip)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dNp_dpresweep_fd[i] = (Npd[i] - Np[i]) / delta
            dTp_dpresweep_fd[i] = (Tpd[i] - Tp[i]) / delta
            dT_dpresweep_fd[i] = (Td - T) / delta
            dQ_dpresweep_fd[i] = (Qd - Q) / delta
            dP_dpresweep_fd[i] = (Pd - P) / delta


        np.testing.assert_allclose(dNp_dpresweep_fd, dNp_dpresweep, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dTp_dpresweep_fd, dTp_dpresweep, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dT_dpresweep_fd, dT_dpresweep, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dpresweep_fd, dQ_dpresweep, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dpresweep_fd, dP_dpresweep, rtol=3e-4, atol=1e-8)


    def test_dprecurveTip(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        P, T, Q, dP_ds, dT_ds, dQ_ds, dP_dv, dT_dv, dQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dprecurveTip = dT_ds[0, 5]
        dQ_dprecurveTip = dQ_ds[0, 5]
        dP_dprecurveTip = dP_ds[0, 5]

        pct = float(precurveTip)
        delta = 1e-6*pct
        pct += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False, precurve=precurve, precurveTip=pct)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dprecurveTip_fd = (Td - T) / delta
        dQ_dprecurveTip_fd = (Qd - Q) / delta
        dP_dprecurveTip_fd = (Pd - P) / delta

        np.testing.assert_allclose(dT_dprecurveTip_fd, dT_dprecurveTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dprecurveTip_fd, dQ_dprecurveTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dprecurveTip_fd, dP_dprecurveTip, rtol=1e-4, atol=1e-8)


    def test_dpresweepTip(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        P, T, Q, dP_ds, dT_ds, dQ_ds, dP_dv, dT_dv, dQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dpresweepTip = dT_ds[0, 6]
        dQ_dpresweepTip = dQ_ds[0, 6]
        dP_dpresweepTip = dP_ds[0, 6]

        pst = float(presweepTip)
        delta = 1e-6*pst
        pst += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False, presweep=presweep, presweepTip=pst)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dpresweepTip_fd = (Td - T) / delta
        dQ_dpresweepTip_fd = (Qd - Q) / delta
        dP_dpresweepTip_fd = (Pd - P) / delta

        np.testing.assert_allclose(dT_dpresweepTip_fd, dT_dpresweepTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dpresweepTip_fd, dQ_dpresweepTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dpresweepTip_fd, dP_dpresweepTip, rtol=1e-4, atol=1e-8)


if __name__ == '__main__':
    unittest.main()
