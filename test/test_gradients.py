"""
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import unittest
import numpy as np
from math import pi
from os import path

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

        afinit = CCAirfoil.initFromAerodynFile  # just for shorthand
        basepath = path.join(path.dirname(path.realpath(__file__)), '5MW_AFFiles') + path.sep

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

        self.CP, self.CT, self.CQ, self.dCP_ds, self.dCT_ds, self.dCQ_ds, self.dCP_dv, self.dCT_dv, \
            self.dCQ_dv = self.rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        self.rotor.derivatives = False
        self.n = len(self.r)

        # X = [r, chord, theta, Rhub, Rtip, presweep, precone, tilt, hubHt]
        # scalars = [precone, tilt, hubHt, Rhub, Rtip, precurvetip, presweeptip]
        # vectors = [r, chord, theta, precurve, presweep]

    def test_dr1(self):

        dNp_dr = self.dNp_dX[0, :]
        dTp_dr = self.dTp_dX[0, :]
        dNp_dr_fd = np.zeros(self.n)
        dTp_dr_fd = np.zeros(self.n)

        for i in range(self.n):
            r = np.array(self.r)
            delta = 1e-6*r[i]
            r[i] += delta

            rotor = CCBlade(r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

            dNp_dr_fd[i] = (Npd[i] - self.Np[i]) / delta
            dTp_dr_fd[i] = (Tpd[i] - self.Tp[i]) / delta


        np.testing.assert_allclose(dNp_dr_fd, dNp_dr, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dr_fd, dTp_dr, rtol=1e-4, atol=1e-8)


    def test_dr2(self):

        dT_dr = self.dT_dv[0, 0, :]
        dQ_dr = self.dQ_dv[0, 0, :]
        dP_dr = self.dP_dv[0, 0, :]
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

            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dT_dr_fd[i] = (Td - self.T) / delta
            dQ_dr_fd[i] = (Qd - self.Q) / delta
            dP_dr_fd[i] = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dr_fd, dT_dr, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dr_fd, dQ_dr, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dr_fd, dP_dr, rtol=3e-4, atol=1e-8)


    def test_dr3(self):

        dCT_dr = self.dCT_dv[0, 0, :]
        dCQ_dr = self.dCQ_dv[0, 0, :]
        dCP_dr = self.dCP_dv[0, 0, :]
        dCT_dr_fd = np.zeros(self.n)
        dCQ_dr_fd = np.zeros(self.n)
        dCP_dr_fd = np.zeros(self.n)

        for i in range(self.n):
            r = np.array(self.r)
            delta = 1e-6*r[i]
            r[i] += delta

            rotor = CCBlade(r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

            dCT_dr_fd[i] = (CTd - self.CT) / delta
            dCQ_dr_fd[i] = (CQd - self.CQ) / delta
            dCP_dr_fd[i] = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dr_fd, dCT_dr, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dCQ_dr_fd, dCQ_dr, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dCP_dr_fd, dCP_dr, rtol=3e-4, atol=1e-8)



    def test_dchord1(self):

        dNp_dchord = self.dNp_dX[1, :]
        dTp_dchord = self.dTp_dX[1, :]
        dNp_dchord_fd = np.zeros(self.n)
        dTp_dchord_fd = np.zeros(self.n)

        for i in range(self.n):
            chord = np.array(self.chord)
            delta = 1e-6*chord[i]
            chord[i] += delta

            rotor = CCBlade(self.r, chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

            dNp_dchord_fd[i] = (Npd[i] - self.Np[i]) / delta
            dTp_dchord_fd[i] = (Tpd[i] - self.Tp[i]) / delta


        np.testing.assert_allclose(dNp_dchord_fd, dNp_dchord, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dchord_fd, dTp_dchord, rtol=5e-5, atol=1e-8)



    def test_dchord2(self):

        dT_dchord = self.dT_dv[0, 1, :]
        dQ_dchord = self.dQ_dv[0, 1, :]
        dP_dchord = self.dP_dv[0, 1, :]
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

            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dT_dchord_fd[i] = (Td - self.T) / delta
            dQ_dchord_fd[i] = (Qd - self.Q) / delta
            dP_dchord_fd[i] = (Pd - self.P) / delta


        np.testing.assert_allclose(dT_dchord_fd, dT_dchord, rtol=5e-6, atol=1e-8)
        np.testing.assert_allclose(dQ_dchord_fd, dQ_dchord, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dchord_fd, dP_dchord, rtol=7e-5, atol=1e-8)

    def test_dchord3(self):

        dCT_dchord = self.dCT_dv[0, 1, :]
        dCQ_dchord = self.dCQ_dv[0, 1, :]
        dCP_dchord = self.dCP_dv[0, 1, :]
        dCT_dchord_fd = np.zeros(self.n)
        dCQ_dchord_fd = np.zeros(self.n)
        dCP_dchord_fd = np.zeros(self.n)

        for i in range(self.n):
            chord = np.array(self.chord)
            delta = 1e-6*chord[i]
            chord[i] += delta

            rotor = CCBlade(self.r, chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

            dCT_dchord_fd[i] = (CTd - self.CT) / delta
            dCQ_dchord_fd[i] = (CQd - self.CQ) / delta
            dCP_dchord_fd[i] = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dchord_fd, dCT_dchord, rtol=5e-6, atol=1e-8)
        np.testing.assert_allclose(dCQ_dchord_fd, dCQ_dchord, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dchord_fd, dCP_dchord, rtol=7e-5, atol=1e-8)




    def test_dtheta1(self):

        dNp_dtheta = self.dNp_dX[2, :]
        dTp_dtheta = self.dTp_dX[2, :]
        dNp_dtheta_fd = np.zeros(self.n)
        dTp_dtheta_fd = np.zeros(self.n)

        for i in range(self.n):
            theta = np.array(self.theta)
            delta = 1e-6*theta[i]
            theta[i] += delta

            rotor = CCBlade(self.r, self.chord, theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

            dNp_dtheta_fd[i] = (Npd[i] - self.Np[i]) / delta
            dTp_dtheta_fd[i] = (Tpd[i] - self.Tp[i]) / delta


        np.testing.assert_allclose(dNp_dtheta_fd, dNp_dtheta, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dtheta_fd, dTp_dtheta, rtol=1e-4, atol=1e-8)


    def test_dtheta2(self):

        dT_dtheta = self.dT_dv[0, 2, :]
        dQ_dtheta = self.dQ_dv[0, 2, :]
        dP_dtheta = self.dP_dv[0, 2, :]
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

            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dT_dtheta_fd[i] = (Td - self.T) / delta
            dQ_dtheta_fd[i] = (Qd - self.Q) / delta
            dP_dtheta_fd[i] = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dtheta_fd, dT_dtheta, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dtheta_fd, dQ_dtheta, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dtheta_fd, dP_dtheta, rtol=7e-5, atol=1e-8)



    def test_dtheta3(self):

        dCT_dtheta = self.dCT_dv[0, 2, :]
        dCQ_dtheta = self.dCQ_dv[0, 2, :]
        dCP_dtheta = self.dCP_dv[0, 2, :]
        dCT_dtheta_fd = np.zeros(self.n)
        dCQ_dtheta_fd = np.zeros(self.n)
        dCP_dtheta_fd = np.zeros(self.n)

        for i in range(self.n):
            theta = np.array(self.theta)
            delta = 1e-6*theta[i]
            theta[i] += delta

            rotor = CCBlade(self.r, self.chord, theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False)

            CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

            dCT_dtheta_fd[i] = (CTd - self.CT) / delta
            dCQ_dtheta_fd[i] = (CQd - self.CQ) / delta
            dCP_dtheta_fd[i] = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dtheta_fd, dCT_dtheta, rtol=5e-6, atol=1e-8)
        np.testing.assert_allclose(dCQ_dtheta_fd, dCQ_dtheta, rtol=7e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dtheta_fd, dCP_dtheta, rtol=7e-5, atol=1e-8)



    def test_dRhub1(self):

        dNp_dRhub = self.dNp_dX[3, :]
        dTp_dRhub = self.dTp_dX[3, :]

        Rhub = float(self.Rhub)
        delta = 1e-6*Rhub
        Rhub += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        dNp_dRhub_fd = (Npd - self.Np) / delta
        dTp_dRhub_fd = (Tpd - self.Tp) / delta

        np.testing.assert_allclose(dNp_dRhub_fd, dNp_dRhub, rtol=1e-5, atol=1e-7)
        np.testing.assert_allclose(dTp_dRhub_fd, dTp_dRhub, rtol=1e-4, atol=1e-7)


    def test_dRhub2(self):

        dT_dRhub = self.dT_ds[0, 3]
        dQ_dRhub = self.dQ_ds[0, 3]
        dP_dRhub = self.dP_ds[0, 3]


        Rhub = float(self.Rhub)
        delta = 1e-6*Rhub
        Rhub += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dRhub_fd = (Td - self.T) / delta
        dQ_dRhub_fd = (Qd - self.Q) / delta
        dP_dRhub_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dRhub_fd, dT_dRhub, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dRhub_fd, dQ_dRhub, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dRhub_fd, dP_dRhub, rtol=5e-5, atol=1e-8)


    def test_dRhub3(self):

        dCT_dRhub = self.dCT_ds[0, 3]
        dCQ_dRhub = self.dCQ_ds[0, 3]
        dCP_dRhub = self.dCP_ds[0, 3]


        Rhub = float(self.Rhub)
        delta = 1e-6*Rhub
        Rhub += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dRhub_fd = (CTd - self.CT) / delta
        dCQ_dRhub_fd = (CQd - self.CQ) / delta
        dCP_dRhub_fd = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dRhub_fd, dCT_dRhub, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dCQ_dRhub_fd, dCQ_dRhub, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dRhub_fd, dCP_dRhub, rtol=5e-5, atol=1e-8)


    def test_dRtip1(self):

        dNp_dRtip = self.dNp_dX[4, :]
        dTp_dRtip = self.dTp_dX[4, :]


        Rtip = float(self.Rtip)
        delta = 1e-6*Rtip
        Rtip += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        dNp_dRtip_fd = (Npd - self.Np) / delta
        dTp_dRtip_fd = (Tpd - self.Tp) / delta

        np.testing.assert_allclose(dNp_dRtip_fd, dNp_dRtip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dRtip_fd, dTp_dRtip, rtol=1e-4, atol=1e-8)


    def test_dRtip2(self):

        dT_dRtip = self.dT_ds[0, 4]
        dQ_dRtip = self.dQ_ds[0, 4]
        dP_dRtip = self.dP_ds[0, 4]

        Rtip = float(self.Rtip)
        delta = 1e-6*Rtip
        Rtip += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dRtip_fd = (Td - self.T) / delta
        dQ_dRtip_fd = (Qd - self.Q) / delta
        dP_dRtip_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dRtip_fd, dT_dRtip, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dRtip_fd, dQ_dRtip, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dRtip_fd, dP_dRtip, rtol=5e-5, atol=1e-8)


    def test_dRtip3(self):

        dCT_dRtip = self.dCT_ds[0, 4]
        dCQ_dRtip = self.dCQ_ds[0, 4]
        dCP_dRtip = self.dCP_ds[0, 4]


        Rtip = float(self.Rtip)
        delta = 1e-6*Rtip
        Rtip += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dRtip_fd = (CTd - self.CT) / delta
        dCQ_dRtip_fd = (CQd - self.CQ) / delta
        dCP_dRtip_fd = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dRtip_fd, dCT_dRtip, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dCQ_dRtip_fd, dCQ_dRtip, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dRtip_fd, dCP_dRtip, rtol=5e-5, atol=1e-8)


    def test_dprecone1(self):

        dNp_dprecone = self.dNp_dX[6, :]
        dTp_dprecone = self.dTp_dX[6, :]

        precone = float(self.precone)
        delta = 1e-6*precone
        precone += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        dNp_dprecone_fd = (Npd - self.Np) / delta
        dTp_dprecone_fd = (Tpd - self.Tp) / delta

        np.testing.assert_allclose(dNp_dprecone_fd, dNp_dprecone, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dprecone_fd, dTp_dprecone, rtol=1e-6, atol=1e-8)



    def test_dprecone2(self):

        dT_dprecone = self.dT_ds[0, 0]
        dQ_dprecone = self.dQ_ds[0, 0]
        dP_dprecone = self.dP_ds[0, 0]


        precone = float(self.precone)
        delta = 1e-6*precone
        precone += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dprecone_fd = (Td - self.T) / delta
        dQ_dprecone_fd = (Qd - self.Q) / delta
        dP_dprecone_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dprecone_fd, dT_dprecone, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dprecone_fd, dQ_dprecone, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dprecone_fd, dP_dprecone, rtol=5e-5, atol=1e-8)


    def test_dprecone3(self):

        dCT_dprecone = self.dCT_ds[0, 0]
        dCQ_dprecone = self.dCQ_ds[0, 0]
        dCP_dprecone = self.dCP_ds[0, 0]

        precone = float(self.precone)
        delta = 1e-6*precone
        precone += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dprecone_fd = (CTd - self.CT) / delta
        dCQ_dprecone_fd = (CQd - self.CQ) / delta
        dCP_dprecone_fd = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dprecone_fd, dCT_dprecone, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dCQ_dprecone_fd, dCQ_dprecone, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dprecone_fd, dCP_dprecone, rtol=5e-5, atol=1e-8)


    def test_dtilt1(self):

        dNp_dtilt = self.dNp_dX[7, :]
        dTp_dtilt = self.dTp_dX[7, :]

        tilt = float(self.tilt)
        delta = 1e-6*tilt
        tilt += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        dNp_dtilt_fd = (Npd - self.Np) / delta
        dTp_dtilt_fd = (Tpd - self.Tp) / delta

        np.testing.assert_allclose(dNp_dtilt_fd, dNp_dtilt, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(dTp_dtilt_fd, dTp_dtilt, rtol=1e-5, atol=1e-8)


    def test_dtilt2(self):

        dT_dtilt = self.dT_ds[0, 1]
        dQ_dtilt = self.dQ_ds[0, 1]
        dP_dtilt = self.dP_ds[0, 1]


        tilt = float(self.tilt)
        delta = 1e-6*tilt
        tilt += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dtilt_fd = (Td - self.T) / delta
        dQ_dtilt_fd = (Qd - self.Q) / delta
        dP_dtilt_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dtilt_fd, dT_dtilt, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dtilt_fd, dQ_dtilt, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dtilt_fd, dP_dtilt, rtol=5e-5, atol=1e-8)


    def test_dtilt3(self):

        dCT_dtilt = self.dCT_ds[0, 1]
        dCQ_dtilt = self.dCQ_ds[0, 1]
        dCP_dtilt = self.dCP_ds[0, 1]


        tilt = float(self.tilt)
        delta = 1e-6*tilt
        tilt += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dtilt_fd = (CTd - self.CT) / delta
        dCQ_dtilt_fd = (CQd - self.CQ) / delta
        dCP_dtilt_fd = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dtilt_fd, dCT_dtilt, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dCQ_dtilt_fd, dCQ_dtilt, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dtilt_fd, dCP_dtilt, rtol=5e-5, atol=1e-8)


    def test_dhubht1(self):

        dNp_dhubht = self.dNp_dX[8, :]
        dTp_dhubht = self.dTp_dX[8, :]

        hubht = float(self.hubHt)
        delta = 1e-6*hubht
        hubht += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            hubht, self.nSector, derivatives=False)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        dNp_dhubht_fd = (Npd - self.Np) / delta
        dTp_dhubht_fd = (Tpd - self.Tp) / delta

        np.testing.assert_allclose(dNp_dhubht_fd, dNp_dhubht, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dTp_dhubht_fd, dTp_dhubht, rtol=1e-5, atol=1e-8)


    def test_dhubht2(self):

        dT_dhubht = self.dT_ds[0, 2]
        dQ_dhubht = self.dQ_ds[0, 2]
        dP_dhubht = self.dP_ds[0, 2]

        hubht = float(self.hubHt)
        delta = 1e-6*hubht
        hubht += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            hubht, self.nSector, derivatives=False)

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dhubht_fd = (Td - self.T) / delta
        dQ_dhubht_fd = (Qd - self.Q) / delta
        dP_dhubht_fd = (Pd - self.P) / delta

        np.testing.assert_allclose(dT_dhubht_fd, dT_dhubht, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dQ_dhubht_fd, dQ_dhubht, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dP_dhubht_fd, dP_dhubht, rtol=5e-5, atol=1e-8)



    def test_dhubht3(self):

        dCT_dhubht = self.dCT_ds[0, 2]
        dCQ_dhubht = self.dCQ_ds[0, 2]
        dCP_dhubht = self.dCP_ds[0, 2]

        hubht = float(self.hubHt)
        delta = 1e-6*hubht
        hubht += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, self.precone, self.tilt, self.yaw, self.shearExp,
            hubht, self.nSector, derivatives=False)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dhubht_fd = (CTd - self.CT) / delta
        dCQ_dhubht_fd = (CQd - self.CQ) / delta
        dCP_dhubht_fd = (CPd - self.CP) / delta

        np.testing.assert_allclose(dCT_dhubht_fd, dCT_dhubht, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dCQ_dhubht_fd, dCQ_dhubht, rtol=5e-5, atol=1e-8)
        np.testing.assert_allclose(dCP_dhubht_fd, dCP_dhubht, rtol=5e-5, atol=1e-8)


    def test_dprecurve1(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        Np, Tp, dNp_dX, dTp_dX, dNp_dprecurve, dTp_dprecurve = \
            rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        dNp_dprecurve_fd = np.zeros((self.n, self.n))
        dTp_dprecurve_fd = np.zeros((self.n, self.n))
        for i in range(self.n):
            pc = np.array(precurve)
            delta = 1e-6*pc[i]
            pc[i] += delta

            rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False, precurve=pc, precurveTip=precurveTip)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

            dNp_dprecurve_fd[i, :] = (Npd - Np) / delta
            dTp_dprecurve_fd[i, :] = (Tpd - Tp) / delta

        np.testing.assert_allclose(dNp_dprecurve_fd, dNp_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dprecurve_fd, dTp_dprecurve, rtol=3e-4, atol=1e-8)


    def test_dprecurve2(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        P, T, Q, dP_ds, dT_ds, dQ_ds, dP_dv, dT_dv, dQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dprecurve = dT_dv[0, 3, :]
        dQ_dprecurve = dQ_dv[0, 3, :]
        dP_dprecurve = dP_dv[0, 3, :]

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

            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dT_dprecurve_fd[i] = (Td - T) / delta
            dQ_dprecurve_fd[i] = (Qd - Q) / delta
            dP_dprecurve_fd[i] = (Pd - P) / delta

        np.testing.assert_allclose(dT_dprecurve_fd, dT_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dprecurve_fd, dQ_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dprecurve_fd, dP_dprecurve, rtol=3e-4, atol=1e-8)


    def test_dprecurve3(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        CP, CT, CQ, dCP_ds, dCT_ds, dCQ_ds, dCP_dv, dCT_dv, dCQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dprecurve = dCT_dv[0, 3, :]
        dCQ_dprecurve = dCQ_dv[0, 3, :]
        dCP_dprecurve = dCP_dv[0, 3, :]


        dCT_dprecurve_fd = np.zeros(self.n)
        dCQ_dprecurve_fd = np.zeros(self.n)
        dCP_dprecurve_fd = np.zeros(self.n)
        for i in range(self.n):
            pc = np.array(precurve)
            delta = 1e-6*pc[i]
            pc[i] += delta

            rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False, precurve=pc, precurveTip=precurveTip)

            CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

            dCT_dprecurve_fd[i] = (CTd - CT) / delta
            dCQ_dprecurve_fd[i] = (CQd - CQ) / delta
            dCP_dprecurve_fd[i] = (CPd - CP) / delta

        np.testing.assert_allclose(dCT_dprecurve_fd, dCT_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dCQ_dprecurve_fd, dCQ_dprecurve, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dCP_dprecurve_fd, dCP_dprecurve, rtol=3e-4, atol=1e-8)


    def test_dpresweep1(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        Np, Tp, dNp_dX, dTp_dX, dNp_dprecurve, dTp_dprecurve = \
            rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        dNp_dpresweep = dNp_dX[5, :]
        dTp_dpresweep = dTp_dX[5, :]

        dNp_dpresweep_fd = np.zeros(self.n)
        dTp_dpresweep_fd = np.zeros(self.n)
        for i in range(self.n):
            ps = np.array(presweep)
            delta = 1e-6*ps[i]
            ps[i] += delta

            rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False, presweep=ps, presweepTip=presweepTip)

            Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

            dNp_dpresweep_fd[i] = (Npd[i] - Np[i]) / delta
            dTp_dpresweep_fd[i] = (Tpd[i] - Tp[i]) / delta

        np.testing.assert_allclose(dNp_dpresweep_fd, dNp_dpresweep, rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(dTp_dpresweep_fd, dTp_dpresweep, rtol=1e-5, atol=1e-8)


    def test_dpresweep2(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        P, T, Q, dP_ds, dT_ds, dQ_ds, dP_dv, dT_dv, dQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dpresweep = dT_dv[0, 4, :]
        dQ_dpresweep = dQ_dv[0, 4, :]
        dP_dpresweep = dP_dv[0, 4, :]


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

            Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

            dT_dpresweep_fd[i] = (Td - T) / delta
            dQ_dpresweep_fd[i] = (Qd - Q) / delta
            dP_dpresweep_fd[i] = (Pd - P) / delta


        np.testing.assert_allclose(dT_dpresweep_fd, dT_dpresweep, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dpresweep_fd, dQ_dpresweep, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dpresweep_fd, dP_dpresweep, rtol=3e-4, atol=1e-8)




    def test_dpresweep3(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        CP, CT, CQ, dCP_ds, dCT_ds, dCQ_ds, dCP_dv, dCT_dv, dCQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dpresweep = dCT_dv[0, 4, :]
        dCQ_dpresweep = dCQ_dv[0, 4, :]
        dCP_dpresweep = dCP_dv[0, 4, :]

        dCT_dpresweep_fd = np.zeros(self.n)
        dCQ_dpresweep_fd = np.zeros(self.n)
        dCP_dpresweep_fd = np.zeros(self.n)
        for i in range(self.n):
            ps = np.array(presweep)
            delta = 1e-6*ps[i]
            ps[i] += delta

            rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
                self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
                self.hubHt, self.nSector, derivatives=False, presweep=ps, presweepTip=presweepTip)

            CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

            dCT_dpresweep_fd[i] = (CTd - CT) / delta
            dCQ_dpresweep_fd[i] = (CQd - CQ) / delta
            dCP_dpresweep_fd[i] = (CPd - CP) / delta


        np.testing.assert_allclose(dCT_dpresweep_fd, dCT_dpresweep, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dCQ_dpresweep_fd, dCQ_dpresweep, rtol=3e-4, atol=1e-8)
        np.testing.assert_allclose(dCP_dpresweep_fd, dCP_dpresweep, rtol=3e-4, atol=1e-8)



    def test_dprecurveTip1(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        Np, Tp, dNp_dX, dTp_dX, dNp_dprecurve, dTp_dprecurve = \
            rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        pct = float(precurveTip)
        delta = 1e-6*pct
        pct += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False, precurve=precurve, precurveTip=pct)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        dNp_dprecurveTip_fd = (Npd - Np) / delta
        dTp_dprecurveTip_fd = (Tpd - Tp) / delta

        np.testing.assert_allclose(dNp_dprecurveTip_fd, 0.0, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dprecurveTip_fd, 0.0, rtol=1e-4, atol=1e-8)


    def test_dprecurveTip2(self):

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

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dprecurveTip_fd = (Td - T) / delta
        dQ_dprecurveTip_fd = (Qd - Q) / delta
        dP_dprecurveTip_fd = (Pd - P) / delta

        np.testing.assert_allclose(dT_dprecurveTip_fd, dT_dprecurveTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dprecurveTip_fd, dQ_dprecurveTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dprecurveTip_fd, dP_dprecurveTip, rtol=1e-4, atol=1e-8)



    def test_dprecurveTip3(self):

        precurve = np.linspace(1, 10, self.n)
        precurveTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, precurve=precurve, precurveTip=precurveTip)

        CP, CT, CQ, dCP_ds, dCT_ds, dCQ_ds, dCP_dv, dCT_dv, dCQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dprecurveTip = dCT_ds[0, 5]
        dCQ_dprecurveTip = dCQ_ds[0, 5]
        dCP_dprecurveTip = dCP_ds[0, 5]

        pct = float(precurveTip)
        delta = 1e-6*pct
        pct += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False, precurve=precurve, precurveTip=pct)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dprecurveTip_fd = (CTd - CT) / delta
        dCQ_dprecurveTip_fd = (CQd - CQ) / delta
        dCP_dprecurveTip_fd = (CPd - CP) / delta

        np.testing.assert_allclose(dCT_dprecurveTip_fd, dCT_dprecurveTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dCQ_dprecurveTip_fd, dCQ_dprecurveTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dCP_dprecurveTip_fd, dCP_dprecurveTip, rtol=1e-4, atol=1e-8)


    def test_dpresweepTip1(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        Np, Tp, dNp_dX, dTp_dX, dNp_dpresweep, dTp_dpresweep = \
            rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)

        pst = float(presweepTip)
        delta = 1e-6*pst
        pst += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False, presweep=presweep, presweepTip=pst)

        Npd, Tpd = rotor.distributedAeroLoads(self.Uinf, self.Omega, self.pitch, self.azimuth)
        dNp_dpresweepTip_fd = (Npd - Np) / delta
        dTp_dpresweepTip_fd = (Tpd - Tp) / delta

        np.testing.assert_allclose(dNp_dpresweepTip_fd, 0.0, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dTp_dpresweepTip_fd, 0.0, rtol=1e-4, atol=1e-8)


    def test_dpresweepTip2(self):

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

        Pd, Td, Qd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=False)

        dT_dpresweepTip_fd = (Td - T) / delta
        dQ_dpresweepTip_fd = (Qd - Q) / delta
        dP_dpresweepTip_fd = (Pd - P) / delta

        np.testing.assert_allclose(dT_dpresweepTip_fd, dT_dpresweepTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dQ_dpresweepTip_fd, dQ_dpresweepTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dP_dpresweepTip_fd, dP_dpresweepTip, rtol=1e-4, atol=1e-8)



    def test_dpresweepTip3(self):

        presweep = np.linspace(1, 10, self.n)
        presweepTip = 10.1
        precone = 0.0
        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=True, presweep=presweep, presweepTip=presweepTip)

        CP, CT, CQ, dCP_ds, dCT_ds, dCQ_ds, dCP_dv, dCT_dv, dCQ_dv = \
            rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dpresweepTip = dCT_ds[0, 6]
        dCQ_dpresweepTip = dCQ_ds[0, 6]
        dCP_dpresweepTip = dCP_ds[0, 6]

        pst = float(presweepTip)
        delta = 1e-6*pst
        pst += delta

        rotor = CCBlade(self.r, self.chord, self.theta, self.af, self.Rhub, self.Rtip,
            self.B, self.rho, self.mu, precone, self.tilt, self.yaw, self.shearExp,
            self.hubHt, self.nSector, derivatives=False, presweep=presweep, presweepTip=pst)

        CPd, CTd, CQd = rotor.evaluate([self.Uinf], [self.Omega], [self.pitch], coefficient=True)

        dCT_dpresweepTip_fd = (CTd - CT) / delta
        dCQ_dpresweepTip_fd = (CQd - CQ) / delta
        dCP_dpresweepTip_fd = (CPd - CP) / delta

        np.testing.assert_allclose(dCT_dpresweepTip_fd, dCT_dpresweepTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dCQ_dpresweepTip_fd, dCQ_dpresweepTip, rtol=1e-4, atol=1e-8)
        np.testing.assert_allclose(dCP_dpresweepTip_fd, dCP_dpresweepTip, rtol=1e-4, atol=1e-8)




if __name__ == '__main__':
    unittest.main()
