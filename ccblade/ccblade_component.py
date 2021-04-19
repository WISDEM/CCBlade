import numpy as np
from openmdao.api import ExplicitComponent
from ccblade.ccblade import CCBlade, CCAirfoil

cosd = lambda x: np.cos(np.deg2rad(x))
sind = lambda x: np.sin(np.deg2rad(x))


class CCBladeGeometry(ExplicitComponent):
    """
    Compute some geometric properties of the turbine based on the tip radius,
    precurve, presweep, and precone.

    Parameters
    ----------
    Rtip : float
        Rotor tip radius.
    precurve_in : numpy array[n_span]
        Prebend distribution along the span.
    presweep_in : numpy array[n_span]
        Presweep distribution along the span.
    precone : float
        Precone angle.

    Returns
    -------
    R : float
        Rotor radius.
    diameter : float
        Rotor diameter.
    precurveTip : float
        Precurve value at the rotor tip.
    presweepTip : float
        Presweep value at the rotor tip.
    """

    def initialize(self):
        self.options.declare("n_span")

    def setup(self):
        n_span = self.options["n_span"]

        self.add_input("Rtip", val=0.0, units="m")
        self.add_input("precurve_in", val=np.zeros(n_span), units="m")
        self.add_input("presweep_in", val=np.zeros(n_span), units="m")
        self.add_input("precone", val=0.0, units="deg")

        self.add_output("R", val=0.0, units="m")
        self.add_output("diameter", val=0.0, units="m")
        self.add_output("precurveTip", val=0.0, units="m")
        self.add_output("presweepTip", val=0.0, units="m")

        self.declare_partials("R", ["Rtip", "precone"])
        self.declare_partials("diameter", ["Rtip", "precone"])

        self.declare_partials(["R", "diameter"], "precurve_in", rows=[0], cols=[n_span - 1])

        self.declare_partials("precurveTip", "precurve_in", val=1.0, rows=[0], cols=[n_span - 1])
        self.declare_partials("presweepTip", "presweep_in", val=1.0, rows=[0], cols=[n_span - 1])

    def compute(self, inputs, outputs):
        Rtip = inputs["Rtip"]
        precone = inputs["precone"]

        outputs["precurveTip"] = inputs["precurve_in"][-1]
        outputs["presweepTip"] = inputs["presweep_in"][-1]

        outputs["R"] = Rtip * cosd(precone) + outputs["precurveTip"] * sind(precone)
        outputs["diameter"] = outputs["R"] * 2

    def compute_partials(self, inputs, J):
        Rtip = inputs["Rtip"]
        precone = inputs["precone"]
        precurveTip = inputs["precurve_in"][-1]

        J["R", "precurve_in"] = sind(precone)
        J["R", "Rtip"] = cosd(precone)
        J["R", "precone"] = (-Rtip * sind(precone) + precurveTip * cosd(precone)) * np.pi / 180.0

        J["diameter", "precurve_in"] = 2.0 * J["R", "precurve_in"]
        J["diameter", "Rtip"] = 2.0 * J["R", "Rtip"]
        J["diameter", "precone"] = 2.0 * J["R", "precone"]


class CCBladeLoads(ExplicitComponent):
    """
    Compute the aerodynamic forces along the blade span given a rotor speed,
    pitch angle, and wind speed.

    This component instantiates and calls a CCBlade instance to compute the loads.
    Analytic derivatives are provided for all inptus except all airfoils*,
    mu, rho, and shearExp.

    Parameters
    ----------
    V_load : float
        Hub height wind speed.
    Omega_load : float
        Rotor rotation speed.
    pitch_load : float
        Blade pitch setting.
    azimuth_load : float
        Blade azimuthal location.
    r : numpy array[n_span]
        Radial locations where blade is defined. Should be increasing and not
        go all the way to hub or tip.
    chord : numpy array[n_span]
        Chord length at each section.
    theta : numpy array[n_span]
        Twist angle at each section (positive decreases angle of attack).
    Rhub : float
        Hub radius.
    Rtip : float
        Tip radius.
    hub_height : float
        Hub height.
    precone : float
        Precone angle.
    tilt : float
        Shaft tilt.
    yaw : float
        Yaw error.
    precurve : numpy array[n_span]
        Precurve at each section.
    precurveTip : float
        Precurve at tip.
    airfoils_cl : numpy array[n_span, n_aoa, n_Re, n_tab]
        Lift coefficients, spanwise.
    airfoils_cd : numpy array[n_span, n_aoa, n_Re, n_tab]
        Drag coefficients, spanwise.
    airfoils_cm : numpy array[n_span, n_aoa, n_Re, n_tab]
        Moment coefficients, spanwise.
    airfoils_aoa : numpy array[n_aoa]
        Angle of attack grid for polars.
    airfoils_Re : numpy array[n_Re]
        Reynolds numbers of polars.
    nBlades : int
        Number of blades
    rho : float
        Density of air
    mu : float
        Dynamic viscosity of air
    shearExp : float
        Shear exponent.
    nSector : int
        Number of sectors to divide rotor face into in computing thrust and power.
    tiploss : boolean
        Include Prandtl tip loss model.
    hubloss : boolean
        Include Prandtl hub loss model.
    wakerotation : boolean
        Include effect of wake rotation (i.e., tangential induction factor is nonzero).
    usecd : boolean
        Use drag coefficient in computing induction factors.

    Returns
    -------
    loads_r : numpy array[n_span]
        Radial positions along blade going toward tip.
    loads_Px : numpy array[n_span]
         Distributed loads in blade-aligned x-direction.
    loads_Py : numpy array[n_span]
         Distributed loads in blade-aligned y-direction.
    loads_Pz : numpy array[n_span]
         Distributed loads in blade-aligned z-direction.
    """

    def initialize(self):
        self.options.declare("modeling_options")

    def setup(self):
        rotorse_options = self.options["modeling_options"]["WISDEM"]["RotorSE"]
        self.n_span = n_span = rotorse_options["n_span"]
        self.n_aoa = n_aoa = rotorse_options["n_aoa"]  # Number of angle of attacks
        self.n_Re = n_Re = rotorse_options["n_Re"]  # Number of Reynolds
        self.n_tab = n_tab = rotorse_options[
            "n_tab"
        ]  # Number of tabulated data. For distributed aerodynamic control this could be > 1

        # inputs
        self.add_input("V_load", val=20.0, units="m/s")
        self.add_input("Omega_load", val=0.0, units="rpm")
        self.add_input("pitch_load", val=0.0, units="deg")
        self.add_input("azimuth_load", val=0.0, units="deg")

        self.add_input("r", val=np.zeros(n_span), units="m")
        self.add_input("chord", val=np.zeros(n_span), units="m")
        self.add_input("theta", val=np.zeros(n_span), units="deg")
        self.add_input("Rhub", val=0.0, units="m")
        self.add_input("Rtip", val=0.0, units="m")
        self.add_input("hub_height", val=0.0, units="m")
        self.add_input("precone", val=0.0, units="deg")
        self.add_input("tilt", val=0.0, units="deg")
        self.add_input("yaw", val=0.0, units="deg")
        self.add_input("precurve", val=np.zeros(n_span), units="m")
        self.add_input("precurveTip", val=0.0, units="m")

        # parameters
        self.add_input("airfoils_cl", val=np.zeros((n_span, n_aoa, n_Re, n_tab)))
        self.add_input("airfoils_cd", val=np.zeros((n_span, n_aoa, n_Re, n_tab)))
        self.add_input("airfoils_cm", val=np.zeros((n_span, n_aoa, n_Re, n_tab)))
        self.add_input("airfoils_aoa", val=np.zeros((n_aoa)), units="deg")
        self.add_input("airfoils_Re", val=np.zeros((n_Re)))

        self.add_discrete_input("nBlades", val=0)
        self.add_input("rho", val=0.0, units="kg/m**3")
        self.add_input("mu", val=0.0, units="kg/(m*s)")
        self.add_input("shearExp", val=0.0)
        self.add_discrete_input("nSector", val=4)
        self.add_discrete_input("tiploss", val=True)
        self.add_discrete_input("hubloss", val=True)
        self.add_discrete_input("wakerotation", val=True)
        self.add_discrete_input("usecd", val=True)

        # outputs
        self.add_output("loads_r", val=np.zeros(n_span), units="m")
        self.add_output("loads_Px", val=np.zeros(n_span), units="N/m")
        self.add_output("loads_Py", val=np.zeros(n_span), units="N/m")
        self.add_output("loads_Pz", val=np.zeros(n_span), units="N/m")

        arange = np.arange(n_span)
        self.declare_partials(
            "loads_Px",
            [
                "Omega_load",
                "Rhub",
                "Rtip",
                "V_load",
                "azimuth_load",
                "chord",
                "hub_height",
                "pitch_load",
                "precone",
                "precurve",
                "r",
                "theta",
                "tilt",
                "yaw",
                "shearExp",
            ],
        )
        self.declare_partials(
            "loads_Py",
            [
                "Omega_load",
                "Rhub",
                "Rtip",
                "V_load",
                "azimuth_load",
                "chord",
                "hub_height",
                "pitch_load",
                "precone",
                "precurve",
                "r",
                "theta",
                "tilt",
                "yaw",
                "shearExp",
            ],
        )
        self.declare_partials("loads_Pz", "*", dependent=False)
        self.declare_partials("loads_r", "r", val=1.0, rows=arange, cols=arange)
        self.declare_partials("*", "airfoils*", dependent=False)

    def compute(self, inputs, outputs, discrete_inputs, discrete_outputs):
        r = inputs["r"]
        chord = inputs["chord"]
        theta = inputs["theta"]
        Rhub = inputs["Rhub"]
        Rtip = inputs["Rtip"]
        hub_height = inputs["hub_height"]
        precone = inputs["precone"]
        tilt = inputs["tilt"]
        yaw = inputs["yaw"]
        precurve = inputs["precurve"]
        precurveTip = inputs["precurveTip"]
        B = discrete_inputs["nBlades"]
        rho = inputs["rho"]
        mu = inputs["mu"]
        shearExp = inputs["shearExp"]
        nSector = discrete_inputs["nSector"]
        tiploss = discrete_inputs["tiploss"]
        hubloss = discrete_inputs["hubloss"]
        wakerotation = discrete_inputs["wakerotation"]
        usecd = discrete_inputs["usecd"]
        V_load = inputs["V_load"]
        Omega_load = inputs["Omega_load"]
        pitch_load = inputs["pitch_load"]
        azimuth_load = inputs["azimuth_load"]

        if len(precurve) == 0:
            precurve = np.zeros_like(r)

        # airfoil files
        af = [None] * self.n_span
        for i in range(self.n_span):
            af[i] = CCAirfoil(
                inputs["airfoils_aoa"],
                inputs["airfoils_Re"],
                inputs["airfoils_cl"][i, :, :, 0],
                inputs["airfoils_cd"][i, :, :, 0],
                inputs["airfoils_cm"][i, :, :, 0],
            )

        ccblade = CCBlade(
            r,
            chord,
            theta,
            af,
            Rhub,
            Rtip,
            B,
            rho,
            mu,
            precone,
            tilt,
            yaw,
            shearExp,
            hub_height,
            nSector,
            precurve,
            precurveTip,
            tiploss=tiploss,
            hubloss=hubloss,
            wakerotation=wakerotation,
            usecd=usecd,
            derivatives=True,
        )

        # distributed loads
        loads, self.derivs = ccblade.distributedAeroLoads(V_load, Omega_load, pitch_load, azimuth_load)
        Np = loads["Np"]
        Tp = loads["Tp"]

        # unclear why we need this output at all
        outputs["loads_r"] = r

        # conform to blade-aligned coordinate system
        outputs["loads_Px"] = Np
        outputs["loads_Py"] = -Tp
        outputs["loads_Pz"][:] = 0.0

    def compute_partials(self, inputs, J, discrete_inputs):
        dNp = self.derivs["dNp"]
        dTp = self.derivs["dTp"]

        J["loads_Px", "r"] = dNp["dr"]
        J["loads_Px", "chord"] = dNp["dchord"]
        J["loads_Px", "theta"] = dNp["dtheta"]
        J["loads_Px", "Rhub"] = np.squeeze(dNp["dRhub"])
        J["loads_Px", "Rtip"] = np.squeeze(dNp["dRtip"])
        J["loads_Px", "hub_height"] = np.squeeze(dNp["dhubHt"])
        J["loads_Px", "precone"] = np.squeeze(dNp["dprecone"])
        J["loads_Px", "tilt"] = np.squeeze(dNp["dtilt"])
        J["loads_Px", "yaw"] = np.squeeze(dNp["dyaw"])
        J["loads_Px", "shearExp"] = np.squeeze(dNp["dshear"])
        J["loads_Px", "V_load"] = np.squeeze(dNp["dUinf"])
        J["loads_Px", "Omega_load"] = np.squeeze(dNp["dOmega"])
        J["loads_Px", "pitch_load"] = np.squeeze(dNp["dpitch"])
        J["loads_Px", "azimuth_load"] = np.squeeze(dNp["dazimuth"])
        J["loads_Px", "precurve"] = dNp["dprecurve"]

        J["loads_Py", "r"] = -dTp["dr"]
        J["loads_Py", "chord"] = -dTp["dchord"]
        J["loads_Py", "theta"] = -dTp["dtheta"]
        J["loads_Py", "Rhub"] = -np.squeeze(dTp["dRhub"])
        J["loads_Py", "Rtip"] = -np.squeeze(dTp["dRtip"])
        J["loads_Py", "hub_height"] = -np.squeeze(dTp["dhubHt"])
        J["loads_Py", "precone"] = -np.squeeze(dTp["dprecone"])
        J["loads_Py", "tilt"] = -np.squeeze(dTp["dtilt"])
        J["loads_Py", "yaw"] = -np.squeeze(dTp["dyaw"])
        J["loads_Py", "shearExp"] = -np.squeeze(dTp["dshear"])
        J["loads_Py", "V_load"] = -np.squeeze(dTp["dUinf"])
        J["loads_Py", "Omega_load"] = -np.squeeze(dTp["dOmega"])
        J["loads_Py", "pitch_load"] = -np.squeeze(dTp["dpitch"])
        J["loads_Py", "azimuth_load"] = -np.squeeze(dTp["dazimuth"])
        J["loads_Py", "precurve"] = -dTp["dprecurve"]



class CCBladeEvaluate(ExplicitComponent):
    """
    Standalone component for CCBlade that is only a light wrapper on CCBlade().

    Currently, this component is not used in any workflow, but it is a
    convenient way to test the derivatives coming out of CCBlade using OpenMDAO's
    check_partials method.

    """

    def initialize(self):
        self.options.declare("modeling_options")

    def setup(self):
        rotorse_init_options = self.options["modeling_options"]["WISDEM"]["RotorSE"]
        self.n_span = n_span = rotorse_init_options["n_span"]
        self.n_aoa = n_aoa = rotorse_init_options["n_aoa"]  # Number of angle of attacks
        self.n_Re = n_Re = rotorse_init_options["n_Re"]  # Number of Reynolds
        self.n_tab = n_tab = rotorse_init_options[
            "n_tab"
        ]  # Number of tabulated data. For distributed aerodynamic control this could be > 1

        # inputs
        self.add_input("V_load", val=20.0, units="m/s")
        self.add_input("Omega_load", val=9.0, units="rpm")
        self.add_input("pitch_load", val=0.0, units="deg")

        self.add_input("r", val=np.zeros(n_span), units="m")
        self.add_input("chord", val=np.zeros(n_span), units="m")
        self.add_input("theta", val=np.zeros(n_span), units="deg")
        self.add_input("Rhub", val=0.0, units="m")
        self.add_input("Rtip", val=0.0, units="m")
        self.add_input("hub_height", val=0.0, units="m")
        self.add_input("precone", val=0.0, units="deg")
        self.add_input("tilt", val=0.0, units="deg")
        self.add_input("yaw", val=0.0, units="deg")
        self.add_input("precurve", val=np.zeros(n_span), units="m")
        self.add_input("precurveTip", val=0.0, units="m")

        # parameters
        self.add_input("airfoils_cl", val=np.zeros((n_span, n_aoa, n_Re, n_tab)))
        self.add_input("airfoils_cd", val=np.zeros((n_span, n_aoa, n_Re, n_tab)))
        self.add_input("airfoils_cm", val=np.zeros((n_span, n_aoa, n_Re, n_tab)))
        self.add_input("airfoils_aoa", val=np.zeros((n_aoa)), units="deg")
        self.add_input("airfoils_Re", val=np.zeros((n_Re)))

        self.add_discrete_input("nBlades", val=0)
        self.add_input("rho", val=0.0, units="kg/m**3")
        self.add_input("mu", val=0.0, units="kg/(m*s)")
        self.add_input("shearExp", val=0.0)
        self.add_discrete_input("nSector", val=4)
        self.add_discrete_input("tiploss", val=True)
        self.add_discrete_input("hubloss", val=True)
        self.add_discrete_input("wakerotation", val=True)
        self.add_discrete_input("usecd", val=True)

        # outputs
        self.add_output("P", val=0.0, units="W")
        self.add_output("T", val=0.0, units="N")
        self.add_output("Q", val=0.0, units="N/m")
        self.add_output("M", val=0.0, units="N/m")

        self.add_output("CP", val=0.0)
        self.add_output("CT", val=0.0)
        self.add_output("CQ", val=0.0)
        self.add_output("CM", val=0.0)

        self.declare_partials("*", "*")
        self.declare_partials("*", "airfoils*", dependent=False)

    def compute(self, inputs, outputs, discrete_inputs, discrete_outputs):
        r = inputs["r"]
        chord = inputs["chord"]
        theta = inputs["theta"]
        Rhub = inputs["Rhub"]
        Rtip = inputs["Rtip"]
        hub_height = inputs["hub_height"]
        precone = inputs["precone"]
        tilt = inputs["tilt"]
        yaw = inputs["yaw"]
        precurve = inputs["precurve"]
        precurveTip = inputs["precurveTip"]
        B = discrete_inputs["nBlades"]
        rho = inputs["rho"]
        mu = inputs["mu"]
        shearExp = inputs["shearExp"]
        nSector = discrete_inputs["nSector"]
        tiploss = discrete_inputs["tiploss"]
        hubloss = discrete_inputs["hubloss"]
        wakerotation = discrete_inputs["wakerotation"]
        usecd = discrete_inputs["usecd"]
        V_load = inputs["V_load"]
        Omega_load = inputs["Omega_load"]
        pitch_load = inputs["pitch_load"]

        if len(precurve) == 0:
            precurve = np.zeros_like(r)

        # airfoil files
        af = [None] * self.n_span
        for i in range(self.n_span):
            af[i] = CCAirfoil(
                inputs["airfoils_aoa"],
                inputs["airfoils_Re"],
                inputs["airfoils_cl"][i, :, :, 0],
                inputs["airfoils_cd"][i, :, :, 0],
                inputs["airfoils_cm"][i, :, :, 0],
            )

        ccblade = CCBlade(
            r,
            chord,
            theta,
            af,
            Rhub,
            Rtip,
            B,
            rho,
            mu,
            precone,
            tilt,
            yaw,
            shearExp,
            hub_height,
            nSector,
            precurve,
            precurveTip,
            tiploss=tiploss,
            hubloss=hubloss,
            wakerotation=wakerotation,
            usecd=usecd,
            derivatives=False,
        )

        loads, _ = ccblade.evaluate(V_load, Omega_load, pitch_load)
        outputs["P"] = loads["P"]
        outputs["T"] = loads["T"]
        outputs["Q"] = loads["Q"]
        outputs["M"] = loads["M"]

        loads, _ = ccblade.evaluate(V_load, Omega_load, pitch_load, coefficients=True)
        outputs["CP"] = loads["CP"]
        outputs["CT"] = loads["CT"]
        outputs["CQ"] = loads["CQ"]
        outputs["CM"] = loads["CM"]

    def compute_partials(self, inputs, J, discrete_inputs):
        r = inputs["r"]
        chord = inputs["chord"]
        theta = inputs["theta"]
        Rhub = inputs["Rhub"]
        Rtip = inputs["Rtip"]
        hub_height = inputs["hub_height"]
        precone = inputs["precone"]
        tilt = inputs["tilt"]
        yaw = inputs["yaw"]
        precurve = inputs["precurve"]
        precurveTip = inputs["precurveTip"]
        B = discrete_inputs["nBlades"]
        rho = inputs["rho"]
        mu = inputs["mu"]
        shearExp = inputs["shearExp"]
        nSector = discrete_inputs["nSector"]
        tiploss = discrete_inputs["tiploss"]
        hubloss = discrete_inputs["hubloss"]
        wakerotation = discrete_inputs["wakerotation"]
        usecd = discrete_inputs["usecd"]
        V_load = inputs["V_load"]
        Omega_load = inputs["Omega_load"]
        pitch_load = inputs["pitch_load"]

        if len(precurve) == 0:
            precurve = np.zeros_like(r)

        # airfoil files
        af = [None] * self.n_span
        for i in range(self.n_span):
            af[i] = CCAirfoil(
                inputs["airfoils_aoa"],
                inputs["airfoils_Re"],
                inputs["airfoils_cl"][i, :, :, 0],
                inputs["airfoils_cd"][i, :, :, 0],
                inputs["airfoils_cm"][i, :, :, 0],
            )

        ccblade = CCBlade(
            r,
            chord,
            theta,
            af,
            Rhub,
            Rtip,
            B,
            rho,
            mu,
            precone,
            tilt,
            yaw,
            shearExp,
            hub_height,
            nSector,
            precurve,
            precurveTip,
            tiploss=tiploss,
            hubloss=hubloss,
            wakerotation=wakerotation,
            usecd=usecd,
            derivatives=True,
        )

        loads, derivs = ccblade.evaluate(V_load, Omega_load, pitch_load)

        dP = derivs["dP"]
        J["P", "r"] = dP["dr"]
        J["P", "chord"] = dP["dchord"]
        J["P", "theta"] = dP["dtheta"]
        J["P", "Rhub"] = np.squeeze(dP["dRhub"])
        J["P", "Rtip"] = np.squeeze(dP["dRtip"])
        J["P", "hub_height"] = np.squeeze(dP["dhubHt"])
        J["P", "precone"] = np.squeeze(dP["dprecone"])
        J["P", "tilt"] = np.squeeze(dP["dtilt"])
        J["P", "yaw"] = np.squeeze(dP["dyaw"])
        J["P", "shearExp"] = np.squeeze(dP["dshear"])
        J["P", "V_load"] = np.squeeze(dP["dUinf"])
        J["P", "Omega_load"] = np.squeeze(dP["dOmega"])
        J["P", "pitch_load"] = np.squeeze(dP["dpitch"])
        J["P", "precurve"] = dP["dprecurve"]
        J["P", "precurveTip"] = dP["dprecurveTip"]

        dT = derivs["dT"]
        J["T", "r"] = dT["dr"]
        J["T", "chord"] = dT["dchord"]
        J["T", "theta"] = dT["dtheta"]
        J["T", "Rhub"] = np.squeeze(dT["dRhub"])
        J["T", "Rtip"] = np.squeeze(dT["dRtip"])
        J["T", "hub_height"] = np.squeeze(dT["dhubHt"])
        J["T", "precone"] = np.squeeze(dT["dprecone"])
        J["T", "tilt"] = np.squeeze(dT["dtilt"])
        J["T", "yaw"] = np.squeeze(dT["dyaw"])
        J["T", "shearExp"] = np.squeeze(dT["dshear"])
        J["T", "V_load"] = np.squeeze(dT["dUinf"])
        J["T", "Omega_load"] = np.squeeze(dT["dOmega"])
        J["T", "pitch_load"] = np.squeeze(dT["dpitch"])
        J["T", "precurve"] = dT["dprecurve"]
        J["T", "precurveTip"] = dT["dprecurveTip"]

        dQ = derivs["dQ"]
        J["Q", "r"] = dQ["dr"]
        J["Q", "chord"] = dQ["dchord"]
        J["Q", "theta"] = dQ["dtheta"]
        J["Q", "Rhub"] = np.squeeze(dQ["dRhub"])
        J["Q", "Rtip"] = np.squeeze(dQ["dRtip"])
        J["Q", "hub_height"] = np.squeeze(dQ["dhubHt"])
        J["Q", "precone"] = np.squeeze(dQ["dprecone"])
        J["Q", "tilt"] = np.squeeze(dQ["dtilt"])
        J["Q", "yaw"] = np.squeeze(dQ["dyaw"])
        J["Q", "shearExp"] = np.squeeze(dQ["dshear"])
        J["Q", "V_load"] = np.squeeze(dQ["dUinf"])
        J["Q", "Omega_load"] = np.squeeze(dQ["dOmega"])
        J["Q", "pitch_load"] = np.squeeze(dQ["dpitch"])
        J["Q", "precurve"] = dQ["dprecurve"]
        J["Q", "precurveTip"] = dQ["dprecurveTip"]

        dM = derivs["dM"]
        J["M", "r"] = dM["dr"]
        J["M", "chord"] = dM["dchord"]
        J["M", "theta"] = dM["dtheta"]
        J["M", "Rhub"] = np.squeeze(dM["dRhub"])
        J["M", "Rtip"] = np.squeeze(dM["dRtip"])
        J["M", "hub_height"] = np.squeeze(dM["dhubHt"])
        J["M", "precone"] = np.squeeze(dM["dprecone"])
        J["M", "tilt"] = np.squeeze(dM["dtilt"])
        J["M", "yaw"] = np.squeeze(dM["dyaw"])
        J["M", "shearExp"] = np.squeeze(dM["dshear"])
        J["M", "V_load"] = np.squeeze(dM["dUinf"])
        J["M", "Omega_load"] = np.squeeze(dM["dOmega"])
        J["M", "pitch_load"] = np.squeeze(dM["dpitch"])
        J["M", "precurve"] = dM["dprecurve"]
        J["M", "precurveTip"] = dM["dprecurveTip"]

        loads, derivs = ccblade.evaluate(V_load, Omega_load, pitch_load, coefficients=True)

        dCP = derivs["dCP"]
        J["CP", "r"] = dCP["dr"]
        J["CP", "chord"] = dCP["dchord"]
        J["CP", "theta"] = dCP["dtheta"]
        J["CP", "Rhub"] = np.squeeze(dCP["dRhub"])
        J["CP", "Rtip"] = np.squeeze(dCP["dRtip"])
        J["CP", "hub_height"] = np.squeeze(dCP["dhubHt"])
        J["CP", "precone"] = np.squeeze(dCP["dprecone"])
        J["CP", "tilt"] = np.squeeze(dCP["dtilt"])
        J["CP", "yaw"] = np.squeeze(dCP["dyaw"])
        J["CP", "shearExp"] = np.squeeze(dCP["dshear"])
        J["CP", "V_load"] = np.squeeze(dCP["dUinf"])
        J["CP", "Omega_load"] = np.squeeze(dCP["dOmega"])
        J["CP", "pitch_load"] = np.squeeze(dCP["dpitch"])
        J["CP", "precurve"] = dCP["dprecurve"]
        J["CP", "precurveTip"] = dCP["dprecurveTip"]

        dCT = derivs["dCT"]
        J["CT", "r"] = dCT["dr"]
        J["CT", "chord"] = dCT["dchord"]
        J["CT", "theta"] = dCT["dtheta"]
        J["CT", "Rhub"] = np.squeeze(dCT["dRhub"])
        J["CT", "Rtip"] = np.squeeze(dCT["dRtip"])
        J["CT", "hub_height"] = np.squeeze(dCT["dhubHt"])
        J["CT", "precone"] = np.squeeze(dCT["dprecone"])
        J["CT", "tilt"] = np.squeeze(dCT["dtilt"])
        J["CT", "yaw"] = np.squeeze(dCT["dyaw"])
        J["CT", "shearExp"] = np.squeeze(dCT["dshear"])
        J["CT", "V_load"] = np.squeeze(dCT["dUinf"])
        J["CT", "Omega_load"] = np.squeeze(dCT["dOmega"])
        J["CT", "pitch_load"] = np.squeeze(dCT["dpitch"])
        J["CT", "precurve"] = dCT["dprecurve"]
        J["CT", "precurveTip"] = dCT["dprecurveTip"]

        dCQ = derivs["dCQ"]
        J["CQ", "r"] = dCQ["dr"]
        J["CQ", "chord"] = dCQ["dchord"]
        J["CQ", "theta"] = dCQ["dtheta"]
        J["CQ", "Rhub"] = np.squeeze(dCQ["dRhub"])
        J["CQ", "Rtip"] = np.squeeze(dCQ["dRtip"])
        J["CQ", "hub_height"] = np.squeeze(dCQ["dhubHt"])
        J["CQ", "precone"] = np.squeeze(dCQ["dprecone"])
        J["CQ", "tilt"] = np.squeeze(dCQ["dtilt"])
        J["CQ", "yaw"] = np.squeeze(dCQ["dyaw"])
        J["CQ", "shearExp"] = np.squeeze(dCQ["dshear"])
        J["CQ", "V_load"] = np.squeeze(dCQ["dUinf"])
        J["CQ", "Omega_load"] = np.squeeze(dCQ["dOmega"])
        J["CQ", "pitch_load"] = np.squeeze(dCQ["dpitch"])
        J["CQ", "precurve"] = dCQ["dprecurve"]
        J["CQ", "precurveTip"] = dCQ["dprecurveTip"]

        dCM = derivs["dCM"]
        J["CM", "r"] = dCM["dr"]
        J["CM", "chord"] = dCM["dchord"]
        J["CM", "theta"] = dCM["dtheta"]
        J["CM", "Rhub"] = np.squeeze(dCM["dRhub"])
        J["CM", "Rtip"] = np.squeeze(dCM["dRtip"])
        J["CM", "hub_height"] = np.squeeze(dCM["dhubHt"])
        J["CM", "precone"] = np.squeeze(dCM["dprecone"])
        J["CM", "tilt"] = np.squeeze(dCM["dtilt"])
        J["CM", "yaw"] = np.squeeze(dCM["dyaw"])
        J["CM", "shearExp"] = np.squeeze(dCM["dshear"])
        J["CM", "V_load"] = np.squeeze(dCM["dUinf"])
        J["CM", "Omega_load"] = np.squeeze(dCM["dOmega"])
        J["CM", "pitch_load"] = np.squeeze(dCM["dpitch"])
        J["CM", "precurve"] = dCM["dprecurve"]
        J["CM", "precurveTip"] = dCM["dprecurveTip"]
