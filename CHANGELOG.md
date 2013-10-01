# CCBlade Changelog

## 1.0 (Sep 26, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[NEW]:

- analytic gradients!
- presweep is added, but only from a blade element perspective---there are no sweep corrections.
- added unit tests for gradients as compared to finite differencing

[CHANGE]:

- redefined the way precurvature is specified so that it was not implicitly coupled with rotor hub precone
- no longer add the 0.0 loads at root and tip
- simplified some of the outputs
- removed direct dependence on coordinate system methods so that it could more easily stand-alone


## 0.3.0 (July 24, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[NEW]:

- blade precurve capability.  This is specified by supplying an array in ``precone`` input rather than a float for just a hub precone.

[CHANGE]:

- consistent imports for using with NREL_WISDEM or as a stand-alone tool


## 0.2.0 (June 12, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[NEW]:

- added precone, tilt, yaw capabilities

- added wind shear option (power-law)

- packaged coordinate system definitions with CCBlade

- added calculations at a particular azimuth (and better azimuthal-averaged calculations)

[CHANGE]:

- better limit handling for these cases which now allow for reversed tangential flow

- improved NREL 5MW test case


## 0.1.0 (May 9, 2013)

Andrew Ning <andrew.ning@nrel.gov>

- initial relase