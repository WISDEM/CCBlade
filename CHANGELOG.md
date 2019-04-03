# CCBlade Changelog

## 1.2.0 (Dec 18, 2019)

Pietro Bortolotti<pietro.bortolotti@nrel.gov>

- Alternative solution formulation using Cl and Cd


## 1.1.1 (Apr 15, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[FIX]: improvement to the behavior of the function and its derivative near a singularity that can exist in the Glauert region.


## 1.1.0 (Jan 10, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[NEW]:

- added the remaining analytic gradients I didn't think I would need (yaw, azimuth, Uinf, Omega, pitch)
- added analytic gradients for parked cases (Omega==0)
- additional unit tests for the new gradients

[CHANGE]:

- gradients are now returned as a dictionary, and full Jacobians are returned.  This makes it much easier to access the gradients, because you don't have to worry about what the ordering is.  Also returning 2D arrays for every value maintains consistency so you don't have to remember which gradients were really just the nonzero diagonal entries of the corresponding Jacobian.  However, this makes the API not backwards compatible.  I apologize if it causes a lot of changes on your end.
- AirfoilPrep.py was moved to a separate repository (https://github.com/NREL-WISDEM/AirfoilPreppy) to more easily accommodate changes and contributions.  The setup.py script should still install from the repository automatically.  For some reason the setup script throws an error saying it could not find AirfoilPrep however it does find it and installs properly.

[FIX]:

- As compared to what I wrote in the docstring, I was returning the transpose of dNp_dprecurve, dTp_dprecurve, dNp_dpresweep, and dTp_dpresweep.  This is now fixed.

## 1.0.2 (Nov 21, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[FIX]:

- same fix as below, but for the dimensional versions (P, T, Q)


## 1.0.1 (Nov 8, 2013)

Andrew Ning <andrew.ning@nrel.gov>

[FIX]:

- derivatives of force coeffcients (CP, CT, CQ) were not handled correctly if array inputs were used

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
