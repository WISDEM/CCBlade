subroutine inductionFactors(r, c, Rhub, Rtip, phi, cl, cd, B, &
    Vx, Vy, useCd, hubLoss, tipLoss, wakerotation, &
    fzero, a, ap)

    implicit none

    integer, parameter :: ReKi = selected_real_kind(15, 307)

    ! in
    real(ReKi), intent(in) :: r, c, Rhub, Rtip, phi, cl, cd
    integer, intent(in) :: B
    real(ReKi), intent(in) :: Vx, Vy
    logical, intent(in) :: useCd, hubLoss, tipLoss, wakerotation
    !f2py logical, optional, intent(in) :: useCd = 1, hubLoss = 1, tipLoss = 1, wakerotation = 1

    ! out
    real(ReKi), intent(out) :: fzero, a, ap

    ! local
    real(ReKi) :: pi, sigma_p, sphi, cphi, lambda_r
    real(ReKi) :: factortip, Ftip, factorhub, Fhub
    real(ReKi) :: k, kp, cn, ct, F
    real(ReKi) :: g1, g2, g3


    ! constants
    pi = 3.1415926535897932
    sigma_p = B/2.0/pi*c/r
    sphi = sin(phi)
    cphi = cos(phi)

    ! resolve into normal and tangential forces
    if ( .not. useCd ) then
        cn = cl*cphi
        ct = cl*sphi
    else
        cn = cl*cphi + cd*sphi
        ct = cl*sphi - cd*cphi
    end if

    ! Prandtl's tip and hub loss factor
    Ftip = 1.0
    if ( tipLoss ) then
        factortip = B/2.0*(Rtip - r)/(r*abs(sphi))
        Ftip = 2.0/pi*acos(exp(-factortip))
    end if

    Fhub = 1.0
    if ( hubLoss ) then
        factorhub = B/2.0*(r - Rhub)/(Rhub*abs(sphi))
        Fhub = 2.0/pi*acos(exp(-factorhub))
    end if

    F = Ftip * Fhub

    ! bem parameters
    k = sigma_p*cn/4/F/sphi/sphi
    kp = sigma_p*ct/4/F/sphi/cphi

    ! compute axial induction factor
    if (phi > 0) then  ! momentum/empirical

        ! update axial induction factor
        if (k <= 2.0/3) then  ! momentum state
            a = k/(1+k)

        else  ! Glauert(Buhl) correction
            if ( abs(k - (25.0/18/F - 1.0)) < 1e-6) then
                k = k + 1.0e-5  ! avoid singularity
            end if

            g1 = 2*F*k - (10.0/9-F)
            g2 = 2*F*k - (4.0/3-F)*F
            g3 = 2*F*k - (25.0/9-2*F)

            a = (g1 - sqrt(g2)) / g3

        end if

    else  ! propeller brake region (a and ap not directly used but update anyway)

        if (k > 1.0) then
            a = k/(k-1)
        else
            a = 0.0  ! dummy value
        end if

    end if

    ! compute tangential induction factor
    ap = kp/(1-kp)

    if (.not. wakerotation) then
        ap = 0.0
        kp = 0.0
    end if

    ! error function
    lambda_r = Vy/Vx
    if (phi > 0) then  ! momentum/empirical
        fzero = sphi/(1-a) - cphi/lambda_r*(1-kp)
    else  ! propeller brake region
        fzero = sphi*(1-k) - cphi/lambda_r*(1-kp)
    end if

end subroutine inductionFactors




subroutine relativeWind(phi, a, ap, Vx, Vy, pitch, &
    chord, theta, rho, mu, alpha, W, Re)

    implicit none

    integer, parameter :: ReKi = selected_real_kind(15, 307)

    ! in
    real(ReKi), intent(in) :: phi, a, ap, Vx, Vy, pitch
    real(ReKi), intent(in) :: chord, theta, rho, mu

    ! out
    real(ReKi), intent(out) :: alpha, W, Re

    ! angle of attack
    alpha = phi - (theta + pitch)

    ! avoid numerical errors when angle is close to 0 or 90 deg
    ! and other induction factor is at some ridiculous value
    ! this only occurs when iterating on Reynolds number
    ! during the phi sweep where a solution has not been found yet
    if ( abs(a) > 10 ) then
        W = Vy*(1+ap)/cos(phi)
    else if ( abs(ap) > 10 ) then
        W = Vx*(1-a)/sin(phi)
    else
        W = sqrt((Vx*(1-a))**2 + (Vy*(1+ap))**2)
    end if

    Re = rho * W * chord / mu


end subroutine relativeWind



subroutine defineCurvature(n, r, precurve, presweep, precone, x_az, y_az, z_az, cone)

    implicit none

    integer, parameter :: ReKi = selected_real_kind(15, 307)

    ! in
    integer, intent(in) :: n
    real(ReKi), dimension(n), intent(in) :: r, precurve, presweep
    real(ReKi), intent(in) :: precone

    ! out
    real(ReKi), dimension(n), intent(out) :: x_az, y_az, z_az, cone


    ! coordinate in azimuthal coordinate system
    ! az_coords = DirectionVector(precurve, presweep, r).bladeToAzimuth(precone)

    x_az = -r*sin(precone) + precurve*cos(precone)
    z_az = r*cos(precone) + precurve*sin(precone)
    y_az = presweep


    ! compute total coning angle for purposes of relative velocity
    cone(1) = atan2(-(x_az(2) - x_az(1)), z_az(2) - z_az(1))
    cone(2:n-1) = 0.5*(atan2(-(x_az(2:n-1) - x_az(1:n-2)), z_az(2:n-1) - z_az(1:n-2)) &
                       + atan2(-(x_az(3:n) - x_az(2:n-1)), z_az(3:n) - z_az(2:n-1)))
    cone(n) = atan2(-(x_az(n) - x_az(n-1)), z_az(n) - z_az(n-1))


end subroutine defineCurvature




subroutine windComponents(n, r, precurve, presweep, precone, yaw, tilt, azimuth, &
    Uinf, OmegaRPM, hubHt, shearExp, Vx, Vy)

    implicit none

    integer, parameter :: ReKi = selected_real_kind(15, 307)

    ! in
    integer, intent(in) :: n
    real(ReKi), dimension(n), intent(in) :: r, precurve, presweep
    real(ReKi), intent(in) :: precone, yaw, tilt, azimuth, Uinf, OmegaRPM, hubHt, shearExp

    ! out
    real(ReKi), dimension(n), intent(out) :: Vx, Vy

    ! local
    real(ReKi) :: sy, cy, st, ct, sa, ca, pi, Omega
    real(ReKi), dimension(n) :: cone, sc, cc, x_az, y_az, z_az
    real(ReKi), dimension(n) :: heightFromHub, V, Vwind_x, Vwind_y, Vrot_x, Vrot_y


    ! rename
    sy = sin(yaw)
    cy = cos(yaw)
    st = sin(tilt)
    ct = cos(tilt)
    sa = sin(azimuth)
    ca = cos(azimuth)
    pi = 3.1415926535897932
    Omega = OmegaRPM * pi/30.0


    call defineCurvature(n, r, precurve, presweep, precone, x_az, y_az, z_az, cone)
    sc = sin(cone)
    cc = cos(cone)


    ! get section heights in wind-aligned coordinate system
    ! heightFromHub = az_coords.azimuthToHub(azimuth).hubToYaw(tilt).z

    heightFromHub = (y_az*sa + z_az*ca)*ct - x_az*st


    ! velocity with shear

    V = Uinf*(1 + heightFromHub/hubHt)**shearExp


    ! transform wind to blade c.s.
    ! Vwind = DirectionVector(V, 0*V, 0*V).windToYaw(yaw).yawToHub(tilt).hubToAzimuth(azimuth).azimuthToBlade(cone)

    Vwind_x = V * ((cy*st*ca + sy*sa)*sc + cy*ct*cc)
    Vwind_y = V * (cy*st*sa - sy*ca)


    ! wind from rotation to blade c.s.
    ! OmegaV = DirectionVector(Omega, 0.0, 0.0)
    ! Vrot = -OmegaV.cross(az_coords)  # negative sign because relative wind opposite to rotation
    ! Vrot = Vrot.azimuthToBlade(cone)

    Vrot_x = -Omega*y_az*sc
    Vrot_y = Omega*z_az


    ! total velocity
    Vx = Vwind_x + Vrot_x
    Vy = Vwind_y + Vrot_y



end subroutine windComponents







subroutine thrustTorque(n, Np, Tp, r, precurve, precone, Rhub, Rtip, T, Q)

    implicit none

    integer, parameter :: ReKi = selected_real_kind(15, 307)

    ! in
    integer, intent(in) :: n
    real(ReKi), dimension(n), intent(in) :: Np, Tp, r, precurve
    real(ReKi), intent(in) :: precone, Rhub, Rtip

    ! out
    real(ReKi), intent(out) :: T, Q

    ! local
    real(ReKi) :: dr
    real(ReKi), dimension(n) :: x_az, y_az, z_az, cone, presweep
    real(ReKi), dimension(n+2) :: rfull, thrust, torque
    integer :: i



    ! get z_az and total cone angle
    presweep = 0.0  ! the value is irrelevant for z_az and cone
    call defineCurvature(n, r, precurve, presweep, precone, x_az, y_az, z_az, cone)

    ! compute sectional contribution of thrust and torque
    rfull(2:n+1) = r
    thrust(2:n+1) = Np*cos(cone)
    torque(2:n+1) = Tp*z_az

    ! add hub/tip for complete integration.  loads go to zero at hub/tip.
    rfull(1) = Rhub
    rfull(n+2) = Rtip
    thrust(1) = 0.0
    thrust(n+2) = 0.0
    torque(1) = 0.0
    torque(n+2) = 0.0

    ! integrate Thrust and Torque (trapezoidal)
    T = 0.0
    do i = 1, n+1
        dr = rfull(i+1) - rfull(i)
        T = T + 0.5*(thrust(i) + thrust(i+1))*dr
        Q = Q + 0.5*(torque(i) + torque(i+1))*dr
    end do


end subroutine thrustTorque