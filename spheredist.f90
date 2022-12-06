    subroutine spheredist(stla0,stlo0,evla0,evlo0,dk)
    implicit none
    real*8, intent(in) :: stla0, stlo0, evla0, evlo0
    
    real*8, intent(out) :: dk
    
    real*8 :: stla, stlo, evla, evlo
    real*8 :: R, pi, pi1, d2r, h
    
    real*8 :: HaverSin
    
    R = 6371.00
    pi = 4.0*atan(1.0)
    pi1 = 4d0*atan(1d0)
    d2r = pi/180.0
    
    stla = stla0*d2r
    stlo = stlo0*d2r
    evla = evla0*d2r
    evlo = evlo0*d2r
    
    h = HaverSin(abs(stla-evla))+cos(stla)*cos(evla)*HaverSin(abs(stlo-evlo))
    dk = 2*R*asin(sqrt(h))
    
    return
    
    end subroutine spheredist
    
    function HaverSin(theta)
    real*8 :: theta
    real*8 :: HaverSin
    
    HaverSin = sin(theta/2.0)*sin(theta/2.0)
    
    return
    end function HaverSin
    