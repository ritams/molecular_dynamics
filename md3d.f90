program md
    implicit none
    real(kind = 8), allocatable :: x(:), y(:), z(:)
    real(kind = 8), allocatable :: vx(:), vy(:), vz(:)
    real(kind = 8), allocatable :: fx(:), fy(:), fz(:), oldfx(:), oldfy(:), oldfz(:)
    real(kind = 8) :: dim, rho, sigma, temp, T, eps, rc, rc2, dt, dt2, ke, pe, etotal, avT, avpe
    real*8 :: avmx, avmy, avmz, fc, vc 
    integer(kind = 8) :: n, nc, dof, i, itr

    n = 108  ! number of particle
    nc = 3
    dof = 3 * n
    rho = 1.0d0
    dim = (dfloat(n) / rho) ** (1.0d0 / 3.0d0)   ! box length

    sigma = 1.0d0
    eps = 1.0d0
    temp = 1.0d0 

    rc = 2.50d0  ! cutt off radius
    rc2 = rc * rc
    fc = 4.0d0 * (12.0d0 / rc ** 13 - 6.0d0 / rc ** 7)
    vc = 4.0d0 * (1.0d0 / rc ** 12 - 1.0d0 / rc ** 6) + fc * rc

    allocate(x(n), y(n), z(n))
    allocate(vx(n), vy(n), vz(n))
    allocate(fx(n), fy(n), fz(n), oldfx(n), oldfy(n), oldfz(n))

    call ini_position(x, y, z)
    call ini_velocity(vx, vy, vz, ke)
    call ini_force(fx, fy, fz, x, y, z, pe)

    avT = 0.0d0
    avpe = 0.0d0
    dt = 0.005d0
    dt2 = dt * dt 
    itr = 50000

    open(unit = 1, file = "energy.txt")
    open(unit = 2, file = "Temperature.txt")
    open(unit = 3, file = "average_momentum.txt")

    do i = 1, itr
        !if (modulo(i, 1000) == 0) write(*, *) i
        call integrate(x, y, z, vx, vy, vz, fx, fy, fz, ke, pe, etotal, T)
        call momentum(avmx, avmy, avmz)


        ! if (i > 50000) then 
        !     avT = avT + T 
        !     avpe = avpe + pe  
        ! end if 
        
        write(3, *) avmx, avmy, avmz
        write(2, *) T
        write(1, *) ke, pe, etotal, i * dt  
    end do

    avT = avT / (itr / 2.0d0)
    avpe = avpe / (itr / 2.0d0)
    write(*, *) "average temperature", avT, avpe 









    contains

    subroutine ini_position(x, y, z)
            implicit none
            real(kind = 8), allocatable :: x(:), y(:), z(:)
            real(kind = 8) :: cell, cell2
            integer(kind = 8) :: i, j, k, m, ip  
        
            cell = 1.0d0 / dfloat(nc)
            cell2  = 0.50d0 * cell
        
            ! bulding the unit cell
        
            x(1) = 0.0d0
            y(1) = 0.0d0
            z(1) = 0.0d0 
        
            x(2) = cell2
            y(2) = cell2
            z(2) = 0.0d0
        
            x(3) = 0.0d0
            y(3) = cell2
            z(3) = cell2
        
            x(4) = cell2
            y(4) = 0.0d0
            z(4) = cell2
        
            m = 0
        
            do i = 1, nc
                do j = 1, nc
                    do k = 1, nc
                        do ip = 1, 4
                            x(ip + m) = x(ip) + (i - 1) * cell
                            y(ip + m) = y(ip) + (j - 1) * cell
                            z(ip + m) = z(ip) + (k - 1) * cell
                        end do
        
                        m = m + 4
                    end do
                end do
            end do
        
            do i = 1, n
                x(i) = x(i) * dim
                y(i) = y(i) * dim
                z(i) = z(i) * dim
            end do
        
    end subroutine

    subroutine ini_velocity(vx, vy, vz, ke)
        implicit none
        real(kind = 8), allocatable, intent(inout) :: vx(:), vy(:), vz(:)
        real(kind = 8), intent(out) :: ke
        real(kind = 8) :: v2, scale, vxav, vyav, vzav

        do i = 1, n 
            vx(i) = (rand() - 0.50d0)
            vy(i) = (rand() - 0.50d0)
            vz(i) = (rand() - 0.50d0)
        end do

        v2 = sum(vx * vx + vy * vy + vz * vz)
        v2 = v2 / dfloat(n)
        scale = dsqrt(3 * temp / v2)

        vxav = sum(vx) / dfloat(n)
        vyav = sum(vy) / dfloat(n)
        vzav = sum(vz) / dfloat(n)

        vx = (vx - vxav) * scale
        vy = (vy - vyav) * scale
        vz = (vz - vzav) * scale

        call momentum(avmx, avmy, avmz)
        ke = 0.50d0 * sum(vx * vx + vy * vy  + vz * vz)
        ke = ke / dfloat(n)

    end subroutine

    subroutine momentum(avmx, avmy, avmz)
        implicit none
        real(kind = 8), intent(inout) :: avmx, avmy, avmz
        
        avmx = sum(vx) / dfloat(n)
        avmy = sum(vy) / dfloat(n)
        avmz = sum(vz) / dfloat(n)

    end subroutine momentum


    subroutine ini_force(fx, fy, fz, x, y, z, pe)
        implicit none
        real(kind = 8), allocatable, intent(inout) :: fx(:), fy(:), fz(:), x(:), y(:), z(:)
        real(kind = 8), intent(inout) :: pe
        real(kind = 8) :: dx, dy, dz, dr2, dr, r2, r6, fr
        integer(kind = 8) :: i, j

        fx = 0.0d0
        fy = 0.0d0 
        fz = 0.0d0
        pe = 0.0d0

        do i = 1, n-1
            do j = i + 1, n 

                dx = x(i) - x(j)
                dx = dx - dim * anint(dx / dim)

                dy = y(i) - y(j)
                dy = dy - dim * anint(dy / dim)
 
                dz = z(i) - z(j)
                dz = dz - dim * anint(dz / dim)

                dr2 = dx * dx + dy * dy + dz * dz
                dr = dsqrt(dr2)

                if(dr2 < rc2) then

                    r2 = 1.0d0 / dr2
                    r6 = r2 ** 3
                    fr = 48.0d0 * dr2 * r6 * (r6 - 0.50d0) - fc
                    pe = pe + 4.0d0 * r6 * (r6 - 1.0d0) + fc * dr - vc

                    fx(i) = fx(i) + fr * dx 
                    fy(i) = fy(i) + fr * dy  
                    fz(i) = fz(i) + fr * dz 

                    fx(j) = fx(j) - fr * dx 
                    fy(j) = fy(j) - fr * dy
                    fz(j) = fz(j) - fr * dz  

                end if 
            end do
        end do
        pe = pe / dfloat(n)
    end subroutine

    subroutine integrate(x, y, z, vx, vy, vz, fx, fy, fz, ke, pe, etotal, T)
        implicit none
        real*8, allocatable, intent(inout) :: x(:), y(:), z(:)
        real*8, allocatable, intent(inout) :: vx(:), vy(:), vz(:)
        real*8, allocatable, intent(inout) :: fx(:), fy(:), fz(:)
        real*8, intent(inout) :: pe, ke, etotal, T

        x = modulo(x + vx * dt + 0.5d0 * fx * dt2, dim)
        y = modulo(y + vy * dt + 0.5d0 * fy * dt2, dim)       ! update in positions
        z = modulo(z + vz * dt + 0.5d0 * fz * dt2, dim)

        oldfx = fx
        oldfy = fy
        oldfz = fz
        call ini_force(fx, fy, fz, x, y, z, pe)

        vx = vx + dt * 0.5d0 * (oldfx + fx)
        vy = vy + dt * 0.5d0 * (oldfy + fy)     ! update in velocities
        vz = vz + dt * 0.5d0 * (oldfz + fz)
        
        ke = 0.5d0 * sum(vx * vx + vy * vy + vz * vz)
        ke = ke / dfloat(n)
        T = (2.0d0 / 3.0d0) * ke

        etotal = ke + pe 
        !call momentum(avmx, avmy, avmz)
    
    
    end subroutine

end program md