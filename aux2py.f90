subroutine aux2nc_liq(nbin,nx,ny,nz,net3d,rho_air,q3d_sbm,n3d_sbm,dm3d)
    use sdsu_io
    implicit none
    integer :: i, j, k, n, nbin, nx, ny, nz
    real,intent(in) :: net3d(nx,ny,nz,nbin), rho_air(nx,ny,nz)
    real,intent(out) :: q3d_sbm(nx,ny,nz)
    real,intent(out) :: n3d_sbm(nx,ny,nz,nbin)
    real,intent(out) :: dm3d(nx,ny,nz)
    n3d_sbm=0.e0 ; q3d_sbm=0.e0; dm3d=0.e0
    do n = 1, nbin
        do k = 1, nz ; do j = 1, ny ; do i = 1, nx
           !call sbm_function(x_sbm(n)%liq, drad_sbm(n)%liq*1.e-2, net3d(i,j,k,n), &
           !                  n3d_sbm(i,j,k,n), q3d_sbm(i,j,k))
           call sbm_function43_dm(rho_air(i,j,k), x_sbm(n)%liq, drad_sbm(n)%liq*1.e-2, rad_sbm(n)%liq, net3d(i,j,k,n), &
                             n3d_sbm(i,j,k,n), q3d_sbm(i,j,k), dm3d(i,j,k))
        enddo ; enddo ; enddo 
     enddo  
end subroutine aux2nc_liq


subroutine set_sdsu_param(nbin)
    use sdsu_io
    implicit none
    integer :: i, ierr, stat, nbin
    real :: const_pi
    character(len=500) :: sdsu_io_file
    
        allocate( &
        x_sbm (1:nbin) ,& ! mass per particle [g]
        den_sbm (1:nbin) ,& ! density per particle [g/cm3]
        rad_sbm (1:nbin) ,& ! radius of particle [cm]
        drad_sbm(1:nbin) ,& ! d (radius) [cm]
        vt_sbm (1:nbin) ,& ! terminal velocity [cm/s]
        brad_sbm(0:nbin) ,& ! boundary of radius bin [cm]
        stat=ierr )
    !
    ! open/read bulk masses of SBM scheme (the wrfout* file)
    sdsu_io_file = 'masses_sbm43.asc'
    const_pi=3.1415927
    open(1,file=trim(sdsu_io_file),status='old' )
    read(1,900) x_sbm%liq, x_sbm%ice_col, x_sbm%ice_pla, x_sbm%ice_den, x_sbm%snow, x_sbm%graupel, x_sbm%hail
    close(1)
    
    sdsu_io_file = 'bulkdens_sbm43.asc'
    open(1,file=trim(sdsu_io_file),status='old' )
    read(1,900) den_sbm%liq, den_sbm%ice_col, den_sbm%ice_pla, den_sbm%ice_den, den_sbm%snow, den_sbm%graupel, den_sbm%hail
    close(1)
    !
    
    
    900 format(6e13.5)
    !
    ! compute bulk radius [cm]
    do i = 1, nbin
    rad_sbm(i)%liq = ((3.e0*x_sbm(i)%liq /4.e0/const_pi/den_sbm(i)%liq )**(1.e0/3.e0))
    rad_sbm(i)%ice_col = ((3.e0*x_sbm(i)%ice_col/4.e0/const_pi/den_sbm(i)%ice_col)**(1.e0/3.e0))
    rad_sbm(i)%ice_pla = ((3.e0*x_sbm(i)%ice_pla/4.e0/const_pi/den_sbm(i)%ice_pla)**(1.e0/3.e0))
    rad_sbm(i)%ice_den = ((3.e0*x_sbm(i)%ice_den/4.e0/const_pi/den_sbm(i)%ice_den)**(1.e0/3.e0))
    rad_sbm(i)%snow = ((3.e0*x_sbm(i)%snow /4.e0/const_pi/den_sbm(i)%snow )**(1.e0/3.e0))
    rad_sbm(i)%graupel = ((3.e0*x_sbm(i)%graupel/4.e0/const_pi/den_sbm(i)%graupel)**(1.e0/3.e0))
    rad_sbm(i)%hail = ((3.e0*x_sbm(i)%hail /4.e0/const_pi/den_sbm(i)%hail )**(1.e0/3.e0))
    enddo
    !
    ! compute boundary of size bin [cm]
    do i = 0, nbin
    if(i==0) then
    brad_sbm(i)%liq = rad_sbm(1)%liq - ( rad_sbm(2)%liq - rad_sbm(1)%liq )/2.e0
    brad_sbm(i)%ice_col = rad_sbm(1)%ice_col - ( rad_sbm(2)%ice_col - rad_sbm(1)%ice_col )/2.e0
    brad_sbm(i)%ice_pla = rad_sbm(1)%ice_pla - ( rad_sbm(2)%ice_pla - rad_sbm(1)%ice_pla )/2.e0
    brad_sbm(i)%ice_den = rad_sbm(1)%ice_den - ( rad_sbm(2)%ice_den - rad_sbm(1)%ice_den )/2.e0
    brad_sbm(i)%snow = rad_sbm(1)%snow - ( rad_sbm(2)%snow - rad_sbm(1)%snow )/2.e0
    brad_sbm(i)%graupel = rad_sbm(1)%graupel - ( rad_sbm(2)%graupel - rad_sbm(1)%graupel )/2.e0
    brad_sbm(i)%hail = rad_sbm(1)%hail - ( rad_sbm(2)%hail - rad_sbm(1)%hail )/2.e0
    elseif(i<nbin) then
    brad_sbm(i)%liq = ( rad_sbm(i)%liq + rad_sbm(i+1)%liq )/2.e0
    brad_sbm(i)%ice_col = ( rad_sbm(i)%ice_col + rad_sbm(i+1)%ice_col )/2.e0
    brad_sbm(i)%ice_pla = ( rad_sbm(i)%ice_pla + rad_sbm(i+1)%ice_pla )/2.e0
    brad_sbm(i)%ice_den = ( rad_sbm(i)%ice_den + rad_sbm(i+1)%ice_den )/2.e0
    brad_sbm(i)%snow = ( rad_sbm(i)%snow + rad_sbm(i+1)%snow )/2.e0
    brad_sbm(i)%graupel = ( rad_sbm(i)%graupel + rad_sbm(i+1)%graupel )/2.e0
    brad_sbm(i)%hail = ( rad_sbm(i)%hail + rad_sbm(i+1)%hail )/2.e0
    elseif(i==nbin) then
    brad_sbm(i)%liq = rad_sbm(nbin)%liq + ( rad_sbm(nbin)%liq - rad_sbm(nbin-1)%liq )/2.e0
    brad_sbm(i)%ice_col = rad_sbm(nbin)%ice_col + ( rad_sbm(nbin)%ice_col - rad_sbm(nbin-1)%ice_col )/2.e0
    brad_sbm(i)%ice_pla = rad_sbm(nbin)%ice_pla + ( rad_sbm(nbin)%ice_pla - rad_sbm(nbin-1)%ice_pla )/2.e0
    brad_sbm(i)%ice_den = rad_sbm(nbin)%ice_den + ( rad_sbm(nbin)%ice_den - rad_sbm(nbin-1)%ice_den )/2.e0
    brad_sbm(i)%snow = rad_sbm(nbin)%snow + ( rad_sbm(nbin)%snow - rad_sbm(nbin-1)%snow )/2.e0
    brad_sbm(i)%graupel = rad_sbm(nbin)%graupel + ( rad_sbm(nbin)%graupel - rad_sbm(nbin-1)%graupel )/2.e0
    brad_sbm(i)%hail = rad_sbm(nbin)%hail + ( rad_sbm(nbin)%hail - rad_sbm(nbin-1)%hail )/2.e0
    endif
    enddo !nbin
    
    do i = 1, nbin
    drad_sbm(i)%liq = brad_sbm(i)%liq - brad_sbm(i-1)%liq
    drad_sbm(i)%ice_col = brad_sbm(i)%ice_col - brad_sbm(i-1)%ice_col
    drad_sbm(i)%ice_pla = brad_sbm(i)%ice_pla - brad_sbm(i-1)%ice_pla
    drad_sbm(i)%ice_den = brad_sbm(i)%ice_den - brad_sbm(i-1)%ice_den
    drad_sbm(i)%snow = brad_sbm(i)%snow - brad_sbm(i-1)%snow
    drad_sbm(i)%graupel = brad_sbm(i)%graupel - brad_sbm(i-1)%graupel
    drad_sbm(i)%hail = brad_sbm(i)%hail - brad_sbm(i-1)%hail
    enddo !nbin
    end subroutine set_sdsu_param