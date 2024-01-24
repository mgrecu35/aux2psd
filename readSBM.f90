module sdsu_io
    type particle_sbm !SBM paticle type
    real :: liq,ice_col,ice_pla,ice_den,snow,graupel,hail
    end type particle_sbm
    
    type ( particle_sbm ), allocatable, dimension(:,:,:) :: & !3D parameters (is:ie,js:je,ks:ke)
    q_sbm    ,& ! particle mixing ratio [g/m3]
    re_sbm      ! particle effective radius [micron]

    type ( particle_sbm ), allocatable, dimension(:,:) :: & !2D parameters (is:ie,js:je)
    qcol_sbm  ! column integrated (equivalent water path) particle amount [kg/m2]

    type ( particle_sbm ), allocatable, dimension(:) :: & !size bin parameters (1:nbin)
    x_sbm   ,& ! mass per particle [g]
    den_sbm ,& ! density per particle [g/cm3]
    rad_sbm ,& ! radius of particle [cm]
    drad_sbm,& ! d (radius) [cm]
    vt_sbm     ! terminal velocity [cm/s]

    type ( particle_sbm ), allocatable, dimension(:) :: & !size bin parameters (0:nbin)
    brad_sbm  ! boundary of radius bin [cm]
end module sdsu_io



!n_sbm%liq = 0.e0 ; q_sbm%liq = 0.e0
!tag_char = 'dr'  ;  sdsu_io_file =trim(sdsu_dir_input)//'auxhist2'//trim(sdsu_inp_name(efile_len-23:efile_len))
!call check( nf90_open(trim(sdsu_io_file), nf90_nowrite, ncidsbm) ) !open nc file

!do n = 1, nbin
!    write(char_bin,"(I2.2)") n  ; para_char = 'ffc'//char_bin//tag_char  
!    call check( nf90_inq_varid(ncidsbm, trim(para_char) , varid ) )       ! SBM PSD function [/ (cm3 g)]
!    call check( nf90_get_var(ncidsbm, varid, net3d(is:ie,js:je,ks:ke), start=(/is,js,ks/), count=(/di,dj,dk/)))
!    do k = 1, mxlyr ; do j = myj_start, myj_end ; do i = myi_start, myi_end
!       call sbm_function(x_sbm(n)%liq, drad_sbm(n)%liq*1.e-2, net3d(i,j,k), &
!                         n_sbm(i,j,k,n)%liq, q_sbm(i,j,k)%liq )
!    enddo ; enddo ; enddo 
! enddo   

subroutine sbm_function43(rho,x,dr,m,n,q)
    implicit none
   !--------------------------------------------------------------------------------------------
   ! Comments:  
   !  Convert HUCM SBM (43bin version) PSD mass mixing ratio to common mass mixing (or number 
   !  concentrations) per volume. 
   !
   ! History:
   !  12/2010  Toshi Matsui@NASA GSFC ; Initial   
   !           
   ! References: 
   !-----------------------------------------------------------------------------------------------------
    real,intent(in) :: rho    ! dry air density [kg/m3]
    real,intent(in) :: x      ! mass per particle [g]
    real,intent(in) :: dr     ! SBM size-bin width [m]
    real,intent(in) :: m      ! SBM PSD mass mixing ratio [kg/kg]
    real,intent(out) :: n     ! particle size density [1/m4]
    real,intent(inout) :: q   ! total mixing ratio [g/m3] 
    real :: q_bin  ! mixing ratio per one bin [g/m3]
    real :: n_bin  ! drop size ditributions [1/m3]
   
   !
   ! initialize
   !
    n_bin = 0.e0 ; q_bin = 0.e0
   
    q_bin = m * rho * 1.e3      ! mass conc per bin [g/m3]
   
    n_bin = q_bin / x / dr      ! # conc per bin [#/m4]
   
    if( isnan(m) ) then
        print*,'MSG sbm_function43, find NaN in m' 
    endif
   
   !
   ! clean numerical noise
   !
    if( n_bin <= 0.e0 ) then
        n_bin = 0.e0  ; q_bin = 0.e0
    endif
   
   !
   ! update PSD and total mixing ratio
   !
    n = n_bin     ! particle size density [1/m4]
    q = q + q_bin ! total mixing ratio [g/m3]
   
   
    return
    end subroutine sbm_function43

    subroutine sbm_function43_dm(rho,x,dr,r,m,n,q,dm)
        implicit none
       !--------------------------------------------------------------------------------------------
       ! Comments:  
       !  Convert HUCM SBM (43bin version) PSD mass mixing ratio to common mass mixing (or number 
       !  concentrations) per volume. 
       !
       ! History:
       !  12/2010  Toshi Matsui@NASA GSFC ; Initial   
       !           
       ! References: 
       !-----------------------------------------------------------------------------------------------------
        real,intent(in) :: rho    ! dry air density [kg/m3]
        real,intent(in) :: x      ! mass per particle [g]
        real,intent(in) :: dr     ! SBM size-bin width [m]
        real,intent(in) :: r     ! SBM size-bin width [cm]
        real,intent(in) :: m      ! SBM PSD mass mixing ratio [kg/kg]
        real,intent(out) :: n     ! particle size density [1/m4]
        real,intent(inout) :: q   ! total mixing ratio [g/m3] 
        real,intent(inout) :: dm   ! total mixing ratio [g/m3]
        real :: q_bin  ! mixing ratio per one bin [g/m3]
        real :: n_bin  ! drop size ditributions [1/m3]
       
       !
       ! initialize
       !
        n_bin = 0.e0 ; q_bin = 0.e0
       
        q_bin = m * rho * 1.e3      ! mass conc per bin [g/m3]
       
        n_bin = q_bin / x / dr      ! # conc per bin [#/m4]
       
        if( isnan(m) ) then
            print*,'MSG sbm_function43, find NaN in m' 
        endif
       
       !
       ! clean numerical noise
       !
        if( n_bin <= 0.e0 ) then
            n_bin = 0.e0  ; q_bin = 0.e0
        endif
       
       !
       ! update PSD and total mixing ratio
       !
        n = n_bin     ! particle size density [1/m4]
        q = q + q_bin ! total mixing ratio [g/m3]
        dm = dm + q_bin * 2*r  ! total weighted  dm [g/m3*cm]
       
       
        return
    end subroutine sbm_function43_dm
    
    
subroutine sbm_function(x,dr,f,n,q)
    implicit none
   !--------------------------------------------------------------------------------------------
   ! Comments:  
   !  Convert HUCM SBM (33bin version) PSD function to common mass mixing ration (or number 
   !  concentrations) per volume. 
   !
   ! History:
   !  12/2010  Toshi Matsui@NASA GSFC ; Initial   
   !           
   ! References: 
   !-----------------------------------------------------------------------------------------------------
    real,intent(in) :: x      ! mass per particle [g]
    real,intent(in) :: dr     ! SBM size-bin width [m]
    real,intent(in) :: f      ! SBM PSD function [/ (cm3 g)]
    real,intent(out) :: n     ! particle size density [1/m4]
    real,intent(inout) :: q   ! total mixing ratio [g/m3] 
    real :: q_bin  ! mixing ratio per one bin [g/m3]
    real :: n_bin  ! drop size ditributions [1/m3]
   
   !
   ! initialize
   !
    n_bin = 0.e0 ; q_bin = 0.e0
   
    q_bin = f*3.0e0*0.23105e0*x*x*1.0e6 ! mass conc per bin [g/m3]
   
    n_bin = q_bin / x / dr              ! # conc per bin [#/m4]
   
   !
   ! clean numerical noise
   !
    if( n_bin <= 0.e0 ) then
        n_bin = 0.e0  ; q_bin = 0.e0
    endif
   
   !
   ! update PSD and total mixing ratio
   !
    n = n_bin     ! particle size density [1/m4]
    q = q + q_bin ! total mixing ratio [g/m3]
   
   
    return
end subroutine sbm_function


!program test_sdsu_param
!integer nbin

!nbin=43

!call set_sdsu_param(nbin)

!end program test_sdsu_param
