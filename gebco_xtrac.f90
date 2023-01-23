PROGRAM gebco_xtrac
  !!======================================================================
  !!                     ***  PROGRAM  gebco_tool  ***
  !!=====================================================================
  !!  ** Purpose : a program to handle GEBCO big files
  !!
  !!  ** Method  : netcdf primitives
  !!
  !! History :  1.0  : 01/2023  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4) :: narg, ijarg
  INTEGER(KIND=4) :: ji,jj, jblk  ! dummy loop index
  INTEGER(KIND=4) :: ij1,ij2
  INTEGER(KIND=4) :: iimin,iimax
  INTEGER(KIND=4) :: ijmin,ijmax
  INTEGER(KIND=4) :: iblksz=100
  INTEGER(KIND=4) :: npi,npj
  INTEGER(KIND=4) :: ncid, id, ierr

  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ielev, imask, ielevcor,ifail

  REAL(KIND=8)   :: dresol=15 ! in second of arc
  REAL(KIND=8)   :: dlonmin, dlonmax, dlatmin,dlatmax
  REAL(KIND=8)   :: dlon0, dlat0

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dl_lon, dl_lat
  REAL(KIND=8)                            :: dl_sign=-1.0d0

  CHARACTER(LEN=80) :: cf_bathy='GEBCO_2022.nc'
  CHARACTER(LEN=80) :: cf_mask='OSM_land_to_gebco.nc'
  CHARACTER(LEN=80) :: cv_bathy='elevation'
  CHARACTER(LEN=80) :: cv_mask='Band1'
  CHARACTER(LEN=80) :: clon='lon'
  CHARACTER(LEN=80) :: clat='lat'
  CHARACTER(LEN=80) :: czon
  CHARACTER(LEN=80) :: cldum

  LOGICAL           :: lglobal    =.FALSE.
  LOGICAL           :: lwij       =.FALSE.
  LOGICAL           :: lwlonlat   =.FALSE.
  !!
  !!----------------------------------------------------------------------
  !! GEBCOTOOLS_1.0 , MEOM 2023
  !! Copyright (c) 2023, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  gebco_xtrac.exe [-z PREDEF-zone ] [-g] [-wij imin imax jmin jmax]'
     PRINT *,'                       [-wlonlat lonmin lonmax latmin latmax] [-nam NAME-zone]'
     PRINT *,'                       [-b BATHY-file] [-neg]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Extract a rectangular zone from '//TRIM(cf_bathy)//' and '//TRIM(cf_mask)
     PRINT *,'        For each extracted zone, this program produce 4 files :'
     PRINT *,"           * <ZONE>.nc  : the extracted bathymetry"
     PRINT *,"           * <ZONE>_mask.nc: the extracted OSM mask"
     PRINT *,"           * <ZONE>_mask_corrected.nc: the extracted masked bathymetry"
     PRINT *,"           * <ZONE>_mask_corrected_fail.nc:  the mask of land points falling "
     PRINT *,"                    into the ocean."
     PRINT *,"                    These latter points are errors in the data base. This file"
     PRINT *,"                    will be used by sosie drowning procedure to fix the errors."
     PRINT *,'      '
     PRINT *,'     ARGUMENTS /OPTIONS:'
     PRINT *,'        -z PREDEF-zone: give a name of a predefined extraction zone. '
     PRINT *,'             Available so for : East-Azov  caribe GOM Bahamas Seychelles Brest'
     PRINT *,'        -g  : take the global domain (huge files ! )'
     PRINT *,'        -wij imin imax jmin jmax : specify an extration windows using I-J index'
     PRINT *,'        -wlonlat lonmin lonmax latmin latmax : specify an extration windows '
     PRINT *,'             using longitude and latitudes. Longitudes are between -180 and 180.'
     PRINT *,'        -nam NAME-zone : give a name for the extracted region (used for '
     PRINT *,'             defining output file names.'
     PRINT *,'        -b BATHY-file : use BATHY-file instead of '//TRIM(cf_bathy)
     PRINT *,'        -neg : do not change sign of bathymetry'
     PRINT *,'      '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-z'   ) ; CALL getarg(ijarg, czon ) ; ijarg=ijarg+1
        lwlonlat = .TRUE.
        SELECT CASE( czon)
        CASE ( 'East-Azov')
           dlonmin=33.2 ; dlonmax=36.2
           dlatmin=45.2 ; dlatmax=46.4
        CASE ( 'caribe')
           dlonmin=-80 ; dlonmax=-55
           dlatmin=5 ; dlatmax=30
        CASE ( 'GOM')
           dlonmin=-105 ; dlonmax=-70
           dlatmin=5 ; dlatmax=40
        CASE ( 'Bahamas')
           dlonmin=-82.2 ; dlonmax=-72.3
           dlatmin=23.1 ; dlatmax=27.3
        CASE ( 'Seychelles')
           dlonmin=42.5 ; dlonmax=62.0
           dlatmin=-26.5 ; dlatmax=-2.0
        CASE ( 'Brest')
           dlonmin=-5.19 ; dlonmax=-3.79
           dlatmin=47.68 ; dlatmax=48.73
        CASE DEFAULT 
           PRINT *, TRIM(czon),' NOT A PREDEFINED ZONE !'
           STOP
        END SELECT
     CASE ( '-g' )
        lglobal= .TRUE.
        lwij   = .TRUE.
     CASE ( '-wij'     ) ; lwij   = .TRUE.
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) iimin
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) iimax
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) ijmin
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) ijmax
     CASE ( '-wlonlat' ) ; lwlonlat = .TRUE.
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) dlonmin
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) dlonmax
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) dlatmin
                         ; CALL getarg(ijarg,cldum) ; ijarg=ijarg+1 ; READ(cldum,*) dlatmax
                         ; PRINT *, dlonmin, dlonmax, dlatmin, dlatmax
     CASE ( '-nam'     ) ; CALL getarg(ijarg,czon ) ; ijarg=ijarg+1
     CASE ( '-b'       ) ; CALL getarg(ijarg,cf_bathy ) ; ijarg=ijarg+1
     CASE ( '-neg'     ) ; dl_sign=1.0d0
           
        ! option
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO




  CALL getsize(cf_bathy,npi,npj)
  ALLOCATE ( dl_lon(npi), dl_lat(npj) ) 
  CALL getlonlat(cf_bathy)

  CALL xtrac (cf_bathy,cv_bathy,czon,ielev) 
  czon=TRIM(czon)//'_mask'
  CALL xtrac (cf_mask, cv_mask, czon,imask)
  WHERE (imask == 1 ) imask=0
  WHERE (imask == -1) imask=1
  czon=TRIM(czon)//'_corrected'

  ielevcor=dl_sign*ielev*imask
  CALL geb_wri(czon,cv_bathy,ielevcor)
  ifail=1
  WHERE (ielevcor < 0  ) ifail = 0
  czon=TRIM(czon)//'_fail'
  CALL geb_wri(czon,cv_bathy,ifail)

CONTAINS

 SUBROUTINE getsize(cd_fil,kpi,kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE getsize  ***
    !!
    !! ** Purpose :   get horizontal size from file 
    !!
    !! ** Method  :   
    !!
    !! References :  
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_fil
    INTEGER(KIND=4), INTENT(out) :: kpi,kpj
    
    INTEGER(KIND=4) :: ierr, icid, id
    ierr = NF90_OPEN(cd_fil,NF90_NOWRITE,icid)
    ierr = NF90_INQ_DIMID(icid,clon,id ) ; ierr = NF90_INQUIRE_DIMENSION(icid, id, len=kpi)
    ierr = NF90_INQ_DIMID(icid,clat,id ) ; ierr = NF90_INQUIRE_DIMENSION(icid, id, len=kpj)
    PRINT *, '  LON = ', kpi
    PRINT *, '  LAT = ', kpj
    ierr = NF90_CLOSE(icid)
    IF (lglobal) THEN
       iimin=1 ; iimax=npi
       ijmin=1 ; ijmax=npj
    ENDIF
      

 END SUBROUTINE getsize
SUBROUTINE getlonlat(cd_fil)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE getlonlat  ***
    !!
    !! ** Purpose :   get lon lat 
    !!
    !! ** Method  :   
    !!
    !! References :  
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_fil

    INTEGER(KIND=4) :: ierr, icid, id
    ierr = NF90_OPEN(cd_fil,NF90_NOWRITE,icid)
    ierr = NF90_INQ_VARID(icid, clon, id) 
    ierr = NF90_GET_VAR(icid, id, dl_lon)
    ierr = NF90_INQ_VARID(icid, clat, id) 
    ierr = NF90_GET_VAR(icid, id, dl_lat)
    ierr = NF90_CLOSE(icid)

    dlon0=dl_lon(1)
    dlat0=dl_lat(1)
    IF ( lwij ) THEN
       dlonmin = dl_lon(iimin)
       dlonmax = dl_lon(iimax)
       dlatmin = dl_lat(ijmin)
       dlatmax = dl_lat(ijmax)
       PRINT *, dlonmin, dlonmax, dlatmin, dlatmax
    ELSE IF (lwlonlat) THEN
       iimin=(dlonmin-dlon0)*60*60/dresol +1
       iimax=(dlonmax-dlon0)*60*60/dresol +1
       ijmin=(dlatmin-dlat0)*60*60/dresol +1
       ijmax=(dlatmax-dlat0)*60*60/dresol +1
       PRINT *, iimin, iimax, ijmin, ijmax
    ELSE
       PRINT *,' Weird ! lwij or lwlonlat shoudl be true ...'
       STOP
    ENDIF
    ALLOCATE ( ielev(iimax-iimin+1, ijmax-ijmin+1) )
    ALLOCATE ( imask(iimax-iimin+1, ijmax-ijmin+1) )
    ALLOCATE ( ielevcor(iimax-iimin+1, ijmax-ijmin+1) )
    ALLOCATE ( ifail(iimax-iimin+1, ijmax-ijmin+1) )

 END SUBROUTINE getlonlat


 SUBROUTINE xtrac( cd_fin, cd_vin, cd_zon,kdta)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE xtrac  ***
    !!
    !! ** Purpose :  xtrac zone from global file
    !!
    !! ** Method  :   
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=2), DIMENSION(:,:), INTENT(inout)  :: kdta
    CHARACTER(LEN=80), INTENT(in) :: cd_fin ! input file
    CHARACTER(LEN=80), INTENT(in) :: cd_vin ! input variable
    CHARACTER(LEN=80), INTENT(in) :: cd_zon ! zone name (for output file 
 !=====
  ierr = NF90_OPEN(cd_fin,NF90_NOWRITE,ncid)

  ierr = NF90_INQ_VARID(ncid, cd_vin, id) 
  ierr = NF90_GET_VAR(ncid, id, kdta(:,:),start=(/iimin,ijmin/), count=(/iimax-iimin+1,ijmax-ijmin+1/) )
  ierr = NF90_CLOSE(ncid)

  CALL geb_wri(cd_zon,cd_vin,kdta)
 END SUBROUTINE xtrac

 SUBROUTINE geb_wri(cd_zon,cd_vin, kelev)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE geb_wri  ***
    !!
    !! ** Purpose :  Write extraction to netcdf file 
    !!
    !! ** Method  :  takes global variables from main
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=2), DIMENSION(:,:), INTENT(in) :: kelev
    CHARACTER(LEN=80), INTENT(in) :: cd_zon
    CHARACTER(LEN=80), INTENT(in) :: cd_vin

    INTEGER(KIND=4)   :: ierr,ncid, id, idlon,idlat,idtim
    INTEGER(KIND=4)   :: idvlon, idvlat,idvtim
    CHARACTER(LEN=80) :: cf_out

    REAL(KIND=8), DIMENSION(1) :: dl_tim

    cf_out=TRIM(cd_zon)//'.nc'
    dl_tim(1) = 0.d0
    ! Create file
    ierr = NF90_CREATE(cf_out,NF90_NETCDF4,ncid)
    ! Create dimensions
    ierr = NF90_DEF_DIM(ncid,clon,iimax-iimin+1,idlon)
    ierr = NF90_DEF_DIM(ncid,clat,ijmax-ijmin+1,idlat)
    ierr = NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,idtim)
    ! Create Variable
    ierr = NF90_DEF_VAR(ncid, clon,     NF90_DOUBLE,(/idlon/),       idvlon, deflate_level=1)
    ierr = NF90_DEF_VAR(ncid, clat,     NF90_DOUBLE,(/idlat/),       idvlat, deflate_level=1)
    ierr = NF90_DEF_VAR(ncid, 'time',   NF90_DOUBLE,(/idtim/),       idvtim                 )
    ierr = NF90_DEF_VAR(ncid, cd_vin,   NF90_SHORT, (/idlon,idlat,idtim/), id, deflate_level=1)
    ! Create Attibute
    ierr = NF90_PUT_ATT(ncid,id,'missing_value',0)
    ! later ...
    ierr = NF90_ENDDEF(ncid)
    ! put variables
    ierr = NF90_PUT_VAR(ncid, idvlon, dl_lon(iimin:iimax) )
    ierr = NF90_PUT_VAR(ncid, idvlat, dl_lat(ijmin:ijmax) )
    ierr = NF90_PUT_VAR(ncid, idvtim, dl_tim(1) )
    ierr = NF90_PUT_VAR(ncid, id,     kelev               )

    ierr = NF90_CLOSE(ncid)


 END SUBROUTINE geb_wri
  
END PROGRAM gebco_xtrac
