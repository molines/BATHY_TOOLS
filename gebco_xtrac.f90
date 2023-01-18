PROGRAM gebco_tool
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
  INTEGER(KIND=4) :: ji,jj, jblk  ! dummy loop index
  INTEGER(KIND=4) :: ij1,ij2
  INTEGER(KIND=4) :: iimin,iimax
  INTEGER(KIND=4) :: ijmin,ijmax
  INTEGER(KIND=4) :: iblksz=100
  INTEGER(KIND=4) :: npi,npj
  INTEGER(KIND=4) :: ncid, id, ierr

  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ielev, imask

  REAL(KIND=8)   :: dresol=15 ! in second of arc
  REAL(KIND=8)   :: dlonmin, dlonmax, dlatmin,dlatmax
  REAL(KIND=8)   :: dlon0, dlat0

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dl_lon, dl_lat

  CHARACTER(LEN=80) :: cf_bathy='GEBCO_2022.nc'
  CHARACTER(LEN=80) :: cf_mask='OSM_land_to_gebco.nc'
  CHARACTER(LEN=80) :: cv_bathy='elevation'
  CHARACTER(LEN=80) :: cv_mask='Band1'
  CHARACTER(LEN=80) :: clon='lon'
  CHARACTER(LEN=80) :: clat='lat'
  CHARACTER(LEN=80) :: czon
  !!
  !!----------------------------------------------------------------------
  !! GEBCOTOOLS_1.0 , MEOM 2023
  !! Copyright (c) 2023, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

  dlonmin=33.2 ; dlonmax=36.2
  dlatmin=45.2 ; dlatmax=46.4
  czon='East-Azov'

  CALL xtrac (cf_bathy,cv_bathy,czon) 
  czon=TRIM(czon)//'_mask'
  CALL xtrac (cf_mask, cv_mask, czon)

CONTAINS
 SUBROUTINE xtrac( cd_fin, cd_vin, cd_zon)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE xtrac  ***
    !!
    !! ** Purpose :  xtrac zone from global file
    !!
    !! ** Method  :   
    !!
    !!----------------------------------------------------------------------

    CHARACTER(LEN=80), INTENT(in) :: cd_fin ! input file
    CHARACTER(LEN=80), INTENT(in) :: cd_vin ! input variable
    CHARACTER(LEN=80), INTENT(in) :: cd_zon ! zone name (for output file 
 !=====
  ierr = NF90_OPEN(cd_fin,NF90_NOWRITE,ncid)
  ierr = NF90_INQ_DIMID(ncid,clon,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npi)
  ierr = NF90_INQ_DIMID(ncid,clat,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npj)
  PRINT *, '  LON = ', npi
  PRINT *, '  LAT = ', npj

! ALLOCATE ( ielev(npi,npj), imask(npi,npj) )
  ALLOCATE ( dl_lon(npi), dl_lat(npj) )
  
  PRINT *, 'ielev, mask and lon/lat allocated'

 ierr = NF90_INQ_VARID(ncid, clon, id) 
 ierr = NF90_GET_VAR(ncid, id, dl_lon)
 ierr = NF90_INQ_VARID(ncid, clat, id) 
 ierr = NF90_GET_VAR(ncid, id, dl_lat)
 dlon0=dl_lon(1)
 dlat0=dl_lat(1)
 print *, dlon0, dlat0
 print *, dlonmin-dlon0, dlatmin-dlat0

 iimin=(dlonmin-dlon0)*60*60/dresol +1
 iimax=(dlonmax-dlon0)*60*60/dresol +1
 ijmin=(dlatmin-dlat0)*60*60/dresol +1
 ijmax=(dlatmax-dlat0)*60*60/dresol +1

 print *, iimin, iimax, ijmin, ijmax

  ALLOCATE ( ielev(iimax-iimin+1, ijmax-ijmin+1) )
  ierr = NF90_INQ_VARID(ncid, cd_vin, id) 
  ierr = NF90_GET_VAR(ncid, id, ielev(:,:),start=(/iimin,ijmin/), count=(/iimax-iimin+1,ijmax-ijmin+1/) )
  ierr = NF90_CLOSE(ncid)

!  print *, ielev

  CALL geb_wri(cd_zon,cd_vin,ielev)
  DEALLOCATE( ielev, dl_lon, dl_lat)
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

    INTEGER(KIND=4)   :: ierr,ncid, id, idlon,idlat
    INTEGER(KIND=4)   :: idvlon, idvlat
    CHARACTER(LEN=80) :: cf_out

    cf_out=TRIM(cd_zon)//'.nc'
    ! Create file
    ierr = NF90_CREATE(cf_out,NF90_NETCDF4,ncid)
    ! Create dimensions
    ierr = NF90_DEF_DIM(ncid,clon,iimax-iimin+1,idlon)
    ierr = NF90_DEF_DIM(ncid,clat,ijmax-ijmin+1,idlat)
    ! Create Variable
    ierr = NF90_DEF_VAR(ncid, clon,     NF90_DOUBLE,(/idlon/),       idvlon, deflate_level=1)
    ierr = NF90_DEF_VAR(ncid, clat,     NF90_DOUBLE,(/idlat/),       idvlat, deflate_level=1)
    ierr = NF90_DEF_VAR(ncid, cd_vin,   NF90_SHORT, (/idlon,idlat/), id,     deflate_level=1)
    ! Create Attibute
    ! later ...
    ierr = NF90_ENDDEF(ncid)
    ! put variables
    ierr = NF90_PUT_VAR(ncid, idvlon, dl_lon(iimin:iimax) )
    ierr = NF90_PUT_VAR(ncid, idvlat, dl_lat(ijmin:ijmax) )
    ierr = NF90_PUT_VAR(ncid, id,     kelev               )

    ierr = NF90_CLOSE(ncid)


 END SUBROUTINE geb_wri
  
END PROGRAM gebco_tool
