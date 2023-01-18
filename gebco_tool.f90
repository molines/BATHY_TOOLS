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
  INTEGER(KIND=4) :: iblksz=100
  INTEGER(KIND=4) :: npi,npj
  INTEGER(KIND=4) :: ncid, id, ierr

  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: ielev, imask

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dl_lon, dl_lat

  CHARACTER(LEN=80) :: cf_bathy='GEBCO_2022.nc'
  CHARACTER(LEN=80) :: cf_mask='OSM_land_to_gebco.nc'
  CHARACTER(LEN=80) :: cv_bathy='elevation'
  CHARACTER(LEN=80) :: cv_mask='band1'
  CHARACTER(LEN=80) :: clon='lon'
  CHARACTER(LEN=80) :: clat='lat'
  !!
  !!----------------------------------------------------------------------
  !! GEBCOTOOLS_1.0 , MEOM 2023
  !! Copyright (c) 2023, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

  ierr = NF90_OPEN(cf_bathy,NF90_NOWRITE,ncid)
  ierr = NF90_INQ_DIMID(ncid,clon,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npi)
  ierr = NF90_INQ_DIMID(ncid,clat,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npj)
  PRINT *, '  LON = ', npi
  PRINT *, '  LAT = ', npj

! ALLOCATE ( ielev(npi,npj), imask(npi,npj) )
! ALLOCATE ( dl_lon(npi), dl_lat(npj) )
  ALLOCATE ( ielev(npi,iblksz) )
  
  PRINT *, 'ielev, mask and lon/lat allocated'

! ierr = NF90_INQ_VARID(ncid, clon, id) 
! ierr = NF90_GET_VAR(ncid, id, dl_lon)
! ierr = NF90_INQ_VARID(ncid, clat, id) 
! ierr = NF90_GET_VAR(ncid, id, dl_lat)
  ierr = NF90_INQ_VARID(ncid, cv_bathy, id) 
  ij1=1
  DO jblk=1,npj/iblksz
    ij2=ij1+iblksz-1
    print *, jblk, ij1,ij2
    ierr = NF90_GET_VAR(ncid, id, ielev(:,:),start=(/1,ij1/), count=(/npi,iblksz/) )
    ij1=ij2+1
  ENDDO
  

  ierr = NF90_CLOSE(ncid)
  
END PROGRAM gebco_tool
