PROGRAM bedmachine_mask
  !!======================================================================
  !!                     ***  PROGRAM   bedmachine_mask ***
  !!=====================================================================
  !!  ** Purpose : extract mask variable from bedmachine file and construct
  !!               a NEMO like mask (0 on land, 1 elsewhere)
  !!
  !!  ** Method  : merge 0 and 3 mask value to ocean
  !!
  !! History :  1.0  : 02/2023  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER(KIND=4) :: narg, ijarg, iargc
  INTEGER(KIND=4) :: npiglo, npjglo
  INTEGER(KIND=4) :: ncid, id, ierr

  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mask_in, mask_out, mask_dra

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat

  CHARACTER(LEN=80) :: cf_in
  CHARACTER(LEN=80) :: cf_out = 'bedmachine_ocean_mask.nc'
  CHARACTER(LEN=80) :: cv_mask = 'mask'
  CHARACTER(LEN=80) :: cv_lon = 'nav_lon'
  CHARACTER(LEN=80) :: cv_lat = 'nav_lat'
  CHARACTER(LEN=80) :: cldum
  

  !!----------------------------------------------------------------------
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  bedmachine_mask.exe -f BEDMACHINE-file [ -o BEDMACHINE-mask]'
     PRINT *,'           [-m MASK-var ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Create a NEMO-like ocean mask (1=ocean, 0=land) from a bedmachine'
     PRINT *,'      mask variable.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f BEDMACHINE-file : specify the name of the bedmachine file.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -o BEDMACHINE-mask : specify the name of the bedmachine mask out.' 
     PRINT *,'       -m MASK-var : specify the name of the mask variable'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' or name given by -o option'
     PRINT *,'         variables : ', TRIM(cv_mask),' ( 1 or 0 )'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE ( '-m'   ) ; CALL getarg(ijarg, cv_mask    ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ierr=NF90_OPEN(cf_in, NF90_NOWRITE, ncid)
print *,NF90_STRERROR(ierr)
  ierr=NF90_INQ_DIMID( ncid, 'x', id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
print *,NF90_STRERROR(ierr)
  ierr=NF90_INQ_DIMID( ncid, 'y', id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)
print *,NF90_STRERROR(ierr)

  ALLOCATE( mask_in(npiglo, npjglo), mask_out(npiglo, npjglo) , mask_dra(npiglo, npjglo))
  ALLOCATE( rlon(npiglo, npjglo), rlat(npiglo, npjglo) )

  ierr = NF90_INQ_VARID(ncid,cv_mask,id) ; ierr = NF90_GET_VAR(ncid,id,mask_in)
print *,NF90_STRERROR(ierr)
  ierr = NF90_INQ_VARID(ncid,cv_lon,id)  ; ierr = NF90_GET_VAR(ncid,id,rlon   )
print *,NF90_STRERROR(ierr)
  ierr = NF90_INQ_VARID(ncid,cv_lat,id)  ; ierr = NF90_GET_VAR(ncid,id,rlat   )
print *,NF90_STRERROR(ierr)
  ierr = NF90_CLOSE(ncid)
print *,NF90_STRERROR(ierr)

  mask_out=0
  mask_dra=0
  WHERE (mask_in == 0 .OR. mask_in == 3 ) mask_out=1  ! Ocean and floating ice
  WHERE (mask_in == 3 )                   mask_dra=1  ! where there is floating ice

  CALL WriteOutput( cf_out)

  CONTAINS
  SUBROUTINE WriteOutput ( cd_fout)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE WriteOutput  ***
    !!
    !! ** Purpose :  Create and write output file
    !!
    !! ** Method  :  Use NF90 primitive
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_fout

    INTEGER(KIND=4) :: ncid, id, idx,idy, ierr, idlon, idlat, iddra

    
    !!----------------------------------------------------------------------
    ierr=NF90_CREATE(cd_fout,NF90_NETCDF4, ncid)
    ! Dimensions
    ierr = NF90_DEF_DIM(ncid,'x',npiglo, idx)
    ierr = NF90_DEF_DIM(ncid,'y',npjglo, idy)
    ! Variables
    ierr = NF90_DEF_VAR(ncid,cv_lon ,NF90_FLOAT, (/idx,idy/) ,idlon, deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,cv_lat ,NF90_FLOAT, (/idx,idy/) ,idlat, deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,cv_mask,NF90_INT  , (/idx,idy/) ,id   , deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,'isfmask',NF90_INT, (/idx,idy/) ,iddra, deflate_level=1,chunksizes=(/1,npjglo/10/) )

    ierr = NF90_ENDDEF(ncid)
    ierr = NF90_PUT_VAR( ncid, idlon, rlon)
    ierr = NF90_PUT_VAR( ncid, idlat, rlat)
    ierr = NF90_PUT_VAR( ncid, id   , mask_out )
    ierr = NF90_PUT_VAR( ncid, iddra, mask_dra )
    ierr = NF90_CLOSE(ncid)

    END SUBROUTINE WriteOutput

END PROGRAM bedmachine_mask
