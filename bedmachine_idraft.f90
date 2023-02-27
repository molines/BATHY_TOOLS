PROGRAM bedmachine_idraft
  !!======================================================================
  !!                     ***  PROGRAM   bedmachine_idraft ***
  !!=====================================================================
  !!  ** Purpose : Compute ice draft from surface and ice thickness variable
  !!
  !!  ** Method  : read compute write ...
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
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mask, mask_oce

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rlon, rlat
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rbed,rsurf,rthick,rdraft

  CHARACTER(LEN=80) :: cf_in
  CHARACTER(LEN=80) :: cf_out = 'bedmachine_ice_draf'
  CHARACTER(LEN=80) :: cv_bed   = 'bed'
  CHARACTER(LEN=80) :: cv_surf  = 'surface'
  CHARACTER(LEN=80) :: cv_thic  = 'thickness'
  CHARACTER(LEN=80) :: cv_mask  = 'mask'
  CHARACTER(LEN=80) :: cv_draft = 'icedraft'
  CHARACTER(LEN=80) :: cv_lon = 'nav_lon'
  CHARACTER(LEN=80) :: cv_lat = 'nav_lat'
  CHARACTER(LEN=80) :: cldum
  

  !!----------------------------------------------------------------------
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  bedmachine_idraft.exe -f BEDMACHINE-file [ -o BEDMACHINE-draft]'
     PRINT *,'         [-bed BEDMACHINE-bed] [-surf BEDMACHINE-surf] [-mask BEDMACHINE-mask]'
     PRINT *,'         [-thick BEDMACHINE-thickness] [-draft BEDMACHINE-draft]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Create a ice draft file from the bedmachine input file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f BEDMACHINE-file : specify the name of the bedmachine file.' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -o BEDMACHINE-draft : specify the name of the bedmachine-draft file.' 
     PRINT *,'       -bed BEDMACHINE-bed : name of variable (in) for bed elevation.'
     PRINT *,'       -surf BEDMACHINE-surf : name of variable (in) for surface elevation.'
     PRINT *,'       -thick BEDMACHINE-thicness : name of variable (in) for ice thickness.'
     PRINT *,'       -mask MASK-var : name of the variable (in) for mask.'
     PRINT *,'       -draft BEDMACHINE-draft:  name of variable (out) for ice draft.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' or name given by -o option'
     PRINT *,'         variables : ', TRIM(cv_draft),' ( 1 or 0 )'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'    ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE ( '-bed'  ) ; CALL getarg(ijarg, cv_bed     ) ; ijarg=ijarg+1
     CASE ( '-surf' ) ; CALL getarg(ijarg, cv_surf    ) ; ijarg=ijarg+1
     CASE ( '-thick') ; CALL getarg(ijarg, cv_thic    ) ; ijarg=ijarg+1
     CASE ( '-mask' ) ; CALL getarg(ijarg, cv_mask    ) ; ijarg=ijarg+1
     CASE ( '-draft') ; CALL getarg(ijarg, cv_draft   ) ; ijarg=ijarg+1
     CASE DEFAULT     ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ierr=NF90_OPEN(cf_in, NF90_NOWRITE, ncid)
  ierr=NF90_INQ_DIMID( ncid, 'x', id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr=NF90_INQ_DIMID( ncid, 'y', id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)

  ALLOCATE( rbed(npiglo, npjglo), rsurf(npiglo, npjglo) ,rthick(npiglo, npjglo) )
  ALLOCATE( rdraft(npiglo, npjglo), rlon(npiglo, npjglo), rlat(npiglo, npjglo) )
  ALLOCATE( mask(npiglo, npjglo), mask_oce(npiglo, npjglo))

  ierr = NF90_INQ_VARID(ncid,cv_bed ,id)  ; ierr = NF90_GET_VAR(ncid,id,rbed   )
  ierr = NF90_INQ_VARID(ncid,cv_surf,id)  ; ierr = NF90_GET_VAR(ncid,id,rsurf  )
  ierr = NF90_INQ_VARID(ncid,cv_thic,id)  ; ierr = NF90_GET_VAR(ncid,id,rthick )
  ierr = NF90_INQ_VARID(ncid,cv_mask,id)  ; ierr = NF90_GET_VAR(ncid,id,mask   )
  ierr = NF90_INQ_VARID(ncid,cv_lon ,id)  ; ierr = NF90_GET_VAR(ncid,id,rlon   )
  ierr = NF90_INQ_VARID(ncid,cv_lat ,id)  ; ierr = NF90_GET_VAR(ncid,id,rlat   )
  ierr = NF90_CLOSE(ncid)

  mask_oce=0
  WHERE (mask == 0 .OR. mask == 3 ) mask_oce=1
  rdraft=(rsurf - rthick) * mask_oce  + rbed * ( 1 - mask_oce) 
  !     isfd = (hsurf - hice) * msk + bed * (1 - msk)
  ! test =  rsurf - rthick -rbed 

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

    INTEGER(KIND=4) :: ncid, id, idx,idy, ierr, idlon, idlat

    
    !!----------------------------------------------------------------------
    ierr=NF90_CREATE(cd_fout,NF90_NETCDF4, ncid)
    ! Dimensions
    ierr = NF90_DEF_DIM(ncid,'x',npiglo, idx)
    ierr = NF90_DEF_DIM(ncid,'y',npjglo, idy)
    ! Variables
    ierr = NF90_DEF_VAR(ncid,cv_lon ,NF90_FLOAT, (/idx,idy/) ,idlon, deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,cv_lat ,NF90_FLOAT, (/idx,idy/) ,idlat, deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,cv_draft,NF90_FLOAT, (/idx,idy/) ,id   , deflate_level=1,chunksizes=(/1,npjglo/10/) )

    ierr = NF90_ENDDEF(ncid)
    ierr = NF90_PUT_VAR( ncid, idlon, rlon)
    ierr = NF90_PUT_VAR( ncid, idlat, rlat)
    ierr = NF90_PUT_VAR( ncid, id   , rdraft )
    ierr = NF90_CLOSE(ncid)

    END SUBROUTINE WriteOutput

END PROGRAM bedmachine_idraft
