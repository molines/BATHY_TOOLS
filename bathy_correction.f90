PROGRAM bathy_correction
  !!======================================================================
  !!                     ***  PROGRAM  bathy_correction  ***
  !!=====================================================================
  !!  ** Purpose :  eliminate -9999.99 value from input file (in fact <0)
  !!
  !!  ** Method  : work on a copy of the input file
  !!
  !! History :  1.0  : 02/2023  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4) :: ji,jj, jii
  INTEGER(KIND=4) :: narg, ijarg, iargc
  INTEGER(KIND=4) :: ncid, id, ierr
  INTEGER(KIND=4) :: npiglo, npjglo

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: bathy
  REAL(KIND=4)                              :: b1, b2
  REAL(KIND=4)                              :: pp_spval=0

  LOGICAL           :: lmean = .FALSE.

  CHARACTER(LEN=80) :: cldum
  CHARACTER(LEN=80) :: cf_in
  CHARACTER(LEN=80) :: cf_out
  CHARACTER(LEN=80) :: cv_in = 'Bathymetry'
  


  !!----------------------------------------------------------------------
  narg=iargc()
  cf_out='none'
  IF ( narg == 0 ) THEN
     PRINT *,' usage : bathy_correction -f BATHY-in [-o BATHY-out] [-mean]'
     PRINT *,'       [-b BATHY-name]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'          This program reset negative value to 0 (in general) or to a '
     PRINT *,'          mean value corresponding to nearest points in the E-W direction'
     PRINT *,'           (case of 180 deg)'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f BATHY-in : specify a NEMO_like bathymetric file' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -o BATHY-out : specify a name for output file (default is <BATHY_in>_corrected)  '
     PRINT *,'       -o BATHY-name : specify a name for the bathy variable [',TRIM(cv_in),']'
     PRINT *,'       -mean : replace -9999.99 by a mean value of the nearest points'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : <BATHY_in>_corrected) unless -o option is used'
     PRINT *,'         variables : same as in input file'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'    ) ; CALL getarg(ijarg, cf_in  ) ; ijarg=ijarg+1
        ! option
     CASE ( '-o'    ) ; CALL getarg(ijarg, cf_out ) ; ijarg=ijarg+1
     CASE ( '-b'    ) ; CALL getarg(ijarg, cv_in  ) ; ijarg=ijarg+1
     CASE ( '-mean' ) ; lmean=.TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  IF ( cf_out == 'none' ) THEN
     cf_out=TRIM(cf_in)//'_corrected'
  ENDIF
  ierr = copyfile(cf_in, cf_out) 

  ierr = NF90_OPEN(cf_out, NF90_WRITE, ncid)
  ierr = NF90_INQ_DIMID(ncid,'x', id) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid,'y', id) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)

  ALLOCATE( bathy(npiglo, npjglo) )
  ierr = NF90_INQ_VARID(ncid, cv_in, id)
  ierr = NF90_GET_VAR(ncid,id, bathy) 

  IF ( lmean) THEN
!   print *, ' Not coded yet !'
!   stop
    DO  jj=npjglo, 1, -1
      DO ji=1,npiglo
        IF ( bathy(ji,jj) < 0 ) THEN
          IF ( ji /= 1 ) THEN
            b1=bathy(ji-1,jj )
            DO jii=ji+1,npiglo
              b2=bathy(jii,jj)
              IF ( b2 >= 0 ) THEN
                bathy(ji,jj)=(b1+b2)/2.
                exit
              ENDIF
            ENDDO
          ENDIF
        ENDIF  
     ENDDO        
    END DO
  ELSE
    WHERE (bathy < 0 ) bathy = 0.
  ENDIF
  ierr = NF90_PUT_VAR(ncid, id, bathy)
  ierr = NF90_CLOSE(ncid)

CONTAINS
   INTEGER FUNCTION copyfile (cd_fin, cd_fout)
    CHARACTER(LEN=*) , INTENT(in) :: cd_fin, cd_fout

    CHARACTER(LEN=80) :: cl_cmd
    cl_cmd='dd bs=1GB if='//TRIM(cd_fin)//' of='//TRIM(cd_fout)
    PRINT *, TRIM( cl_cmd)
    CALL SYSTEM ( cl_cmd)
    copyfile=0
   END FUNCTION copyfile
  


  !!----------------------------------------------------------------------
END PROGRAM bathy_correction
