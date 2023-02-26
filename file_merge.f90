PROGRAM file_merge_AA
  !!======================================================================
  !!                     ***  PROGRAM  <module>  ***
  !!=====================================================================
  !!  ** Purpose : One shot program to reconstruct splitted files for
  !!               NEMOBAT in AA (not standard splitting )
  !!
  !!  ** Method  :  NF90 primitive
  !!
  !! History :  1.0  : 02/2023  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4) :: ji, jj, jf
  INTEGER(KIND=4) :: i1,i2,ig1,ig2
  INTEGER(KIND=4) :: narg, ijarg, iargc
  INTEGER(KIND=4) :: ncid, id, ierr, ncoo
  INTEGER(KIND=4) :: idx, idy, idvlon, idvlat, idv
  INTEGER(KIND=4) :: nfile
  INTEGER(KIND=4) :: npi, npj, npiglo, npjglo
  INTEGER(KIND=4),DIMENSION(:,:), ALLOCATABLE :: ilimit

  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rglo, rwk
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: rloc

  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: cf_in
  CHARACTER(LEN=255)                           :: cf_list
  CHARACTER(LEN=255)                           :: cl_lim
  CHARACTER(LEN=255)                           :: cf_glob='none'
  CHARACTER(LEN=80)                            :: cf_out='none'
  CHARACTER(LEN=80)                            :: cv_in='none'
  CHARACTER(LEN=80)                            :: clcut
  CHARACTER(LEN=255)                           :: cldum

  LOGICAL                                      :: lneg=.FALSE. ! flag for changing sign

  !!----------------------------------------------------------------------
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : file_merge_AA.exe -l FILE-list -o FILE-out'
     PRINT *,'    -ilim LIMIT-list -g GLOBAL-file -v VAR-name -neg'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Reconstruct the global file from the splitted files '
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l FILE-list : a comma separated list of  files to merge' 
     PRINT *,'       -ilim LIMIT-list : a comma separated list of  limit pairs'
     PRINT *,'            describing the position of FILES in global domain'
     PRINT *,'           Example -imim 1,730,700,1430,1400,2143 '
     PRINT *,'             In this example : 3 pairs (corresponding to 3 files)'
     PRINT *,'       -g GLOBAL-file :name of global coordinate file'
     PRINT *,'       -v VAR-name :name of variable to recombine'
     PRINT *,'       -o FILE-out : name of output file'
     PRINT *,'       -neg : change sign of var in the output file'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       ' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) 
     PRINT *,'         variables : same as input files'
     PRINT *,'      ' 
     PRINT *,'     SEE ALSO :'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-l'   ) ; CALL getarg(ijarg, cf_list ) ; ijarg=ijarg+1
        CALL ParseFile(cf_list)
        ALLOCATE (ilimit(nfile,2) )
     CASE ( '-ilim') ; CALL getarg(ijarg, cl_lim  ) ; ijarg=ijarg+1
        CALL ParseLim(cl_lim)
     CASE ( '-g'   ) ; CALL getarg(ijarg, cf_glob ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cv_in   ) ; ijarg=ijarg+1
     CASE ( '-neg' ) ; lneg=.true.
        ! option
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ierr = 0
  IF ( cf_glob == 'none') THEN
     ierr = ierr +1
     PRINT *, "  Need to give the name of global coordinates with -g"
  ELSE
     PRINT *, ' Global coordinates file : ', TRIM(cf_glob)
  ENDIF
  IF ( cf_out == 'none') THEN
     ierr = ierr +1
     PRINT *, "  Need to give the name of output file with -o"
  ELSE
     PRINT *, ' Output file : ', TRIM(cf_out)
  ENDIF
  IF ( cv_in == 'none') THEN
     ierr = ierr +1
     PRINT *, "  Need to give the name of variable  with -v"
  ELSE
     PRINT *, ' Working variables : ', TRIM(cv_in)
  ENDIF

  IF  ( ierr /= 0 ) THEN
     PRINT *, ierr,' error(s) detected --> stop !'
     STOP
  ENDIF

  PRINT *,nfile,'  files to merge:'
  DO ji =1, nfile
     PRINT *, TRIM(cf_in(ji) ), ilimit(ji,1),' - ',ilimit(ji,2)
  ENDDO

  ! look for global domain size
  ierr=NF90_OPEN(cf_glob,NF90_NOWRITE,ncid)
  ierr=NF90_INQ_DIMID(ncid,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr=NF90_INQ_DIMID(ncid,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)
  PRINT *, ' GLOBAL SIZE : ', npiglo,' x ',npjglo
  ierr=NF90_CLOSE(ncid)
  ALLOCATE ( rglo(npiglo,npjglo) )

  DO jf=1,nfile
     ierr=NF90_OPEN(cf_in(jf),NF90_NOWRITE,ncid)
     ierr=NF90_INQ_DIMID(ncid,'x',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npi)
     ierr=NF90_INQ_DIMID(ncid,'y',id) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npj)
     IF (jf == 1 ) THEN
        ierr=0
        IF ( npi == npiglo ) THEN
           clcut='V'
           IF ( ilimit (1,1) /= 1 ) THEN
              ierr=ierr+1
              PRINT *,' ilimit (1,1) != 1 !'
           ENDIF
           IF ( ilimit (nfile,2) /= npjglo ) THEN
              ierr=ierr+1
              PRINT *,' ilimit (nfile,2) != npjglo !'
           ENDIF
        ELSEIF ( npj == npjglo ) THEN
           clcut='H'
           IF ( ilimit (1,1) /= 1 ) THEN
              ierr=ierr+1
              PRINT *,' ilimit (1,1) != 1 !'
           ENDIF
           IF ( ilimit (nfile,2) /= npiglo ) THEN
              ierr=ierr+1
              PRINT *,' ilimit (nfile,2) != npiglo !'
           ENDIF

        ELSE
           PRINT *,' cannot infer cutting direction ...PROBLEM!'
           STOP
        ENDIF
        PRINT *,' Cut direction :', TRIM (clcut)
        IF (ierr /= 0 ) STOP 
     ENDIF
     ALLOCATE( rloc(npi,npj) )
     ierr= NF90_INQ_VARID(ncid,cv_in,id) ; ierr=NF90_GET_VAR(ncid,id,rloc)
     ierr=NF90_CLOSE(ncid)
     IF (clcut == 'H') THEN
        IF ( jf == 1     ) ig1=1
        IF ( jf == nfile ) THEN 
          ig2=npiglo
        ELSE
          ig2=(ilimit(jf+1,1) + ilimit(jf,2) ) / 2
        ENDIF

        i1=(ig1-ilimit(jf,1)) +1
        i2 = i1+( ig2 - ig1  )
        
        PRINT *, jf, ig1, ig2, i1,i2, npi, npj

        rglo(ig1:ig2,:)=rloc(i1:i2,:)
        ig1=ig2+1
     ELSE ! clcut='V'
        IF ( jf == 1     ) ig1=1
        IF ( jf == nfile ) THEN
          ig2=npjglo
        ELSE
          ig2=(ilimit(jf+1,1) + ilimit(jf,2) ) / 2
        ENDIF

        i1=(ig1-ilimit(jf,1)) +1
        i2 = i1+( ig2 - ig1  )

        PRINT *, jf, ig1, ig2, i1,i2, npi, npj

        rglo(:,ig1:ig2)=rloc(:,i1:i2)
        ig1=ig2+1
     ENDIF


     DEALLOCATE(rloc)
  ENDDO
  IF  (lneg ) rglo=-1.0*rglo
  ! violent ... assume positive variable (keep -9999.99 flag value )
  WHERE ( rglo < 0 .AND. rglo > -8000 ) rglo = 0
  IF ( clcut == 'H' ) THEN
     ! case of AA, reset periodic condition:
     rglo(1     ,:) = rglo(npiglo-1,:)
     rglo(npiglo,:) = rglo(2       ,:)
  ENDIF
  ! Write global file (taking nav_lon nav_lat from global coord)
  ierr= NF90_CREATE(cf_out,NF90_NETCDF4,ncid)
  ierr= NF90_DEF_DIM( ncid,'x',npiglo,idx)
  ierr= NF90_DEF_DIM( ncid,'y',npjglo,idy)

  ierr=NF90_DEF_VAR( ncid,'nav_lon',NF90_FLOAT,(/idx,idy/), idvlon)
  ierr=NF90_DEF_VAR( ncid,'nav_lat',NF90_FLOAT,(/idx,idy/), idvlat)
  ierr=NF90_DEF_VAR( ncid,cv_in,NF90_FLOAT,(/idx,idy/), idv)

  ierr=NF90_ENDDEF(ncid)

! retrieve nav_lon, nav_lat
  ALLOCATE( rwk(npiglo,npjglo) )
  ierr = NF90_OPEN(cf_glob,NF90_NOWRITE,ncoo)
  ierr = NF90_INQ_VARID(ncoo,'glamt',id) ; ierr=NF90_GET_VAR(ncoo,id,rwk)
  ierr=NF90_PUT_VAR( ncid, idvlon, rwk)
  ierr = NF90_INQ_VARID(ncoo,'gphit',id) ; ierr=NF90_GET_VAR(ncoo,id,rwk)
  ierr = NF90_CLOSE(ncoo)
  ierr=NF90_PUT_VAR( ncid, idvlat, rwk)
  DEALLOCATE (rwk)

  ierr=NF90_PUT_VAR( ncid, idv, rglo)
  ierr=NF90_CLOSE(ncid)


CONTAINS
  SUBROUTINE ParseFile (cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseFile  ***
    !!
    !! ** Purpose :  Decode -l option
    !!
    !! ** Method  :  look for , in the argument string and set the number of
    !!         variable (nvaro), allocate cv_fix array and fill it with the
    !!         decoded  names.
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdum

    CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1
    !!----------------------------------------------------------------------

    inchar= LEN(TRIM(cdum))
    ! scan the input string and look for ',' as separator
    nfile=1
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cl_dum(nfile) = cdum(i1:ji-1)
          i1=ji+1
          nfile=nfile+1
       ENDIF
    ENDDO

    ! last name of the list does not have a ','
    cl_dum(nfile) = cdum(i1:inchar)

    ALLOCATE ( cf_in(nfile) )
    DO ji=1, nfile
       cf_in(ji) = cl_dum(ji)
    ENDDO
  END SUBROUTINE ParseFile

  SUBROUTINE ParseLim (cdum)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ParseLim  ***
    !!
    !! ** Purpose :  Decode -ilim option'
    !!
    !! ** Method  :  look for , in the argument string and set the number of
    !!         variable (nvaro), allocate cv_fix array and fill it with the
    !!         decoded  names.
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cdum

    CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
    INTEGER  :: ji
    INTEGER  :: inchar,  i1=1, imax
    !!----------------------------------------------------------------------

    inchar= LEN(TRIM(cdum))
    ! scan the input string and look for ',' as separator
    imax=1
    DO ji=1,inchar
       IF ( cdum(ji:ji) == ',' ) THEN
          cl_dum(imax) = cdum(i1:ji-1)
          i1=ji+1
          imax=imax+1
       ENDIF
    ENDDO

    ! last name of the list does not have a ','
    cl_dum(imax) = cdum(i1:inchar)
    IF ( imax /= nfile * 2 ) THEN
       PRINT *,' Uncorrect number of files or ilimits...'
       STOP
    ENDIF

    imax=1
    DO ji=1, nfile
       READ(cl_dum(imax),*) ilimit(ji,1) ; imax=imax+1
       READ(cl_dum(imax),*) ilimit(ji,2) ; imax=imax+1
    ENDDO
  END SUBROUTINE ParseLim

END PROGRAM file_merge_AA

