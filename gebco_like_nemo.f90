PROGRAM bedmachine_time
  !!======================================================================
  !!                     ***  PROGRAM   bedmachine_time ***
  !!=====================================================================
  !!  ** Purpose : Compute add time_counter dimension for NEMO compat
  !!
  !!  ** Method  : read compute write ...
  !!
  !! History :  1.0  : 03/2023  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER(KIND=4) :: narg, ijarg, iargc
  INTEGER(KIND=4) :: npiglo, npjglo
  INTEGER(KIND=4) :: ncid, id, ierr, ncidi
  INTEGER(KIND=4) :: nvar    ! sumber of variables to process

  REAL(KIND=4)                              :: spval=-9999.99

  CHARACTER(LEN=80) :: cf_in
  CHARACTER(LEN=80) :: cf_out = 'GEBCO_like_nemo.nc'
  CHARACTER(LEN=80) :: cv_lon = 'nav_lon'
  CHARACTER(LEN=80) :: cv_lat = 'nav_lat'
  CHARACTER(LEN=80) :: cv_x = 'lon'
  CHARACTER(LEN=80) :: cv_y = 'lat'
  CHARACTER(LEN=80) :: cv_time = 'time_counter'
  CHARACTER(LEN=80) :: cv_lst                    ! list of input variables
  CHARACTER(LEN=80) :: cv_prec                   ! list of input variables
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE  :: cv_in , cv_inprec
  CHARACTER(LEN=80) :: cldum
  

  !!----------------------------------------------------------------------
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : gebco_like_nemo.exe -f GEBCO-file [ -o GEBCO_like_nemo]'
     PRINT *,'         -v VAR-lst -prec PREC-lst'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        Create a GEBCO file with the same data than the input file, but'
     PRINT *,'        with 2D nav_lon, nav_lat as well as a time_counter, for copatibility'
     PRINT *,'        with tools for NEMO.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f GEBCO-file : specify the name of the GEBCO input file'
     PRINT *,'       -v VAR-lst : Comma separated list of variables to process in the'
     PRINT *,'             input file.'
     PRINT *,'       -prec PREC-lst : Comma separated list of the precision of the variables'
     PRINT *,'           given in the -v argument. Allowed are I2 I4 R4 R8 B (B=BYTE)'
     PRINT *,'     OPTIONS : '
     PRINT *,'       -o GEBCO-like-nemo : name of output file'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        none'
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cf_out) ,' or name given by -o option'
     PRINT *,'         variables : The same than in input file'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in ) ; ijarg=ijarg+1
     CASE ( '-v'   ) ; CALL getarg(ijarg, cv_lst) ; ijarg=ijarg+1
                     ; CALL ParseVars(cv_lst,cv_in) 
        ! option
     CASE ( '-prec' ) ; CALL getarg(ijarg, cv_prec) ; ijarg=ijarg+1
                      ; CALL ParseVars(cv_prec,cv_inprec) 
     CASE ( '-o'    ) ; CALL getarg(ijarg, cf_out     ) ; ijarg=ijarg+1
     CASE DEFAULT     ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  IF (.NOT. ALLOCATED(cv_inprec)  ) THEN
     ALLOCATE (cv_inprec(nvar) )
     cv_inprec(1:nvar) = 'R4'  ! set default precision
  ENDIF

  CALL WriteOutput( cf_out)

  CONTAINS
   SUBROUTINE ParseVars (cdum, cd_tab)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ParseVars  ***
      !!
      !! ** Purpose :  Decode variables names to be used
      !!
      !! ** Method  :  look for , in the argument string and set the number of
      !!         variable (nvaro), allocate cv_fix array and fill it with the
      !!         decoded  names.
      !!
      !!----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: cdum
      CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(out) :: cd_tab

      CHARACTER(LEN=80), DIMENSION(100) :: cl_dum  ! 100 is arbitrary
      INTEGER  :: ji
      INTEGER  :: inchar,  i1, ivar
      !!----------------------------------------------------------------------

      ivar=1
      i1=1
      inchar= LEN(TRIM(cdum))
      ! scan the input string and look for ',' as separator
      DO ji=1,inchar
         IF ( cdum(ji:ji) == ',' ) THEN
            cl_dum(ivar) = cdum(i1:ji-1)
            i1=ji+1
            ivar=ivar+1
         ENDIF
      ENDDO
      nvar=ivar

      ! last name of the list does not have a ','
      cl_dum(nvar) = cdum(i1:inchar)

      ALLOCATE ( cd_tab(nvar) )
      DO ji=1, nvar
         cd_tab(ji) = cl_dum(ji)
      ENDDO
   END SUBROUTINE ParseVars

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

    INTEGER(KIND=4) :: jvar, jatt, ji, jj
    INTEGER(KIND=4) :: ncid, id, idx,idy,idt, ierr, idlon, idlat
    INTEGER(KIND=4) :: ncidi
    INTEGER(KIND=4) :: idvx, idvy, idvtim, iatt, idmap
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: idv
    INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: iprec
    REAL(KIND=4), DIMENSION(:)  , ALLOCATABLE :: rvlon, rvlat


    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dltab
    REAL(KIND=8), DIMENSION(1)                :: dltime

    CHARACTER(LEN=255) :: clnam
    
    !!----------------------------------------------------------------------
    ierr=NF90_OPEN(cf_in, NF90_NOWRITE, ncidi)
    ierr=NF90_INQ_DIMID( ncidi, 'lon', id) ; ierr = NF90_INQUIRE_DIMENSION(ncidi,id,len=npiglo)
    ierr=NF90_INQ_DIMID( ncidi, 'lat', id) ; ierr = NF90_INQUIRE_DIMENSION(ncidi,id,len=npjglo)

    PRINT *, 'NPIGLO = ', npiglo
    PRINT *, 'NPJGLO = ', npjglo

    ALLOCATE(dltab(npiglo, npjglo) )
    ALLOCATE(rvlon(npiglo) )
    ALLOCATE(rvlat(npjglo) )
    dltime=0.d0

    ALLOCATE( idv(nvar), iprec(nvar) )
    DO jvar = 1, nvar
       SELECT CASE (cv_inprec(jvar)) 
       CASE ( 'I2' ) ; iprec(jvar) = NF90_SHORT
       CASE ( 'I4' ) ; iprec(jvar) = NF90_INT
       CASE ( 'R4' ) ; iprec(jvar) = NF90_FLOAT
       CASE ( 'R8' ) ; iprec(jvar) = NF90_DOUBLE
       CASE ( 'B' )  ; iprec(jvar) = NF90_BYTE
       CASE DEFAULT 
           PRINT *,' Precision ',TRIM(cv_inprec(jvar)),' not supported...'
           STOP
       END SELECT
    ENDDO
    ierr=NF90_CREATE(cd_fout,NF90_NETCDF4, ncid)
    ! Dimensions
    ierr = NF90_DEF_DIM(ncid,'x',npiglo, idx)
    ierr = NF90_DEF_DIM(ncid,'y',npjglo, idy)
    ierr = NF90_DEF_DIM(ncid,'time_counter', NF90_UNLIMITED,idt)
    ! Variables
    ierr = NF90_DEF_VAR(ncid,'crs' , NF90_CHAR,  (/1/) ,idmap )

    ierr = NF90_DEF_VAR(ncid,cv_lon , NF90_FLOAT,  (/idx,idy/) ,idlon, deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,cv_lat , NF90_FLOAT,  (/idx,idy/) ,idlat, deflate_level=1,chunksizes=(/1,npjglo/10/) )
    ierr = NF90_DEF_VAR(ncid,cv_x   , NF90_DOUBLE,    (/idx/)     ,idvx,  deflate_level=1,chunksizes=(/npiglo/)      )
    ierr = NF90_DEF_VAR(ncid,cv_y   , NF90_DOUBLE,    (/idy/)     ,idvy,  deflate_level=1,chunksizes=(/npjglo/)      )
    ierr = NF90_DEF_VAR(ncid,cv_time ,NF90_DOUBLE, (/idt/)     ,idvtim                                            )
    DO jvar =1, nvar
       ierr = NF90_DEF_VAR(ncid,cv_in(jvar),iprec(jvar), (/idx,idy,idt/) ,idv(jvar), deflate_level=1 , chunksizes=(/1,npjglo/10,1/) )
!       ierr = NF90_PUT_ATT(ncid,idv(jvar) ,'missing_value',spval)
    ENDDO
    ! Attribute (copy them from input file )
    !crs
    ierr = NF90_INQ_VARID( ncidi,'crs',id) 
    ierr = NF90_INQUIRE_VARIABLE(ncidi,id,nAtts=iatt)
    DO jatt = 1, iatt
      ierr = NF90_INQ_ATTNAME(ncidi,id,jatt,clnam)
      ierr = NF90_COPY_ATT(ncidi,id,clnam,ncid,idmap)
    ENDDO
    ! lon
    ierr = NF90_INQ_VARID( ncidi,cv_x,id) 
    ierr = NF90_INQUIRE_VARIABLE(ncidi,id,nAtts=iatt)
    DO jatt = 1, iatt
      ierr = NF90_INQ_ATTNAME(ncidi,id,jatt,clnam)
      ierr = NF90_COPY_ATT(ncidi,id,clnam,ncid,idvx)
    ENDDO
    ! lat
    ierr = NF90_INQ_VARID( ncidi,cv_y,id) 
    ierr = NF90_INQUIRE_VARIABLE(ncidi,id,nAtts=iatt)
    DO jatt = 1, iatt
      ierr = NF90_INQ_ATTNAME(ncidi,id,jatt,clnam)
      ierr = NF90_COPY_ATT(ncidi,id,clnam,ncid,idvy)
    ENDDO
    ! Variables to workwith
    DO jvar=1,nvar
       ierr = NF90_INQ_VARID( ncidi,cv_in(jvar),id) 
       ierr = NF90_INQUIRE_VARIABLE(ncidi,id,nAtts=iatt)
       DO jatt = 1, iatt
         ierr = NF90_INQ_ATTNAME(ncidi,id,jatt,clnam)
         ierr = NF90_COPY_ATT(ncidi,id,clnam,ncid,idv(jvar))
       ENDDO
    ENDDO
    ! nav_lon, nav_lat time_counter ...
    ierr = NF90_PUT_ATT( ncid,idlon,'standard_name','longitude' )
    ierr = NF90_PUT_ATT( ncid,idlon,'long_name'    ,'Longitude' )
    ierr = NF90_PUT_ATT( ncid,idlon,'units'        ,'degrees_east')

    ierr = NF90_PUT_ATT( ncid,idlat,'standard_name','latitude' )
    ierr = NF90_PUT_ATT( ncid,idlat,'long_name'    ,'Latitude' )
    ierr = NF90_PUT_ATT( ncid,idlat,'units'        ,'degrees_north')
    
    ierr = NF90_PUT_ATT( ncid,idvtim,'axis'         ,'T'         )
    ierr = NF90_PUT_ATT( ncid,idvtim,'standard_name','time'      )
    ierr = NF90_PUT_ATT( ncid,idvtim,'long_name'    ,'Time_axis' )
    ierr = NF90_PUT_ATT( ncid,idvtim,'calendar'    ,'gregorian'  )
    ierr = NF90_PUT_ATT( ncid,idvtim,'units'       ,'seconds since 2022-01-01 00:00:00'  )
    ierr = NF90_PUT_ATT( ncid,idvtim,'time_origin' ,'2022-01-01 00:00:00'  )

    ! Global attributes
    id=NF90_GLOBAL
    ierr = NF90_INQUIRE_VARIABLE(ncidi,id,nAtts=iatt)
    DO jatt = 1, iatt
      ierr = NF90_INQ_ATTNAME(ncidi,id,jatt,clnam)
      ierr = NF90_COPY_ATT(ncidi,id,clnam,ncid,NF90_GLOBAL)
    ENDDO

    ierr = NF90_ENDDEF(ncid)
    ! now write :
    ierr = NF90_PUT_VAR(ncid,idmap,'')   ! crs
    ierr = NF90_INQ_VARID (ncidi,cv_x, id ) ; ierr = NF90_GET_VAR(ncidi,id,rvlon(:) ) ; ierr = NF90_PUT_VAR(ncid,idvx,rvlon)  ! lon
    DO jj = 1, npjglo
       dltab(:,jj) = rvlon(:)
    ENDDO
    ierr = NF90_PUT_VAR(ncid,idlon,dltab)  ! nav_lon
    ierr = NF90_INQ_VARID (ncidi,cv_y, id ) ; ierr = NF90_GET_VAR(ncidi,id,rvlat(:) ) ; ierr = NF90_PUT_VAR(ncid,idvy,rvlat)  ! lat
    DO ji = 1, npiglo
       dltab(ji,:) = rvlat(:)
    ENDDO
    ierr = NF90_PUT_VAR(ncid,idlat,dltab)  ! nav_lat

    ierr = NF90_PUT_VAR(ncid,idvtim,dltime, start=(/1/), count=(/1/) )      ! time counter

    DO jvar = 1, nvar
      PRINT *, 'working for ', TRIM(cv_in(jvar))
      ierr = NF90_INQ_VARID (ncidi,cv_in(jvar), id ) 
      ierr = NF90_GET_VAR(ncidi,id,dltab(:,:) ) 
      ierr = NF90_PUT_VAR(ncid,idv(jvar),dltab(:,:), start=(/1,1,1/), count=(/npiglo,npjglo,1/))  ! 
    ENDDO

    ierr = NF90_CLOSE(ncid)
    ierr = NF90_CLOSE(ncidi)

    END SUBROUTINE WriteOutput

END PROGRAM bedmachine_time
