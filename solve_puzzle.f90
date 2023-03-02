PROGRAM solve_puzzle
  !!======================================================================
  !!                     ***  PROGRAM  solve_puzzle  ***
  !!=====================================================================
  !!  ** Purpose : This program is dedicated at merging different files
  !!             corresponding of pieces of bathymetry, computed separately
  !!             on differents NEMO grid.  The important point is that all grid
  !!             share common points.
  !!             The final target grid is eORCA36, NEMO.4.0 style (whith halos)
  !!             The pieces to assemble :
  !!                * corrected bathymetry on NEMO 4.2 style grid (GEBCO_2022)
  !!                * POLE36 bathy corresponding to the NPole area with halo (4.0)
  !!                * FOLD36 bathy corresponding to the E-W periodic halo (4.0)
  !!                * AUSTRAL36 corresponding to bedmachine-antarctica 3 original data (4.0)
  !!                * GREENLAND36 corresponding to bedmachine-greenland 5 original data (4.0)
  !!              All subdomains have their own coordinate file.
  !!
  !!  ** Method  : Create a bathy variable at the size of the target grid
  !!               Fill the bathy with different sources, checking the grid matching
  !!               Do to the unique usage of this program file names are hard coded...
  !!
  !! History :  1.0  : 02/2023  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4) :: ji42, jj42, jigre, jjgre
  INTEGER(KIND=4) :: ncid, id, ierr

  INTEGER(KIND=4) :: npi40, npj40,  npi_pol,npj_pol,   npi_fol,npj_fol,   npi_gre,npj_gre
  INTEGER(KIND=4) :: npi42, npj42,  npi_aus,npj_aus
  INTEGER(KIND=4) :: npi,npj
  INTEGER(KIND=4) :: ii40, ij40, iifol,ijfol,  iiloc, ijloc,  iiaus,ijaus,   iigre, ijgre
  INTEGER(KIND=4) :: ni40, nj40
  

  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rbat40, rdra40
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rbat42, rbat_pol, rbat_fol, rbat_gre, rbat_aus
  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: rma_gre, rdra_gre
  REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: rco40, rco42, rco_pol, rco_fol, rco_gre, rco_aus

  ! coordinates files
  CHARACTER(LEN=255) :: cf_coor_40     = 'eORCA_R36_coordinates_v3.3.nc4'
  CHARACTER(LEN=255) :: cf_coor_42     = 'eORCA_R36_coordinates_v3.0.noz.nc4'
  CHARACTER(LEN=255) :: cf_coor_40_pol = 'POLE36_coordinates_3.3.nc'
  CHARACTER(LEN=255) :: cf_coor_40_fol = 'FOLD36_coordinates_3.3.nc'
  CHARACTER(LEN=255) :: cf_coor_40_aus = 'AUSTRAL36_coordinates_v3.3.nc'
  CHARACTER(LEN=255) :: cf_coor_40_gre = 'GREENLAND36_coordinates_v3.3.nc'
  ! Bathymetric files
  CHARACTER(LEN=255) :: cf_bat_40     = 'eORCA_R36_Bathymetry_v4.0.2.nc'
  CHARACTER(LEN=255) :: cf_bat_42     = 'eORCA36_GEBCO_2022_PM-JMM_bathy_v2_4.2.nc'
  CHARACTER(LEN=255) :: cf_bat_40_pol = 'POLE36_GEBCO2022_Bathymetry_time_corrected.nc'
  CHARACTER(LEN=255) :: cf_bat_40_fol = 'FOLD36_GEBCO2022_Bathymetry_clean.nc'
  CHARACTER(LEN=255) :: cf_bat_40_aus = 'AUSTRAL36_bedmachine3_Bathymetry.nc'
  CHARACTER(LEN=255) :: cf_bat_40_gre = 'GREENLAND36_bedmachine5_Bathymetry.nc'
  ! Icedraft
  CHARACTER(LEN=255) :: cf_dra_40_aus = 'AUSTRAL36_bedmachine3_Draft.nc'
  CHARACTER(LEN=255) :: cf_dra_40_gre = 'GREENLAND36_bedmachine5_Draft.nc'
  ! Mask
  CHARACTER(LEN=255) :: cf_msk_40_gre = 'GREENLAND36_bedmachine5_Mask2.nc'

  LOGICAL   :: lflag
  !!----------------------------------------------------------------------
  ! read coord for each grid
  CALL GetCoor(cf_coor_40    , rco40  , npi40, npj40    ) ; PRINT *,TRIM(cf_coor_40)//' done.'
  CALL GetCoor(cf_coor_42    , rco42  , npi42, npj42    ) ; PRINT *,TRIM(cf_coor_42)//' done.'
  CALL GetCoor(cf_coor_40_pol, rco_pol, npi_pol, npj_pol) ; PRINT *,TRIM(cf_coor_40_pol)//' done.'
  CALL GetCoor(cf_coor_40_fol, rco_fol, npi_fol, npj_fol) ; PRINT *,TRIM(cf_coor_40_fol)//' done.'
  CALL GetCoor(cf_coor_40_gre, rco_gre, npi_gre, npj_gre) ; PRINT *,TRIM(cf_coor_40_gre)//' done.'
  CALL GetCoor(cf_coor_40_aus, rco_aus, npi_aus, npj_aus) ; PRINT *,TRIM(cf_coor_40_aus)//' done.'

  ! Deal with 4.2.0 grid
  ! look for a common point with target grid (4.0)
  iiloc=5600; ijloc=5600 
  CALL NearestPoint (rco42(iiloc,ijloc,1), rco42(iiloc,ijloc,2), npi40,npj40, rco40(:,:,1), rco40(:,:,2), &
                     ii40,ij40, lflag,.true.)
  PRINT *, iiloc, ijloc, ii40, ij40, lflag
  PRINT *, rco42(iiloc,ijloc,1), rco40(ii40,ij40,1)
  PRINT *, rco42(iiloc,ijloc,2), rco40(ii40,ij40,2)

 ALLOCATE( rbat40(npi40,npj40) )
 CALL GetBat(cf_bat_42,'Bathymetry',rbat42,npi,npj)
 IF ( npi /= npi42 .OR. npj /= npj42 ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_bat_42)
 ENDIF
 ! compute limit of patch
 DO ji42=1,npi42
    ni40=ji42-iiloc+ii40 
   DO jj42=1,npj42
      nj40=jj42-ijloc+ij40 
      IF ( ni40  >= 1 .AND. ni40 <= npi40 ) THEN
         IF ( nj40  >= 1 .AND. nj40 <= npj40 ) THEN
            rbat40(ni40,nj40)=rbat42(ji42,jj42)
         ENDIF
      ENDIF
   ENDDO
 ENDDO
 DEALLOCATE (rbat42)

 ! Now deal with periodic boundaries
 iiloc=1; ijloc=5800
 CALL NearestPoint (rco40(iiloc,ijloc,1), rco40(iiloc,ijloc,2), npi_fol,npj_fol, rco_fol(:,:,1), rco_fol(:,:,2), &
                     iifol,ijfol, lflag,.true.)
  PRINT *, iiloc, ijloc, iifol, ijfol, lflag
  PRINT *, rco40(iiloc,ijloc,1), rco_fol(iifol,ijfol,1)
  PRINT *, rco40(iiloc,ijloc,2), rco_fol(iifol,ijfol,2)

 CALL GetBat(cf_bat_40_fol,'Bathymetry',rbat_fol,npi,npj)
 IF ( npi /= npi_fol .OR. npj /= npj_fol ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_bat_40_fol)
 ENDIF
  rbat40(1:2,:) = rbat_fol(iifol:iifol+1,:)
  rbat40(npi40-1:npi40,:)=rbat40(1:2,:)
  DEALLOCATE(rbat_fol)

 ! now deal with the POLE
 CALL GetBat(cf_bat_40_pol,'Bathymetry',rbat_pol,npi,npj)
 IF ( npi /= npi_pol .OR. npj /= npj_pol ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_bat_40_pol)
 ENDIF

 rbat40(:,npj40-5:npj40) = rbat_pol(:,npj_pol-5:npj_pol)
  DEALLOCATE(rbat_pol)

! Now patch the Antarctica j=1 --> j=3336
 iiloc=6000; ijloc=3336
 ! check coherency :
 IF ( rco40(iiloc,ijloc,1) /= rco_aus(iiloc,ijloc,1) ) THEN
    PRINT *, 'PB in AUSTRAL GRID I'
 ENDIF
 IF ( rco40(iiloc,ijloc,2) /= rco_aus(iiloc,ijloc,2) ) THEN
    PRINT *, 'PB in AUSTRAL GRID J'
 ENDIF

 ALLOCATE( rdra40(npi40,npj40) )

 CALL GetBat(cf_bat_40_aus,'Bathymetry',rbat_aus,npi,npj)
 IF ( npi /= npi_aus .OR. npj /= npj_aus ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_bat_40_aus)
 ENDIF
 rbat40(:,1:ijloc)=rbat_aus(:,1:ijloc)

 CALL GetBat(cf_dra_40_aus,'icedraft',rbat_aus,npi,npj)
 IF ( npi /= npi_aus .OR. npj /= npj_aus ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_dra_40_aus)
 ENDIF
 rdra40(:,:) = 0.0  ! icedraft for aus (and gre later):0
 rdra40(:,1:ijloc)=rbat_aus(:,1:ijloc)
 DEALLOCATE(rbat_aus)

! now patch Greenland Bathymetry and icedraft
 iiloc=701 ; ijloc=1211
 ! CALL NearestPoint (rco_gre(iiloc,ijloc,1), rco_gre(iiloc,ijloc,2), npi40,npj40, rco40(:,:,1), rco40(:,:,2), &
 !                    ii40,ij40, lflag,.true.)
 ! For some reason the above NearestPoint search fail to find the common point. Set it manually
 ii40=9100  ; ij40=9810
  PRINT *, iiloc, ijloc, ii40, ij40, lflag
  PRINT *, rco_gre(iiloc,ijloc,1), rco40(ii40,ij40,1)
  PRINT *, rco_gre(iiloc,ijloc,2), rco40(ii40,ij40,2)

 CALL GetBat(cf_msk_40_gre,'tmask',rma_gre,npi,npj)
 IF ( npi /= npi_gre .OR. npj /= npj_gre ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_msk_40_gre)
 ENDIF
 CALL GetBat(cf_bat_40_gre,'Bathymetry',rbat_gre,npi,npj)
 IF ( npi /= npi_gre .OR. npj /= npj_gre ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_bat_40_gre)
 ENDIF
 CALL GetBat(cf_dra_40_gre,'icedraft',rdra_gre,npi,npj)
 IF ( npi /= npi_gre .OR. npj /= npj_gre ) THEN
   PRINT *, 'dimension error in ',TRIM(cf_dra_40_gre)
 ENDIF
 DO jigre=2,npi_gre-1
    ni40=jigre-iiloc+ii40
   DO jjgre=2,npj_gre-1
      nj40=jjgre-ijloc+ij40
      IF ( ni40  >= 1 .AND. ni40 <= npi40 ) THEN
         IF ( nj40  >= 1 .AND. nj40 <= npj40 ) THEN
            IF ( rma_gre (jigre,jjgre) >= 0 ) THEN
            rbat40(ni40,nj40)=rbat_gre(jigre,jjgre)
            rdra40(ni40,nj40)=rdra_gre(jigre,jjgre)
            ENDIF
         ENDIF
      ENDIF
   ENDDO
 ENDDO
 DEALLOCATE( rbat_gre, rdra_gre )

!! DONE :)
 CALL WriteOutput(cf_bat_40)
    
  
  
CONTAINS
 SUBROUTINE GetCoor(cd_fcoor, pco , kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetCoor  ***
    !!
    !! ** Purpose :   Open coordinates file, set size and read glamt, gphit
    !!
    !! ** Method  :   
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                            INTENT(in   ) :: cd_fcoor
    REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE, INTENT(  out) :: pco
    INTEGER(KIND=4),                             INTENT(  out) :: kpi, kpj

    INTEGER(KIND=4) :: icid, id, ierr
    !!----------------------------------------------------------------------
    ierr = NF90_OPEN(cd_fcoor, NF90_NOWRITE, icid)
    ierr = NF90_INQ_DIMID( icid, 'x', id) ; ierr = NF90_INQUIRE_DIMENSION( icid, id, len=kpi)
    ierr = NF90_INQ_DIMID( icid, 'y', id) ; ierr = NF90_INQUIRE_DIMENSION( icid, id, len=kpj)
    ALLOCATE( pco(kpi, kpj,2) )
    ierr = NF90_INQ_VARID( icid,'glamt', id) ; ierr = NF90_GET_VAR(icid, id, pco(:,:,1 ) )
    ierr = NF90_INQ_VARID( icid,'gphit', id) ; ierr = NF90_GET_VAR(icid, id, pco(:,:,2 ) )
    ierr = NF90_CLOSE(icid)
 END SUBROUTINE GetCoor

 SUBROUTINE GetBat(cd_fbat, cd_vbat, pba , kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetCoor  ***
    !!
    !! ** Purpose :   Open coordinates file, set size and read glamt, gphit
    !!
    !! ** Method  :   
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),                            INTENT(in   ) :: cd_fbat
    CHARACTER(LEN=*),                            INTENT(in   ) :: cd_vbat
    REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE,   INTENT(  out) :: pba
    INTEGER(KIND=4),                             INTENT(  out) :: kpi, kpj

    INTEGER(KIND=4) :: icid, id, ierr
    !!----------------------------------------------------------------------
    ierr = NF90_OPEN(cd_fbat, NF90_NOWRITE, icid)
    ierr = NF90_INQ_DIMID( icid, 'x', id) ; ierr = NF90_INQUIRE_DIMENSION( icid, id, len=kpi)
    ierr = NF90_INQ_DIMID( icid, 'y', id) ; ierr = NF90_INQUIRE_DIMENSION( icid, id, len=kpj)
    ALLOCATE( pba(kpi, kpj) )
    ierr = NF90_INQ_VARID( icid,cd_vbat, id) ; ierr = NF90_GET_VAR(icid, id, pba(:,: ) )
    ierr = NF90_CLOSE(icid)
 END SUBROUTINE GetBat

  SUBROUTINE NearestPoint(ddlon, ddlat, kpi, kpj, ddlam, ddphi, kpiloc, kpjloc, ld_bnd,ld_reset)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE NearestPoint  ***
    !!
    !! ** Purpose : Compute the positions of the nearest i,j in the grid
    !!              from the given longitudes and latitudes
    !!
    !! ** Method  : Starts on the middle of the grid, search in a 20x20 box,
    !!              and move the box in the direction where the distance 
    !!              between the box and the point is minimum.
    !!              Iterates ...
    !!              Stops when the point is outside the grid.
    !!
    !! References : P.A. Darbon and A. de Miranda acknowledged for this 
    !!              clever algorithm developped in CLIPPER.
    !!----------------------------------------------------------------------
    REAL(KIND=4),                     INTENT(in) :: ddlon, ddlat      !: lon and lat of target point
    INTEGER(KIND=4),                 INTENT (in) :: kpi, kpj          !: grid size
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ddlam, ddphi      !: model grid layout
    INTEGER(KIND=4),              INTENT (inout) :: kpiloc, kpjloc    !: nearest point location
    LOGICAL                                      :: ld_bnd            !: reach boundary flag
    LOGICAL                                      :: ld_reset          !: reach boundary flag

    INTEGER(KIND=4)                              :: ji, jj
    INTEGER(KIND=4), PARAMETER                   :: jp_blk=10
    INTEGER(KIND=4)                              :: ii0, ij0
    INTEGER(KIND=4)                              :: ii1, ij1
    REAL(KIND=4)                                 :: zdist
    REAL(KIND=4)                                 :: zdistmin, zdistmin0
    LOGICAL, SAVE                                :: ll_bndcell, ll_first=.TRUE.
    !!----------------------------------------------------------------------
    IF ( ld_reset ) ll_first = .TRUE.
    IF ( ll_first ) THEN
       kpiloc = kpi/2 ; kpjloc = kpj/2    ! seek from the middle of domain
       ll_first=.FALSE.
       PRINT *, 'SEEK point : ', kpiloc, kpjloc
    ENDIF

    zdistmin=1000000. ; zdistmin0=1000000.
    ii0 = kpiloc      ; ij0 = kpjloc
    ll_bndcell=.TRUE. ; ld_bnd=.FALSE.

    ! loop until found or boundary reach
    DO  WHILE ( ll_bndcell .AND. .NOT. ld_bnd )
       ii0 = kpiloc - jp_blk ;  ii1 = kpiloc + jp_blk
       ij0 = kpjloc - jp_blk ;  ij1 = kpjloc + jp_blk

       ! search only the inner domain
       IF (ii0 <= 0 ) ii0 = 2
       IF (ii1 > kpi) ii1 = kpi - 1
       IF (ij0 <= 0 ) ij0 = 2
       IF( ij1 > kpj) ij1 = kpj - 1

       ! within a block jp_blk+1 x jp_blk+1:
       DO jj=ij0,ij1
          DO ji=ii0,ii1
             ! compute true distance (orthodromy) between target point and grid point
             zdist    = dist(ddlon, ddlam(ji,jj), ddlat, ddphi(ji,jj) )
             zdistmin = MIN(zdistmin, zdist)
             ! update kpiloc, kpjloc if distance decreases
             IF (zdistmin /=  zdistmin0 ) THEN
                kpiloc=ji
                kpjloc=jj
             ENDIF
             zdistmin0=zdistmin
          END DO
       END DO

       ll_bndcell=.FALSE.
       ! if kpiloc, kpjloc belong to block boundary proceed to next block, centered on kpiloc, kpjloc
       IF (kpiloc == ii0 .OR. kpiloc == ii1) ll_bndcell=.TRUE.
       IF (kpjloc == ij0 .OR. kpjloc == ij1) ll_bndcell=.TRUE.

       ! boundary reach ---> not found
       IF (kpiloc == 2  .OR. kpiloc ==kpi-1) ld_bnd=.TRUE.
       IF (kpjloc == 2  .OR. kpjloc ==kpj-1) ld_bnd=.TRUE.
    END DO

  END SUBROUTINE  NearestPoint

  REAL(KIND=4) FUNCTION dist(ddlona, ddlonb, ddlata, ddlatb)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION dist  ***
    !!
    !! ** Purpose : Compute the distance (km) between
    !!              point A (lona, lata) and B (lonb, latb)  
    !!
    !! ** Method  : Use of double precision is important. Compute the 
    !!              distance along the orthodromy
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), INTENT(in) :: ddlata, ddlona, ddlatb, ddlonb

    REAL(KIND=4), SAVE :: dl_latar, dl_latbr, dl_lonar, dl_lonbr
    REAL(KIND=4)       :: dl_pds
    REAL(KIND=4), SAVE :: dl_ux, dl_uy, dl_uz
    REAL(KIND=4)       :: dl_vx, dl_vy, dl_vz
    REAL(KIND=4), SAVE :: dl_prevlat=-1000.d0
    REAL(KIND=4), SAVE :: dl_prevlon=-1000.d0
    REAL(KIND=4), SAVE :: dl_r, dl_pi, dl_conv

    LOGICAL :: ll_first=.TRUE.
    !!----------------------------------------------------------------------
    ! initialise some values at first call
    IF ( ll_first ) THEN
       ll_first = .FALSE.
       ! constants
       dl_pi   = ACOS(-1.d0)
       dl_conv = dl_pi/180.d0  ! for degree to radian conversion
       ! Earth radius
       dl_r    = (6378.137d0+6356.7523d0)/2.0d0 ! km
    ENDIF

    ! compute these term only if they differ from previous call
    IF ( ddlata /= dl_prevlat .OR. ddlona /= dl_prevlon) THEN
       dl_latar   = ddlata*dl_conv
       dl_lonar   = ddlona*dl_conv
       dl_ux      = COS(dl_lonar)*COS(dl_latar)
       dl_uy      = SIN(dl_lonar)*COS(dl_latar)
       dl_uz      = SIN(dl_latar)
       dl_prevlat = ddlata
       dl_prevlon = ddlona
    ENDIF

    dl_latbr = ddlatb*dl_conv
    dl_lonbr = ddlonb*dl_conv
    dl_vx    = COS(dl_lonbr)*COS(dl_latbr)
    dl_vy    = SIN(dl_lonbr)*COS(dl_latbr)
    dl_vz    = SIN(dl_latbr)

    dl_pds   = dl_ux*dl_vx + dl_uy*dl_vy + dl_uz*dl_vz

    IF (dl_pds >= 1.) THEN
       dist = 0.
    ELSE
       dist = dl_r*ACOS(dl_pds)
    ENDIF

  END FUNCTION dist

  SUBROUTINE WriteOutput( cd_fout) 
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE WriteOutput ***
    !!
    !! ** Purpose :  Create final bathymetry file with Bathymetry and icedraft 
    !!
    !! ** Method  :  NF90 standard 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_fout

    INTEGER(KIND=4) :: icid, idx,idy, idvlon, idvlat, idvbat, idvdra,ierr
    !!----------------------------------------------------------------------
    ierr = NF90_CREATE( cd_fout, NF90_NETCDF4, icid)
    ! dimension
    ierr = NF90_DEF_DIM( icid,'x',npi40, idx)
    ierr = NF90_DEF_DIM( icid,'y',npj40, idy)
    ! variables
    ierr = NF90_DEF_VAR (icid,'nav_lon',NF90_FLOAT,(/idx,idy/), idvlon, &
       &                 chunksizes=(/npi40/10,npj40/10/), deflate_level=1  )
    ierr = NF90_DEF_VAR (icid,'nav_lat',NF90_FLOAT,(/idx,idy/), idvlat, &
       &                 chunksizes=(/npi40/10,npj40/10/), deflate_level=1  )
    ierr = NF90_DEF_VAR (icid,'Bathymetry',NF90_FLOAT,(/idx,idy/), idvbat, &
       &                 chunksizes=(/npi40/10,npj40/10/), deflate_level=1  )
    ierr = NF90_DEF_VAR (icid,'icedraft',NF90_FLOAT,(/idx,idy/), idvdra, &
       &                 chunksizes=(/npi40/10,npj40/10/), deflate_level=1  )
    ! Attributes
    ierr = NF90_PUT_ATT(icid,idvlon,'standard_name','longitude')
    ierr = NF90_PUT_ATT(icid,idvlon,'long_name','Longitude')
    ierr = NF90_PUT_ATT(icid,idvlon,'units','degrees_east')

    ierr = NF90_PUT_ATT(icid,idvlat,'standard_name','latitude')
    ierr = NF90_PUT_ATT(icid,idvlat,'long_name','Latitude')
    ierr = NF90_PUT_ATT(icid,idvlat,'units','degrees_north')

    ierr = NF90_PUT_ATT(icid,idvbat,'standard_name','bathymetry')
    ierr = NF90_PUT_ATT(icid,idvbat,'long_name','Ocean Bathymetry')
    ierr = NF90_PUT_ATT(icid,idvbat,'units','meters')
    ierr = NF90_PUT_ATT(icid,idvbat,'missing_value',0.)

    ierr = NF90_PUT_ATT(icid,idvdra,'standard_name','icedraft')
    ierr = NF90_PUT_ATT(icid,idvdra,'long_name','Ice Draft')
    ierr = NF90_PUT_ATT(icid,idvdra,'units','meters')

    ierr = NF90_PUT_ATT(icid, NF90_GLOBAL, 'NEMO_compatibility','Suitable for NEMO version <= 4.0.x. Resizing required for version >= 4.2.x' )
    ierr = NF90_PUT_ATT(icid, NF90_GLOBAL, 'Authors',' J.-M. Molines, P. Mathiot, IGE - Grenoble Feb. 2023' )
    ierr = NF90_PUT_ATT(icid, NF90_GLOBAL, 'Description', &
     & 'eORCA36 Bathymetry and icedraft version 4.0. This Bathymetry has been set up using '//&
     & 'GEBCO_2022 database for the global ocean, BedMachine v3 for Antarctica and '// &
     & 'BedMachine-Greenland v5 for Greenland.  Coast line correction (hand made) was performed '// &
     & 'using BMGTOOLS. Details of the making of are available in https::/ add github reference' )
    ierr = NF90_PUT_ATT(icid, NF90_GLOBAL, 'Reference_GEBCO_2022', 'GEBCO Compilation Group (2022) GEBCO 2022 Grid (doi:10.5285/e0f0bb80-ab44-2739-e053-6c86abc0289c) ' )
    ierr = NF90_PUT_ATT(icid, NF90_GLOBAL, 'Reference_BedMachine_Antarctica_v3', &
     &' - Morlighem, M. (2022). MEaSUREs BedMachine Antarctica, Version 3 [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/FPSU0V1MWUB6. Date Accessed 01-11-2023.'//&
     &' - Morlighem, M., E. Rignot, T. Binder, D. D. Blankenship, R. Drews, G. Eagles, O. Eisen, F. Ferraccioli, R. Forsberg, P. Fretwell, V. Goel, J. S. Greenbaum, H. Gudmundsson, J. Guo, V. Helm, C. Hofstede, I. Howat, A. Humbert, W. Jokat, N. B. Karlsson, W. Lee, K. Matsuoka, R. Millan, J. Mouginot, J. Paden, F. Pattyn, J. L. Roberts, S. Rosier, A. Ruppel, H. Seroussi, E. C. Smith, D. Steinhage, B. Sun, M. R. van den Broeke, T. van Ommen, M. van Wessem, and D. A. Young. 2020. Deep glacial troughs and stabilizing ridges unveiled beneath the margins of the Antarctic ice sheet. Nature Geoscience. 13. DOI: 10.1038/s41561-019-0510-8.' )
    ierr = NF90_PUT_ATT(icid, NF90_GLOBAL, 'Reference_BedMachine_Greenland_v5', &
     &' - Morlighem, M. et al. (2022). IceBridge BedMachine Greenland, Version 5 [Data Set]. Boulder, Colorado USA. NASA National Snow and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/GMEVBWFLWA7X. Date Accessed 01-11-2023.'//&
     &' - Morlighem, M., C. Williams, E. Rignot, L. An, J. E. Arndt, J. Bamber, G. Catania, N. Chauché, J. A. Dowdeswell, B. Dorscheel, I. Fenty, K. Hogan, I. Howat, A. Hubbard, M. Jakobsson, T. M. Jordan, K. K. Kjeldsen, R. Millan, L. Mayer, J. Mouginot, B. Noel, C. O Cofaigh, S. J. Palmer, S. Rysgaard, H. Seroussi, M. J. Siegert, P. Slabon, F. Straneo, M. R. van den Broeke, W. Weinrebe, M. Wood, and K. Zinglersen. 2017. BedMachine v3: Complete bed topography and ocean bathymetry mapping of Greenland from multi-beam echo sounding combined with mass conservation. Geophysical Research Letters. 44. DOI: 10.1002/2017GL074954.' )


    ierr = NF90_ENDDEF(icid)
    ! write variables
    ierr = NF90_PUT_VAR(icid, idvlon, rco40(:,:,1) )
    ierr = NF90_PUT_VAR(icid, idvlat, rco40(:,:,2) )
    ierr = NF90_PUT_VAR(icid, idvbat, rbat40(:,:) )
    ierr = NF90_PUT_VAR(icid, idvdra, rdra40(:,:) )
    ierr = NF90_CLOSE(icid)

  END SUBROUTINE WriteOutput



END PROGRAM solve_puzzle
