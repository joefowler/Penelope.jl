C  ****  Files included to simplify compilation.
      INCLUDE 'penelope.f'
      INCLUDE 'pengeom.f'
      INCLUDE 'penvared.f'
      INCLUDE 'rita.f'
      INCLUDE 'timer.f'

C  *********************************************************************
C                       SUBROUTINE PINITW
C  *********************************************************************
      SUBROUTINE PINITW(EPMAX,NMATR,INFO,PMFILE)
c  Wrap the Penelope subroutine PEINIT so that the caller doesn't need
c  to handle Fortran IO units. Good for Julia to call.
      USE PENELOPE_mod ! Has MAXMAT
      USE PENGEOM_mod ! I/O of the geometry subroutines
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      CHARACTER PMFILE(MAXMAT)*20 ! Material data filenames
      print *, 'Here are the pmfile arguments', MAXMAT
      do i=1,nmatr
            print *,i,pmfile(i)
      enddo
      OPEN (UNIT=16) ! Output file
      CALL PEINIT(EPMAX,NMATR,16,INFO,PMFILE) ! Initialises PENELOPE
      CLOSE(UNIT=16)
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE GINITW
C  *********************************************************************
      SUBROUTINE GINITW(PARINP,NPINP,NMATG,NBOD)
c  Wrap the Penelope subroutine GEOMIN so that the caller doesn't need
c  to handle Fortran IO units. Good for Julia to call.
      USE PENGEOM_mod ! I/O of the geometry subroutines
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION PARINP(20)

      PRINT *, 'Reading geometry...'
      OPEN(17,FILE='epma1.geo') ! Geometry definition file
      OPEN(18,FILE='geometry.rep') ! Geometry report filename
      CALL GEOMIN(PARINP,NPINP,NMATG,NBOD,17,18) ! Initialises PENGEOM
      CLOSE(UNIT=17)
      CLOSE(UNIT=18)
      RETURN
      END

C  *********************************************************************
C                       SUBROUTINE MINITW
C  *********************************************************************
      SUBROUTINE MINITW(NM,PARAM)
c  Julia wrapper to set the materials parameters in the PENELOPE_mod module.
      USE PENELOPE_mod
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION PARAM(7*MAXMAT)

      if (NM.gt.MAXMAT) then
          PRINT *, 'NM must not exceed MAXMAT', NM, MAXMAT
          STOP 'NM must not exceed MAXMAT'
      endif

      do i=1,NM
          j = (i-1)*7
          EABS(1,i) = param(j+1)
          EABS(2,i) = param(j+2)
          EABS(3,i) = param(j+3)
          C1(i) = param(j+4)
          C2(i) = param(j+5)
          wcc(i) = param(j+6)
          wcr(i) = param(j+7)
      enddo
      RETURN
      END
