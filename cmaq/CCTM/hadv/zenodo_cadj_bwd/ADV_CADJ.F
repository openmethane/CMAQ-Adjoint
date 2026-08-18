      MODULE ADV_CADJ
      IMPLICIT NONE
      REAL, ALLOCATABLE, SAVE :: ADV_RHOJ(:, :, :, :) ! advected RHOJ
      
      CONTAINS
       SUBROUTINE ADV_CADJ_INIT(JDATE, JTIME)
       USE GRID_CONF
       IMPLICIT NONE
       INCLUDE SUBST_IOPARMS     ! I/O parameters definitions
       INCLUDE SUBST_IODECL      ! I/O definitions and declarations
       CHARACTER( 16 ) :: PNAME='ADV_CADJ'
       CHARACTER( 96 ) :: XMSG = ' '
       INTEGER ALLOCSTAT
       INTEGER JDATE, JTIME
       logical, save :: firstime = .true.

c protect allocation from being called twice 
       if (firstime) then
         firstime = .false.
         ALLOCATE(
     &        ADV_RHOJ(1:NCOLS, 1:NROWS, 1:NLAYS, 1:1),
     &        STAT = ALLOCSTAT)
         IF ( ALLOCSTAT .NE. 0 ) THEN
            XMSG = 'Failure allocating ADV_RHOJ'
            CALL M3EXIT(PNAME, JDATE, JTIME, XMSG, XSTAT1)
         END IF
       end if
       END SUBROUTINE ADV_CADJ_INIT
      END MODULE
