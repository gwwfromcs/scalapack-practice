! Matrix-Matrix Multiply, AB = C using LibSci ScaLAPACK
! filename:  test_scalapack.f90
! compile on NERSC edison: ftn -fast test_scalapack.f90 -o scalpack-test2.exe

! input:     ABCp.dat
!              prow    number of rows in proc grid
!              pcol    number of columns in proc grid
!              n       number of rows/columns in matrix A
!              nb      matrix distribution block size   
! oputput:   fort.u, where u=1000+processor number, and stdout

      implicit none

      integer :: n, m, k, nb, mb, kb    ! problem size and block size
      integer :: myunit                 ! local output unit number
      integer :: myArows, myAcols       ! size of local subset of global array
      integer :: myBrows, myBcols       ! size of local subset of global array
      integer :: myCrows, myCcols       ! size of local subset of global array
      integer :: i,j, igrid,jgrid, iproc,jproc, myi,myj, p
      real*8, dimension(:,:), allocatable :: myA, myB, myC
      real*8, dimension(:,:), allocatable :: A, B, C

      integer :: numroc   ! blacs routine
      integer :: me, procs, icontxt, prow, pcol, myrow, mycol  ! blacs data
      integer :: info     ! scalapack return value
      integer, dimension(9)   :: ides_a, ides_b, ides_c ! scalapack array desc

! Read problem description

      open(unit=1,file="ABCp2.dat",status="old",form="formatted")
      read(1,*) prow
      read(1,*) pcol
      read(1,*) n
      read(1,*) m
      read(1,*) k
      read(1,*) nb
      read(1,*) mb
      read(1,*) kb
      close(1)

      if (((n/nb) < prow) .or. ((n/nb) < pcol)) then
         print *,"Problem size too small for processor set!"
         stop 100
      endif

! Initialize blacs processor grid

      call blacs_pinfo   (me,procs)
      call blacs_get     (0, 0, icontxt)
      call blacs_gridinit(icontxt, 'R', prow, pcol)
      call blacs_gridinfo(icontxt, prow, pcol, myrow, mycol)
      myunit = 1000 + me
      write(myunit,*) "--------"
      write(myunit,*) "Output for processor ",me," to unit ",myunit
      write(myunit,*) "Proc ",me,": myrow, mycol in p-array is ", &
         myrow, mycol

! First call the serial version dgemm do the matrix multiplication 
      allocate (A(n, k))
      allocate (B(m, k))
      allocate (C(n, m))
      do i=1,n
        do j=1,k
          A(i,j) = real(i+j)/real(n*k)
        enddo
      enddo
      do i=1,m
        do j=1,k
          B(i,j) = real(i-j)/real(m*k)
        enddo
      enddo
      do i=1,n
        do j=1,m
          C(i,j) = 0.0
        enddo
      enddo

      call dgemm('N','T',n,m,k,1.0d0,A,n,B,m,0.d0,C,n)

      do i=1,n
        write(myunit,*) (C(i,j), j=1,m)
      enddo

      deallocate(A)
      deallocate(B)
      deallocate(C)

! Construct local arrays
! Global structure:  matrix A of n rows and n columns
!                    matrix B of n rows and n column
!                    matrix C of n rows and n column

      myArows = numroc(n, nb, myrow, 0, prow)
      myAcols = numroc(k, kb, mycol, 0, pcol)

      myBrows = numroc(m, mb, myrow, 0, prow)
      myBcols = numroc(k, kb, mycol, 0, pcol)

      myCrows = numroc(n, nb, myrow, 0, prow) 
      myCcols = numroc(m, mb, mycol, 0, pcol)

      write(myunit,*) "Size of global array is ", n, " x ", n
      write(myunit,*) "Size of block is        ", nb, " x ", nb
      write(myunit,*) "Size of local array is  ", myArows, " x ", myAcols
      write(myunit,*) "Size of local array is  ", myBrows, " x ", myBcols
      write(myunit,*) "Size of local array is  ", myCrows, " x ", myCcols

! Initialize local arrays    

      allocate(myA(myArows,myAcols)) 
      allocate(myB(myBrows,myBcols)) 
      allocate(myC(myCrows,myCcols)) 

      do i=1,n
         call g2l(i,n,prow,nb,iproc,myi)
         if (myrow==iproc) then
            do j=1,k
            call g2l(j,m,pcol,mb,jproc,myj)
               if (mycol==jproc) then
                  myA(myi,myj) = real(i+j)/real(n*k)
                  write(myunit,*)"A(",i,",",j,")", &
                                 " --> myA(",myi,",",myj,")=",myA(myi,myj), &
                                 "on proc(",iproc,",",jproc,")"

               endif
            enddo
         endif
      enddo

      do i=1,m
         call g2l(i,m,prow,mb,iproc,myi)
         if (myrow==iproc) then
            do j=1,k
            call g2l(j,k,pcol,kb,jproc,myj)
               if (mycol==jproc) then
                  myB(myi,myj) = real(i-j)/real(m*k)
                  write(myunit,*)"B(",i,",",j,")", &
                                 " --> myB(",myi,",",myj,")=",myB(myi,myj), &
                                 "on proc(",iproc,",",jproc,")"

               endif
            enddo
         endif
      enddo

      do i=1,n
         call g2l(i,n,prow,nb,iproc,myi)
         if (myrow==iproc) then
            do j=1,m
            call g2l(j,m,pcol,mb,jproc,myj)
               if (mycol==jproc) then
                  myC(myi,myj) = 0.0
!                 write(myunit,*)"C(",i,",",j,")", &
!                                " --> myC(",myi,",",myj,")=",myC(myi,myj), &
!                                "on proc(",iproc,",",jproc,")"

               endif
            enddo
         endif
      enddo

! Prepare array descriptors for ScaLAPACK 

      ides_a(1) = 1         ! descriptor type
      ides_a(2) = icontxt   ! blacs context
      ides_a(3) = n         ! global number of rows
      ides_a(4) = k         ! global number of columns
      ides_a(5) = nb        ! row block size
      ides_a(6) = kb        ! column block size
      ides_a(7) = 0         ! initial process row
      ides_a(8) = 0         ! initial process column
      ides_a(9) = myArows   ! leading dimension of local array

      ides_b(1) = 1         ! descriptor type
      ides_b(2) = icontxt   ! blacs context
      ides_b(3) = m         ! global number of rows
      ides_b(4) = k         ! global number of columns
      ides_b(5) = mb        ! row block size
      ides_b(6) = kb        ! column block size
      ides_b(7) = 0         ! initial process row
      ides_b(8) = 0         ! initial process column
      ides_b(9) = myBrows   ! leading dimension of local array

      ides_c(1) = 1         ! descriptor type
      ides_c(2) = icontxt   ! blacs context
      ides_c(3) = n         ! global number of rows
      ides_c(4) = m         ! global number of columns
      ides_c(5) = nb        ! row block size
      ides_c(6) = mb        ! column block size
      ides_c(7) = 0         ! initial process row
      ides_c(8) = 0         ! initial process column
      ides_c(9) = myCrows   ! leading dimension of local array

! Call ScaLAPACK library routine

      call pdgemm('N','T',n,m,k,1.0d0, myA,1,1,ides_a,  &
                    myB,1,1,ides_b,0.d0, &
                    myC,1,1,ides_c )

! Print results

      ! call g2l(n,n,prow,nb,iproc,myi)
      ! call g2l(n,n,pcol,nb,jproc,myj)
      ! if ((myrow==iproc) .and. (mycol==jproc))  &
      !    write(*,*) 'C(',n,n,')=',myC(myi,myj)
      do i = 1, myCrows
         write(myunit, '((e10.3))') (myC(i,j), j = 1, myCcols)
      enddo

! Deallocate the local arrays

      deallocate(myA, myB, myC)

! End blacs for processors that are used

      call blacs_gridexit(icontxt)
      call blacs_exit(0)

      end

! convert global index to local index in block-cyclic distribution

   subroutine g2l(i,n,np,nb,p,il)

   implicit none
   integer :: i    ! global array index, input
   integer :: n    ! global array dimension, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: p    ! processor array index, output
   integer :: il   ! local array index, output
   integer :: im1   

   im1 = i-1
   p   = mod((im1/nb),np)
   il  = (im1/(np*nb))*nb + mod(im1,nb) + 1

   return
   end

! convert local index to global index in block-cyclic distribution

   subroutine l2g(il,p,n,np,nb,i)

   implicit none
   integer :: il   ! local array index, input
   integer :: p    ! processor array index, input
   integer :: n    ! global array dimension, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: i    ! global array index, output
   integer :: ilm1   

   ilm1 = il-1
   i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1

   return
   end

