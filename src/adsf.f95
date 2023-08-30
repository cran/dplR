!*******************************************************************
! calculates a cubic spline function to filter a time series.
! the 'stiffness' of spline is specified by parameter 'ls'
!
! n:  number of values in time series
! y:  time series array to be modeled
! res:  array of cubic spline function values computed
! stiffness: stiffness vector (50% frequency response) of spline
!
!*******************************************************************
!
! AGB Nov 2021 
! Modifying Ed Cook's spline2 subroutine to be a little less ugly
! This included removing all the implicit declarations and setting
! implicit none. Also removing archaic calls to dabs, dcos, dsqrt.
! Changed to follow f95 in terms of I/O variable decalration and
! c binding as well to play nicely with C and .Call
!
!*******************************************************************
      subroutine ads_f(y, n, stiffness, res) bind(C, name = "ads_f_")
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind = c_double), intent(in)  :: y !input, nyrs
      integer(kind = c_int), intent(in)  :: n,stiffness !input
      real(kind = c_double), intent(out) :: res  !output
      integer, parameter :: dp = kind(1d0)
      real(dp) :: v,arg,rn,d1,d2,a,c1,c2,p,pct,pi,sum
      integer :: i,j,k,l,m,imncp1,i1,i2,iw,k1,kl,jm1,n1,nm2,nc,nc1,ncp1

      dimension a(9001,4),c1(4),c2(3),p(9001),stiffness(9001),res(n),y(n)
      data c1 /1.d0,-4.d0,6.d0,-2.d0/
      data c2 /0.d0,.33333333333333d0,1.33333333333333d0/
      data pct,pi / .50d0, 3.1415926535897935d0 /
      if(n .lt. 4)then
        res(1)=-9998. !too few points. bail out.
        return
      endif
	  
      nm2=n-2
! here is where the stiffness vector stiffness is converted into a p vector
! that the spline routine actually needs for calculating the spline
      do i=1,nm2
        v=dble(stiffness(i))                                   
        arg=(2.d0*pi)/v
        p(i)=(6.d0*(cos(arg)-1.d0)**2.d0)/(cos(arg)+2.d0)
      enddo

      do i=1,nm2
         do j=1,3
            a(i,j)=c1(j)+p(i)*c2(j)
            a(i,4)=dble(y(i))+c1(4)*dble(y(i+1))+dble(y(i+2))
         end do
      end do
      a(1,1)=c2(1)
      a(1,2)=c2(1)
      a(2,1)=c2(1)
      nc=2
! begin ludapb
      rn=dble(1.d0/(dble(nm2)*16.d0))
      d1=1.d0
      d2=0.d0
      ncp1=nc+1
      if(nc .ne. 0)then
! initialize zero elements
        do i=1,nc
           do j=i,nc
              k=ncp1-j
              a(i,k)=0.d0
           end do
        end do
! 'i' is row index of element being computed
! 'j' is column index of element being computed
! 'l' is row index of previously computed vector
!       being used to compute inner product
! 'm' is column index
      endif
      do i=1,nm2
         imncp1=i-ncp1
         i1=max0(1,1-imncp1)
         do j=i1,ncp1
            l=imncp1+j
            i2=ncp1-j
            sum=a(i,j)
            jm1=j-1
            if(jm1 .gt. 0)then
              do k=1,jm1
                 m=i2+k
                 sum=sum-a(i,k)*a(l,m)
              end do
            endif
            if(j .eq. ncp1)then
              if(a(i,j)+sum*rn .le. a(i,j))then
                res(1)=-9999.
                return
              endif
              a(i,j)=1.d0/sqrt(sum)
! update determinant
              d1=d1*sum
   35         if(abs(d1) .gt. 1.d0)then
                d1=d1*.0625d0
                d2=d2+4.d0
                goto 35
              endif
   47         if(abs(d1) .le. .0625d0)then
                d1=d1*16.d0
                d2=d2-4.d0
                goto 47
                else
                  goto 60
              endif
            endif
            a(i,j)=sum*a(l,ncp1)
   60       continue
         end do
      end do
! end ludapb / begin luelpb
! solution ly = b
      nc1=nc+1
      iw=0
      l=0
      do i=1,nm2
         sum=a(i,4)
         if(nc .gt. 0)then
           if(iw .ne. 0)then
             l=l+1
             if(l .gt. nc) l=nc
               k=nc1-l
               kl=i-l
               do j=k,nc
                 sum=sum-a(kl,4)*a(i,j)
                 kl=kl+1
               end do
               else
               if(sum .ne. 0.d0)iw=1
               endif
             endif
             a(i,4)=sum*a(i,nc1)
      end do
! solution ux = y
      a(nm2,4)=a(nm2,4)*a(nm2,nc1)
      n1=nm2+1
      do i=2,nm2
         k=n1-i
         sum=a(k,4)
         if(nc .gt. 0)then
           kl=k+1
           k1=min0(nm2,k+nc)
           l=1
           do j=kl,k1
             sum=sum-a(j,4)*a(j,nc1-l)
             l=l+1
           end do 
         endif
         a(k,4)=sum*a(k,nc1)
      end do
! end luelpb
      do i=3,nm2
         res(i)=a(i-2,4)+c1(4)*a(i-1,4)+a(i,4)
      end do
      res(1)=a(1,4)
      res(2)=c1(4)*a(1,4)+a(2,4)
      res(n-1)=a(nm2-1,4)+c1(4)*a(nm2,4)
      res(n)=a(nm2,4)
      do i=1,n
         res(i)=y(i)-res(i)
      end do
      return
      end subroutine ads_f
