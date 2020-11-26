module nrecip_module
  
contains

!--------------------------------------------------------------------------
!                       NUMERICAL RECIPES ROUTINES
!--------------------------------------------------------------------------
FUNCTION ran2(idum)
  INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  doubleprecision ran2,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
     NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  !$OMP THREADPRIVATE(iv,iy,idum2)
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
END FUNCTION ran2
! C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

SUBROUTINE hunt(xx,n,x,jlo)
  INTEGER jlo,n
  doubleprecision x,xx(n)
  INTEGER inc,jhi,jm
  LOGICAL ascnd
  ascnd=xx(n).gt.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt
! C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  double precision a(np,np),b(n)
  INTEGER i,ii,j,ll
  double precision sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        enddo
     else if (sum.ne.0.) then
        ii=i
     endif
     b(i)=sum
  enddo
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     enddo
     b(i)=sum/a(i,i)
  enddo
  return
END SUBROUTINE lubksb
! (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

SUBROUTINE ludcmp(a,n,np,indx,d,success)
  INTEGER n,np,indx(n),NMAX
  double precision d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  double precision aamax,dum,sum,vv(NMAX)
  logical,optional :: success
  if(present(success)) success = .true.
  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     enddo
     if (aamax.eq.0.) then
        if(present(success)) then
           success=.false.
           return
        else
           stop 'singular matrix in ludcmp'
        endif
     endif
     vv(i)=1./aamax
  enddo
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
     enddo
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     enddo
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        enddo
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.)a(j,j)=TINY
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        enddo
     endif
  enddo
  return
END SUBROUTINE ludcmp
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

!     --------------------------------------------------------------
!                Brent's algorithm for root finding
!     --------------------------------------------------------------
FUNCTION zbrent(func,x1,x2,tol)
  INTEGER ITMAX
  DOUBLEPRECISION :: zbrent,tol,x1,x2,func,EPSS
  EXTERNAL func
  PARAMETER (ITMAX=100,EPSS=3.d-8)
  INTEGER :: iter
  DOUBLEPRECISION :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  a=x1
  b=x2
  fa=func(a)
  fb=func(b)
  if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
     stop 'root must be bracketed for zbrent'
  endif
  c=b
  fc=fb
  do iter=1,ITMAX
     if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
        c=a
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc).lt.abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.*EPSS*abs(b)+0.5*tol
     xm=.5*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
        zbrent=b
        return
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa
        if(a.eq.c) then
           p=2.*xm*s
           q=1.-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
           q=(q-1.)*(r-1.)*(s-1.)
        endif
        if(p.gt.0.) q=-q
        p=abs(p)
        if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        endif
     else
        d=xm
        e=d
     endif
     a=b
     fa=fb
     if(abs(d) .gt. tol1) then
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb=func(b)
  enddo
  stop 'zbrent exceeding maximum iterations'
  zbrent=b
  return
END FUNCTION zbrent
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

!     --------------------------------------------------------------
!                   NUMERICAL RECIPES ROUTINE: TRIDAG
!     --------------------------------------------------------------
SUBROUTINE tridag(a,b,c,r,u,n)
  INTEGER :: n
  DOUBLEPRECISION :: a(n),b(n),c(n),r(n),u(n)
  INTEGER :: j
  DOUBLEPRECISION :: bet,gam(n+2)
  if(b(1).eq.0.d0)stop 'tridag: rewrite equations'
  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     if(bet.eq.0.d0)stop 'tridag failed'
     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo
  return
END SUBROUTINE tridag
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

!-------------------------------------------------------------------
!                   FUNCTION: TEST VALIDITY OF NUMBER
!
!     0 = Number is okay
!     1 = Number is INF
!     2 = Number is NAN
!-------------------------------------------------------------------
function number_invalid(a)
  implicit none
  doubleprecision a,b,c
  integer number_invalid
  logical div,sub
  !
  b=a*2.d0
  b=b/2.d0
  c=a-1.d100
  !
  div = (b.eq.a)
  sub = (c.lt.a)
  !
  if(div.and.sub) then
     number_invalid = 0
  elseif(div) then
     number_invalid = 1
  else
     number_invalid = 2
  endif
  !
  return
end function number_invalid

end module nrecip_module
