subroutine Gsweep( M, NumK, k, thresh,  nr, LC , tol , nSSQ,SSQ,DF,VC,SD,C,Ci,Var,info)
implicit none
integer, intent(in) ::  nr
double precision, intent(inout), dimension(nr,nr)::M
integer, intent(in):: NumK
integer, intent(in), dimension(NumK)::k
double precision, intent(in)::thresh
integer, intent(out), dimension(nr):: LC
double precision, intent(in)::tol
integer, intent(in) ::  nSSQ
double precision, intent(out), dimension(nSSQ)::SSQ
double precision, dimension(NumK,NumK):: Mtemp
double precision, dimension(NumK,NumK,nSSQ):: ZAZ
double precision, intent(out), dimension(nSSQ,nSSQ):: C
double precision, intent(out), dimension(nSSQ,nSSQ):: Ci
double precision, dimension(nSSQ,nSSQ):: W
double precision, intent(out), dimension(nSSQ,nSSQ):: Var
integer,intent(out), dimension(nSSQ):: DF 
double precision,intent(out), dimension(nSSQ):: VC 
double precision,intent(out), dimension(nSSQ):: SD 
double precision, dimension(NumK):: VClong 
integer, intent(out) :: info
integer ::  i, j, n, l ! l runs over variables
double precision :: B, CSS, D, ESS, Rsq 
integer :: inter,rows,n1

external dtrtri

!open(unit=11, file="bla.txt")
!write(11,*)k
!close(unit=11)
!return
if (k(1)==0) then ! check for intercept in model
    inter=1
else
    inter=0
end if
rows=int(M(1,1))
if (inter==0) then
    M(1,1)=1.0d0
endif
DF=0
CSS=0.0d0
ZAZ=0.0d0
W=0
l=0
do n=1,NumK  !over columns k to be swept
  if(abs(M(n,n)) < tol)  then 
      LC(n)=1 ! indicate linear dependency of current column
      if ((n==NumK) .OR. (k(min(n,NumK)) /= k(min(n+1,NumK)))) then
          SSQ(l+1)=CSS
          Mtemp =0.0d0
          if(n<NumK) Mtemp((n+1):NumK,(n+1):NumK)=M((n+1):NumK,(n+1):NumK)
          if ( l<nSSQ) ZAZ(:,:,l+1)= Mtemp
          if (l>0) ZAZ(:,:,l)=ZAZ(:,:,l)-Mtemp
          l=l+1
      end if
      cycle !F95 equivalent of c continue
  end if
  if (l>0) DF(l)=DF(l)+1
  ESS=M(n,n)
  D=ESS
  ! sweep 
  do i=n,nr
      if (i==n) then
          M(n:nr,n)= M(n:nr,n)/D
      else
          B=M(n,i)
          if (dabs(B) < tol) cycle
          M(n:nr,i)=M(n:nr,i)-B*M(n:nr,n)
          M(n,i)=-B/D
      end if
  end do
  M(n,n)=1.0d0/D

  CSS=M(nr,nr)
  Rsq=(CSS-ESS)/CSS

  if (Rsq > (1-thresh)) then
      LC(n)=1
  end if

  if ((n==NumK) .OR. (k(min(n,NumK)) /= k(min(n+1,Numk)))) then
      SSQ(l+1)=CSS
      ! calculate the ZAZ matrices (difference of unswept blocks
      ! of M = XTX  after sweeping out one Z at a time)
      Mtemp =0
      if(n<NumK) Mtemp((n+1):NumK,(n+1):NumK)=M((n+1):NumK,(n+1):NumK)
      if ( l<nSSQ) ZAZ(:,:,l+1) = Mtemp
      if (l>0) ZAZ(:,:,l)=ZAZ(:,:,l)-Mtemp
      l=l+1
  end if
end do
DF(nSSQ)=rows-sum(DF(1:(nSSQ-1)))-inter

! calculate the matrix C (SSQ = C * VC) 
l=1
C=0.0d0
do n=2,NumK
if (n==NumK) then 
    n1=n
else
    n1=n+1
end if
    C(:,l)=C(:,l)+ZAZ(n,n,:)
    if ((n==NumK) .OR. (k(n) /= k(n1))) then
        l=l+1
    end if
end do
C(:,nSSQ)=DF


do l=1,nSSQ-1
    SSQ(l)=SSQ(l)-SSQ(l+1)
end do
Ci=C
call dtrtri('U','N',nSSQ,Ci,nSSQ,info) !invert C matrix
! calculate variances and SD's 
VC = matmul(Ci,SSQ)
do l=1,nSSQ
  SD(l)=sqrt(max(0.0d0,VC(l)))
enddo 
! generate a vector containing the variance for each 1:NumK 
VClong(1)=0.0d0
do l=2,NumK
  VClong(l)=VC(k(l))
enddo

! calculate the covariance matrix of the variances
do  i= 1,(nSSQ-1)
    do  j= i,(nSSQ-1)
    W(j,i) =   dot_product(matmul(VClong,(ZAZ(:,:,i)*ZAZ(:,:,j))),VClong)
      W(i,j) = W(j,i)
    enddo
    do j = 1,(nSSQ-1)
      W(i,i) = W(i,i)+2*VC(j)*VC(nSSQ) *C(i,j)
    enddo
    W(i,i) = W(i,i)+VC(nSSQ)*VC(nSSQ)*DF(i)
enddo
  
W(nSSQ,nSSQ) = VC(nSSQ)*VC(nSSQ)*DF(nSSQ)
Var = 2*matmul(matmul(Ci,W),transpose(Ci)) 

     

     
!close(unit=11)

end subroutine Gsweep


