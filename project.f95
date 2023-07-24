program chorin
implicit none
real*4::L,H,delta_x,delta_y,delta_t,simtime,delta_t2,start,finish,a_2,Re,ro
real*4,dimension(:,:),allocatable::u,v,u0,v0,p,p0,x,y
integer::i,j,m,n,tn,nts
open(1,file='input.txt',status='old')
open(2,file='output.txt',status='replace')
read(1,10)m,n
read(1,11)H,L,Re,ro,simtime,delta_t
allocate(x(m,n))
allocate(y(m,n))
allocate(u(m,n))
allocate(u0(m,n))
allocate(v(m,n))
allocate(v0(m,n))
allocate(p(m,n))
allocate(p0(m,n))
delta_x=L/(m-1)
delta_y=H/(n-1)
tn=floor(simtime/delta_t)
delta_t2=simtime-delta_t*tn
a_2=12./ro
!***************************************
do i=1,m
   do j=1,n
      x(i,j)=(i-1)*delta_x
      y(i,j)=(j-1)*delta_y
   end do
end do
!**********************************************
!*************shart avalie******************
u0=0.
v0=0.
p0=0.
    write(2,*) 'VARIABLES =X,Y,U,V,P'
    write(2,*) 'ZONE I=',n,' J=',m
    write(2,*) 'SOLUTIONTIME=',0.0
 do i=1,m
    do j=1,n
    write(2,13)x(i,j),y(i,j),u0(i,j),v0(i,j),p0(i,j)
   end do
end do


!***********left boundary*****
do j=2,n-1
 u0(1,j)=1
 v0(1,j)=0.
 p0(1,j)=10.
 end do
!********* right boundary******
do j=2,n-1
 u0(m,j)=u0(m-1,j)
 v0(m,j)=0.
 p0(m,j)=p0(m-1,j)
end do
!*********top boundary********* 
do i=2,m-1
  u0(i,n)=0.
  v0(i,n)=0.
  p0(i,n)=p0(i,n-1)
end do
!*********bottom boundary******
do i=2,m-1
  u0(i,1)=0.
  v0(i,1)=0.
  p0(i,1)=p0(i,1)
end do
!*****************************
u0=u
v0=v
p0=p
call cpu_time(start)
do nts=1,tn
  write(2,*) 'ZONE I=',n,' J=',m
  write(2,*) 'SOLUTIONTIME=',delta_t*nts
  do i=2,m-1
   do j=3,n-2
!*******U***************
u(i,j)=u0(i,j)+delta_t*(((u0(i+1,j)-2.*u0(i,j)+u0(i-1,j))/delta_x**2+(u0(i,j+1)-2.*u0(i,j)+u0(i,j-1))/delta_y**2)/re&
&-(u0(i,j)*(u0(i+1,j)-u0(i-1,j))*(p0(i+1,j)-p0(i-1,j))/delta_x**2+(u0(i,j+1)-u0(i,j-1)*(v0(i,j+1)-v0(i,j-1)))/4./delta_y**2))
!*******V**************
v(i,j)=v0(i,j)+delta_t*(((v0(i+1,j)-2.*v0(i,j)+v0(i-1,j))/delta_x**2+(v0(i,j+1)-2.*v0(i,j)+v0(i,j-1))/delta_y**2)/re&
&-(v0(i,j)*(v0(i,j+1)-v0(i,j-1))*(p0(i,j+1)-p0(i,j-1))/delta_y**2+(u(i+1,j)-u(i-1,j)*(v0(i+1,j)-v0(i-1,j)))/4./delta_x**2))   
!********P**************
p(i,j)=p0(i,j)-((a_2*delta_t)/(4.*delta_x*delta_y))*(u0(i+1,j)-u0(i-1,j)+v0(i,j+1)-v0(i,j-1))
end do
   end do

 do i=1,m
  do j=1,n
   write(2,13)x(i,j),y(i,j),u(i,j),v(i,j),p(i,j)
 end do
   end do
end do
    if (delta_t2 /= 0.0) then
 write(2,*) 'ZONE I=',n,' J=',m
 write(2,*) 'SOLUTIONTIME=',simtime
do i=2,m-1
   do j=2,n-1
!**********u**********!
u(i,j)=u0(i,j)+delta_t2*(((u0(i+1,j)-2.*u0(i,j)+u0(i-1,j))/delta_x**2+(u0(i,j+1)-2.*u0(i,j)+u0(i,j-1))/delta_y**2)/re&
&-(u0(i,j)*(u0(i+1,j)-u0(i-1,j))*(p0(i+1,j)-p0(i-1,j))/delta_x**2+(u0(i,j+1)-u0(i,j-1)*(v0(i,j+1)-v0(i,j-1)))/4./delta_y**2))
!**********v**********!
v(i,j)=v0(i,j)+delta_t2*(((v0(i+1,j)-2.*v0(i,j)+v0(i-1,j))/delta_x**2+(v0(i,j+1)-2.*v0(i,j)+v0(i,j-1))/delta_y**2)/re&
&-(v0(i,j)*(v0(i,j+1)-v0(i,j-1))*(p0(i,j+1)-p0(i,j-1))/delta_y**2+(u(i+1,j)-u(i-1,j)*(v0(i+1,j)-v0(i-1,j)))/4./delta_x**2))
!************P***********
p(i,j)=p0(i,j)-((a_2*delta_t2)/(4*delta_x*delta_y))*(u0(i+1,j)-u0(i-1,j)+v0(i,j+1)-v0(i,j-1))

end do
  end do

   do i=1,m
    do j=1,n
      u0=u
      write(2,13)x(i,j),y(i,j),u(i,j),v(i,j),p(i,j)
    end do
   end do
 end if

     
!***********************************************!
    call cpu_time(finish)
    write(*,*) "CPU TIME IS= ",finish - start,"S"
    10 format(T9,I6)
    11 format(T9,F15.8)
    13 format(2X,3F25.10,22X,3F25.10,42X,3F25.10,52X,3F25.10,82X,3F25.5)
end program chorin