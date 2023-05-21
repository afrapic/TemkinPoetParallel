module baseLag_module 
    use sturmians_data
    use sturmians_base 
    use complex_module
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !vers: 8 Jun. 2009
    !Modificado para matrices reales resuelve con DSYGVD para un sistema de la forma A*x=(beta)*B*x
    !con A y B reales simetricas y B definida positiva. 
    !Para onda saliente el sistema es de la misma forma con A(N,N) complejo. Se resuelve aplicando
    !factorizacion de Cholesky en B=U**T*U, luego Ap=U**(-T)*A*U**(-1) y se llega a un sistema
    !de la forma Ap*y=beta*y con y=U*x. Luego de esta transformacion, la matriz Ap es simetrica,
    !y el unico elemento complejo sigue siendo Ap(N-1,N-1). Se resuelve con un método general
    !ya que la matriz es compleja no hermítica
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    subroutine base_Laguerre(n,b_n,ele)
        use sturmians_data
        use sturmians_base 
	implicit none

	type(base_lag) :: b_n
	integer i,j,ele,n,info,order,lwork
	integer,allocatable :: iwork(:)
	real*8 hcplx,abnrm
        real*8,allocatable :: rgrid(:)
	real*8,allocatable :: h(:,:),v(:,:),w(:),work(:),temp(:,:),temp1(:,:)
	real*8,allocatable :: scale(:),rconde(:),rcondv(:)
	complex*16 En
	complex*16,allocatable :: normP(:), vl(:,:),eigenvalues(:),an(:,:),hp(:,:),zwork(:)

 !Grilla radial
        allocate(rgrid(0:nrmax(n)-1))
	do i=0,nrmax(n)-1
            rgrid(i)=i*rmax(n)/(nrmax(n)*1d0-1d0)
        enddo
        bs(n)%r=rgrid
        deallocate(rgrid)

        !Energía compleja
	En=E(n)*c1

 !Tamaño de la base
       	order=norder(n,ele)

	!Construye la matriz del hamiltoniano (h) y del potencial (v)
        allocate(v(0:order-1,0:order-1),h(0:order-1,0:order-1))
	call matrix_v(n,ele,order,v)
	call matrix_h(n,En,ele,order,h,hcplx)

 !Guarda la matriz del potencial
	allocate(temp(0:order-1,0:order-1))
	temp=v
 !Y multiplica el potencial *(-1)
	v=(-1d0)*v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !CALCULO DE AUTOVALORES Y AUTOVECTORES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Para energia negativa y energia positiva con condicion estacionaria o caja las matrices son
 !las dos reales, y se utiliza DSYGVD para resolver el sistema. 
 !Al finalizar guarda los autovalores en eigenvalues (complex) y los autovectores en an (complex)
 	if(wave(n).ne.'out')then 
            lwork=1+6*order+2*order**2
            allocate(w(0:order-1),work(lwork),iwork(3+5*order))    
            call DSYGVD(1,'V','U',order,h,order,v,order,w,work,lwork,iwork,3+5*order,info)
            deallocate(v,work,iwork)
            if(info.eq.0)then
                allocate(eigenvalues(0:order-1))
                eigenvalues=w+ci*0d0
                allocate(an(0:order-1,0:order-1))
                an=h+ci*0d0
                deallocate(h,w)
            else
                write(*,*) 'Error al calcular los autovalores con DSYGVD'
                stop
            endif
	endif

 !Para energia onda saliente la matriz (-v) es real y definida positiva y h es real simétrica 
 !excepto en h(N-1,N-1) donde vale h(n-1,n-1)+im*hcplx. 
 !Al finalizar guarda los autovalores en eigenvalues (complex) y los autovectores en an(complex)
 	if(wave(n).eq.'out')then 
      !Factorizacion de Cholesky de -v=u**T*u.
            call DPOTRF('U',order,v,order,info)
            if(info.ne.0)then
                write(*,*) 'Fallo factorizacion Cholesky'
                stop
            endif
            !Transforma h en inv(u**T)*h*inv(u)
            call DSYGST(1,'U',order,h,order,v,order,info)
            if(info.ne.0)then
                write(*,*) 'Fallo transformacion h'
                stop
            endif
            !Transforma la parte compleja
            allocate(temp1(0:order-1,0:order-1))
            temp1=0d0
            temp1(order-1,order-1)=hcplx
            call DSYGST(1,'U',order,temp1,order,v,order,info)
            allocate(hp(0:order-1,0:order-1))
            !Suma parte real y compleja y ahora resuelve hp*y=beta*y
            hp=h+ci*temp1
            deallocate(h,temp1)
            do i=0,order-1
                do j=i,order-1
                    hp(j,i)=hp(i,j)
                enddo
            enddo
            allocate(eigenvalues(0:order-1),vl(0:order-1,0:order-1),an(0:order-1,0:order-1),scale(order))
            allocate(rconde(order),rcondv(order),zwork(2*order+2),work(2*order))
            call ZGEEVX('B','V','V','N',order,hp,order,eigenvalues,vl,order,an,order,i,j,scale,abnrm,rconde,&
                 rcondv,zwork,2*order+2,work,info)
            deallocate(hp,vl,scale,rconde,rcondv,zwork,work)
            if(info.ne.0)then
                write(*,*) 'Error en zgeevx al resolver sistema matricial'
                stop
            endif
            !Transforma los autovectores al sistema original x=U**-1*y
            allocate(hp(0:order-1,0:order-1))
            hp=cmplx(v)
            call ZTRSM('l','U','N','N',order,order,c1,hp,order,an,order)
            deallocate(v,hp)
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!
	if(out.eq.'y')then
	allocate(zwork(0:norder(n,ele)+nasint))
	call coef_Coulomb(norder(n,ele)+nasint,ele,lambda(n,ele),En,Z(n),zwork)
	bs(n)%coef_asint(:,ele)=zwork(:)
	deallocate(zwork)
	endif
!!!!!!!!!!

        !Ordena los autovalores en orden creciente de la parte real
        call Orden(n,ele,order,eigenvalues,an) 
	
        !Calcula la normalizacion
	allocate(normP(0:norder(n,ele)-1))
	if(nrm.eq.'v')then
 	call normPotential(n,ele,order,an,temp,normP)
	endif
	if(nrm.eq.'o')then
	call normOverlap(n,ele,order,an,normP)
	endif
	if(nrm.eq.'n')then
	allocate(zwork(0:norder(n,ele)))
  	call coef_Coulomb(norder(n,ele),ele,lambda(n,ele),En,Z(n),zwork)
	do i=0,norder(n,ele)-1
	normP(i)=zwork(norder(n,ele)-1)/an(norder(n,ele)-1,i)
	enddo
	deallocate(zwork)
	endif
	deallocate(temp)

     !Guarda en bs los autovalores y autovectores normalizados
	do j=1,norder(n,ele)
	bs(n)%beta(j,ele)=eigenvalues(j-1)
		do i=0,norder(n,ele)-1
        	bs(n)%a(i,j,ele)=an(i,j-1)*normP(j-1)
		enddo
        enddo
        deallocate(an,eigenvalues,normP)

	return
    end subroutine base_Laguerre

    !----------------------------------------------------------------------------------------------------------
    subroutine matrix_v(n,ele,ord,matrixv)
	!<fi(jp)|V|fi(j)>, real
        use sturmians_data
        use sturmians_base 
	use gaulag_module
	use gauleg_module
	use inf_nan_detection
        implicit none
	character ich
        integer n,ele,j,jp,ord,k,npoints,info
	real*8,allocatable :: x(:),w(:),v(:,:),Lx(:)
     	real*8 matrixv(0:ord-1,0:ord-1),b
        complex*16 gam1,gam2,one
	one=(1d0,0d0)

        matrixv=0d0

	if(pot(n).eq.'yuk')then   !Yukawa -e^(-alfa*r)/r
            b=1d0+alfa(n,ele)/2d0/lambda(n,ele)
            do jp=0,ord-1
                call gammaln(one*(jp*1d0+2d0*ele+2d0),gam1)
                call gammaln(one*(jp*1d0+1d0),gam2)
                matrixv(jp,0)=real((exp(gam1-gam2)*(b-1d0)**jp)/(b)**(jp+2*ele+2))
                matrixv(0,jp)=matrixv(jp,0)
            enddo
            do jp=1,ord-1
                do j=jp,ord-1
                    matrixv(jp,j)=((1d0*jp)/(1d0*j)+1d0+2d0*ele/(1d0*j)+1d0/(1d0*j))*(b-1d0)*matrixv(jp,j-1)/b &
                         -((1d0*jp)/(1d0*j)+2d0*ele/(1d0*j)+1d0/(1d0*j))*(b-2d0)*matrixv(jp-1,j-1)/b
                    matrixv(j,jp)=matrixv(jp,j)
                    if(isnan(matrixv(jp,j)))then
                        write(*,'(a)') 'Algun elemento de v es NaN'
                        stop
                    endif
                enddo
            enddo
            matrixv=(-1d0)*matrixv
	endif !Yukawa -e^(-alfa*r)/r

	if(pot(n).eq.'exp')then  !exponencial -exp(-alfa*r)
     !Armo la del Yukawa hasta N
            allocate(v(0:ord,0:ord))
            v=(0d0,0d0)
            b=1d0+alfa(n,ele)/2d0/lambda(n,ele)
            do jp=0,ord
                call gammaln(one*(jp*1d0+2d0*ele+2d0),gam1)
                call gammaln(one*(jp*1d0+1d0),gam2)
                v(jp,0)=real((exp(gam1-gam2)*(b-1d0)**jp)/(b)**(jp+2*ele+2))
                v(0,jp)=v(jp,0)
            enddo
            do jp=1,ord
                do j=jp,ord
                    v(jp,j)=((1d0*jp)/(1d0*j)+1d0+2d0*ele/(1d0*j)+1d0/(1d0*j))*(b-1d0)*v(jp,j-1)/b &
                         -((1d0*jp)/(1d0*j)+2d0*ele/(1d0*j)+1d0/(1d0*j))*(b-2d0)*v(jp-1,j-1)/b
                    v(j,jp)=v(jp,j)
                enddo
            enddo
            !Construyo la del exponencial
            do jp=0,ord-1
                matrixv(jp,0)=(ele+1d0)*v(jp,0)/lambda(n,ele)-v(jp,1)/2d0/lambda(n,ele)
                matrixv(0,jp)=matrixv(jp,0)
            enddo
            do jp=1,ord-1
                do j=jp,ord-1
                    matrixv(jp,j)=(j+ele+1d0)*v(jp,j)/lambda(n,ele)-(j+1d0)*v(jp,j+1)/2d0/lambda(n,ele)-(j+2d0*ele+1d0)*v(jp,j-1)/2d0/lambda(n,ele)
                    matrixv(j,jp)=matrixv(jp,j)
                    if(isnan(matrixv(jp,j)))then
                        write(*,'(a)') 'Algun elemento de v es NaN'
                        stop
                    endif
                enddo
            enddo
            deallocate(v)
            matrixv=matrixv*(-1d0)
	endif

	if(pot(n).eq.'cou')then   !Coulombiano -Zi/r
            do jp=0,ord-1
                call gammaln(one*(jp*1d0+2d0*ele+2d0),gam1)
                call gammaln(one*(jp*1d0+1d0),gam2)
                matrixv(jp,jp)=real(-Zi(n,ele)*exp(gam1-gam2))
                if(isnan(matrixv(jp,j)))then
                    write(*,'(a)') 'Algun elemento de v es NaN'
                    stop
                endif
            enddo
	endif

	if(pot(n).eq.'hul')then   !Hulten -exp(-r/alfa)/(1-exp(-r/alfa))
            npoints=200
            allocate(x(npoints),w(npoints))
            call gaulag(x,w,2d0*ele+1d0)
            allocate(Lx(0:ord))
            b=2d0*lambda(n,ele)*alfa(n,ele)
            do jp=0,ord-1
                do j=jp,ord-1
                    do k=1,npoints
                        call laguerre_general(ord,2d0*ele+1d0,x(k),Lx)
                        matrixv(jp,j)=matrixv(jp,j)+w(k)*x(k)*Lx(jp)*Lx(j)*exp(-x(k)/b)/(1d0-exp(-x(k)/b))/2d0/lambda(n,ele)
                    enddo
                    matrixv(j,jp)=matrixv(jp,j)
                    if(isnan(matrixv(jp,j)))then
                        write(*,'(a)') 'Algun elemento de v es NaN'
                        stop
                    endif
                enddo
            enddo
            matrixv=matrixv*(-1d0)
            deallocate(x,w,Lx)
	endif

	if(pot(n).eq.'cbx')then   !Coulombiano cortado -1/r si r<r0; 0 si r>r0
            npoints=200
            allocate(x(npoints),w(npoints))
            call gauleg(0d0,r0(n,ele),x,w)
            allocate(Lx(0:ord-1))
            do jp=0,ord-1
                do j=jp,ord-1
                    do k=1,npoints
                        call basis(ord,ele,lambda(n,ele),x(k),Lx)
                        matrixv(jp,j)=matrixv(jp,j)+w(k)*Lx(jp)*Lx(j)*(-1d0)/x(k)
                    enddo
                    matrixv(j,jp)=matrixv(jp,j)
                    if(isnan(matrixv(jp,j)))then
                        write(*,'(a)') 'Algun elemento de v es NaN'
                        stop
                    endif
                enddo
            enddo
            deallocate(x,w,Lx)
	endif

	if(pot(n).eq.'vbx')then   !pozo esférico -1d0 si r<r0; 0 si r>r0
            npoints=200
            allocate(x(npoints),w(npoints))
            call gauleg(0d0,r0(n,ele),x,w)
            allocate(Lx(0:ord-1))
            do jp=0,ord-1
                do j=jp,ord-1
                    do k=1,npoints
                        call basis(ord,ele,lambda(n,ele),x(k),Lx)
                        matrixv(jp,j)=matrixv(jp,j)+w(k)*Lx(jp)*Lx(j)*(-1d0)
                    enddo
                    matrixv(j,jp)=matrixv(jp,j)
                    if(isnan(matrixv(jp,j)))then
                        write(*,'(a)') 'Algun elemento de v es NaN'
                        stop
                    endif
                enddo
            enddo
            deallocate(x,w,Lx)
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Chequea que la matriz sea definida positiva
	allocate(v(0:ord-1,0:ord-1),w(0:ord-1),x(3*ord-1))
	v=(-1d0)*matrixv
	call DSYEV('N','U',ord,v,ord,w,x,3*ord-1,info)
	deallocate(v,x)
	if(info.eq.0)then
            if(minval(w).lt.0d0)then
                write(*,*) 'La matriz (-V) del potencial no es definida positiva.'
                stop
            endif
	else
            write(*,'(a5,1x,i5,1x,a)') 'info=',info,'.No se pudo chequear si (-V) es definida positiva.'
            stop
	endif
 	deallocate(w)

	return
    end subroutine matrix_v

    !-----------------------------------------------------------------------------------------------
    subroutine matrix_h(i,En,ele,ord,matrixm,hcomplex)
	!<fi(jp)|H0-E|fi(j)>, con el ultimo elemento de la matriz complejo
	!si la condicion es de onda saliente. La matriz matrizm guarda la parte
	!real en todos los casos y hcomplex es la parte imaginaria de matrixm(ord-1,ord-1),
	!para energia negativa y onda estacionaria es cero.
        use sturmians_data
        use sturmians_base 
	implicit none
	integer i,ele,j,jp,ord
	real*8 x,matrixm(0:ord-1,0:ord-1),hcomplex
	real*8,allocatable :: jac(:,:)
	complex*16,allocatable :: hn(:)
	complex*16 En

	hcomplex=0d0
	allocate(jac(0:ord,0:ord))
	call matrix_jacobi(ord,ele,lambda(i,ele),En,Z(i),jac)

 !Llena la matriz con los elementos correspondientes
	matrixm=0d0
	do jp=0,ord-1
            do j=jp,ord-1
                matrixm(jp,j)=jac(jp,j)
                matrixm(j,jp)=matrixm(jp,j)
            enddo
	enddo

	if(E(i).gt.0d0)then
            if(wave(i).eq.'out')then
                allocate(hn(0:ord))
                call coef_Coulomb(ord,ele,lambda(i,ele),En,Z(i),hn)
                matrixm(ord-1,ord-1)=real(matrixm(ord-1,ord-1)+jac(ord-1,ord)*hn(ord)/hn(ord-1))
                hcomplex=aimag(jac(ord-1,ord-1)+jac(ord-1,ord)*hn(ord)/hn(ord-1))
                deallocate(hn)
            endif
	endif

	if(out.eq.'y')then
            if(E(i).lt.0d0)then
                allocate(hn(0:ord))
                call coef_Coulomb(ord,ele,lambda(i,ele),En,Z(i),hn)
                matrixm(ord-1,ord-1)=matrixm(ord-1,ord-1)+real(jac(ord-1,ord)*hn(ord)/hn(ord-1))
                deallocate(hn)
            endif
	endif

	if(out.eq.'y')then
            if(E(i).gt.0d0)then
                if(wave(i).eq.'sta')then
                    allocate(hn(0:ord))
                    call coef_Coulomb(ord,ele,lambda(i,ele),En,Z(i),hn)
                    matrixm(ord-1,ord-1)=matrixm(ord-1,ord-1)+jac(ord-1,ord)*aimag(hn(ord))/aimag(hn(ord-1))
                    deallocate(hn)
                endif
            endif
	endif

	deallocate(jac)

	return
    end subroutine matrix_h

    !----------------------------------------------------------------------------------------------------------------
    subroutine matrix_jacobi(ord,ele,lda,En,Zi,jac)
	implicit none
	integer i,ord,ele,j,jp
	real*8 Zi,lda,jac(0:ord,0:ord)
	complex*16 En,x,gam1,gam2,one,im

	one=(1d0,0d0)
	im=(0d0,1d0)
	jac=(0d0,0d0)

	x=(En-lda*lda/2d0)/(lda*lda/2d0+En)

	do jp=0,ord
            do j=jp,ord
                If(jp.eq.j)then
                    jac(jp,j)=real(2d0*(j+ele+1d0)*x+2d0*lda*Zi/(-lda*lda/2d0-En))
                endif
                If(jp.eq.j+1)then
                    jac(jp,j)=(-1d0)*(2d0*ele+2d0+j)
                endif
                If(jp.eq.j-1)then
                    jac(jp,j)=(-1d0)*j
                endif
                call gammaln(one*(j*1d0+2d0*ele+2d0),gam1)
                call gammaln(one*(j*1d0+1d0),gam2)
                jac(jp,j)=real((-1d0)*jac(jp,j)*exp(gam1-gam2)*(lda*lda/2d0+En)/2d0/lda)
                jac(j,jp)=jac(jp,j)
            enddo
	enddo

	return
    end subroutine matrix_jacobi
    !--------------------------------------------------------------------------------------------------------------------------
    subroutine coef_Coulomb(ord,l,lda,En,Zi,hn)
	use pi_module
	implicit none
	integer i,ord,l,n
	real*8 Zi,lda
	real*8,allocatable :: jac(:,:)
	complex*16 En,hn(0:ord),k,t,x,nu,tita,sigma,hyp0,hyp1,gama,one,im,aj,bj
	complex*16,external :: f21
	complex*16,allocatable :: rn(:)
	one=(1d0,0d0)
	im=(0d0,1d0)
	hn=(0d0,0d0)

	if(real(En).gt.0d0)then
            k=sqrt(2d0*En)
            t=-Zi/k
            nu=k/lda/2d0
            x=(-1d0)*(0.5d0-im*nu)/(0.5d0+im*nu)
            tita=atan2(aimag(x),real(x))

            call gammaln(one*(l+1d0-im*t),gama)
            sigma=atan2(aimag(exp(gama)),real(exp(gama)))	
            hyp0=f21(one*(-l-im*t),one,one*(l+2d0-im*t),exp((-2d0)*im*tita))
            hyp1=f21(one*(-l-im*t),one*(2d0),one*(1d0+l+2d0-im*t),exp((-2d0)*im*tita))

            call gammaln(one*(l+2d0-im*t),gama)	
            hn(0)=(-1d0)*exp((-1d0)*im*tita)*hyp0*exp(-gama)
            call gammaln(one*(1d0+l+2d0-im*t),gama)
            hn(1)=(-1d0)*exp((-2d0)*im*tita)*hyp1*exp(-gama)

            allocate(jac(0:ord,0:ord))
            call matrix_jacobi(ord,l,lda,En,Zi,jac)

            do n=1,ord-1
                hn(n+1)=-exp(log(jac(n,n))-log(jac(n,n+1)))*hn(n)-exp(log(jac(n,n-1))-log(jac(n,n+1)))*hn(n-1)
            enddo
            deallocate(jac)

            hn=hn*exp(im*sigma)*exp(pi*t/2d0)*exp(-tita*t)/(2d0*sin(tita))**l
	endif

	if(real(En).lt.0d0)then 
            k=im*sqrt(-2d0*En)
            t=-Zi/k
            nu=k/lda/2d0
            x=(-1d0)*(0.5d0-im*nu)/(0.5d0+im*nu)
            tita=-im*log(x)

            call gammaln(one*(l+1d0-im*t),gama)
            sigma=atan2(aimag(exp(gama)),real(exp(gama)))	
            hyp0=f21(one*(-l-im*t),one,one*(l+2d0-im*t),exp((-2d0)*im*tita))
            hn(0)=(-1d0)*exp((-1d0)*im*tita)*hyp0*exp(-gama)*exp(im*sigma)*exp(pi*t/2d0)*exp(-tita*t)/(2d0*sin(tita))**l 

            allocate(rn(0:ord))
            call continued_fraction_hypgeo(-l-im*t,one*(ord-1d0+1d0),one*(ord-1d0+l+2d0-im*t),exp(-2d0*im*tita),rn(ord))
            rn(ord)=rn(ord)*exp(-im*tita)/(1d0+l/(ord)+1d0/(ord*1d0)-im*t/(ord*1d0))

            allocate(jac(0:ord,0:ord))
            call matrix_jacobi(ord,l,lda,En,Zi,jac)
            do n=ord-1,1,-1
                rn(n)=-one/(exp(log(jac(n,n+1))-log(jac(n,n-1)))*rn(n+1)+exp(log(jac(n,n))-log(jac(n,n-1))))
            enddo

            do n=1,ord
                hn(n)=rn(n)*hn(n-1)
            enddo

            if(Zi.ne.0d0)then
                hn=-im**(-im*t)*hn
            else
                if(abs(k).gt.lda)then
                    hn=hn*(-1d0)
                endif
            endif

            deallocate(rn,jac)
            hn=real(hn)+im*0d0

	endif

	return
    end subroutine coef_Coulomb

    !-------------------------------------------------------------------------------------------------------------
    Subroutine normPotential(n_el,l,ord,an,matrixv,normP)
        use sturmians_data
        use sturmians_base 
	implicit none
 !Normaliza <y(bi)|V|y(bi)>=1
	integer n_el,i,n,m,ord,l
	real*8 matrixv(0:ord-1,0:ord-1)
	complex*16 an(0:ord-1,0:ord-1)
	complex*16 normP(0:norder(n_el,l)-1)

	normP=(0d0,0d0)

	do i=0,norder(n_el,l)-1
            do n=0,norder(n_el,l)-1
                do m=0,norder(n_el,l)-1
                    normP(i)=normP(i)+an(n,i)*an(m,i)*matrixv(n,m)
                enddo
            enddo
	enddo

        normP=1d0/normP
	normP=sqrt(normP)

	return
    end subroutine  normPotential

    !------------------------------------------------------------------------------------------
    Subroutine normOverlap(n_el,l,ord,an,normP)
        use sturmians_data
        use sturmians_base 
	implicit none
 !Normaliza <y(bi)|y(bi)>=1
	integer n_el,i,j,n,m,ord,l
	complex*16,allocatable :: temp(:,:)
	complex*16 an(0:ord-1,0:ord-1),normP(0:ord-1),one,gam1,gam2
	one=(1d0,0d0)

	allocate(temp(0:ord-1,0:ord-1))
	temp=(0d0,0d0)
 	do i=0,ord-1
            do j=i,ord-1
                if(i.eq.j)then
                    call gammaln(one*(i+2*l+2),gam1)
                    call gammaln(one*(i+1),gam2)
                    temp(i,i)=exp(gam1-gam2)*(i+l+1d0)/(lambda(n_el,l))
                endif
                if(i.eq.j+1)then
                    call gammaln(one*(j+2*l+3),gam1)
                    call gammaln(one*(j+2),gam2)
                    temp(i,j)=(-1d0)*(j+1d0)*exp(gam1-gam2)/2d0/lambda(n_el,l)
                endif
                if(i.eq.j-1)then
                    call gammaln(one*(j+2*l+1),gam1)
                    call gammaln(one*(j),gam2)
                    temp(i,j)=(-1d0)*(j+2d0*l+1d0)*exp(gam1-gam2)/2d0/lambda(n_el,l)
                endif
                temp(j,i)=temp(i,j)
            enddo
 	enddo

	normP=(0d0,0d0)
	do i=0,norder(n_el,l)-1
            do n=0,norder(n_el,l)-1
                do m=0,norder(n_el,l)-1
                    normP(i)=normP(i)+an(n,i)*an(m,i)*temp(n,m)
                enddo
            enddo
	enddo
	deallocate(temp)

        normP=1d0/normP
	normP=sqrt(normP)

	return
    end subroutine  normOverlap

    !-------------------------------------------------------------------------------------------------------------
    subroutine Orden(n,l,ord,eigenvalues,eigenvectors)
        use sturmians_data
        use sturmians_base 
	use inf_nan_detection
	implicit none
	integer i,n,ord,l,nsize
	integer,allocatable :: indx(:)
	real*8,allocatable :: temp(:)
	logical,allocatable :: mask(:)
	complex*16 eigenvalues(0:ord-1),eigenvectors(0:ord-1,0:ord-1),im
	complex*16,allocatable :: eigenvect_temp(:,:),eigenval_temp(:)

	allocate(temp(0:ord-1))
	do i=0,ord-1
            temp(i)=abs(eigenvalues(i))
	enddo

	allocate(indx(0:ord-1))
	call indexx(ord,temp,indx)
	deallocate(temp)
	allocate(eigenvect_temp(0:ord-1,0:ord-1),eigenval_temp(0:ord-1))
	eigenvect_temp=eigenvectors
	do i=0,ord-1
            eigenval_temp(i)=eigenvalues(indx(i)-1)
            eigenvectors(:,i)=eigenvect_temp(:,indx(i)-1)
	enddo
	eigenvalues=eigenval_temp
 	deallocate(indx) 
	deallocate(eigenvect_temp,eigenval_temp)

	return
    end subroutine Orden

    !------------------------------------------------------------------------------------------
    subroutine basis(ord,l,lda,rr,fi)
	use complex_module
	implicit none

	integer ord,l,n
	real*8 rr,lda,fi(0:ord-1)
	real*8,allocatable :: Lx(:)

	allocate(Lx(0:ord))
	call laguerre_general(ord,2d0*l+1d0,2d0*lda*rr,Lx)

	do n=0,ord-1
            fi(n)=exp(-lda*rr)*Lx(n)*(2d0*lda*rr)**(l+1d0)
	enddo

	deallocate(Lx)

	return
    end subroutine basis

    !------------------------------------------------------------------------------------------

end module baseLag_module


