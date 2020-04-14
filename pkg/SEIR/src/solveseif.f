      subroutine rk4(f, parf, nparf, nb, n, k, a, b, y0,
     *     rzerob, atx, type, rzero, x, y, h)
      integer n, i, k, nparf, nb, type
      double precision a, b, y0(k), h, x(n+1), y(n+1,k)
      double precision f1(k), f2(k), f3(k), f4(k), parf(nparf)
      double precision rzerob(nb), atx(nb), rzero(n+1), tmp
      external f
 
      h = (b-a)/n
      y(1,:) = y0
      x(1) = a
      rzero(1) = rzerob(1)
      do i=2,(n+1)
         x(i) = a+h*(i-1)
         call f(x(i-1), y(i-1,:),parf,nb,rzerob, atx, type,
     *        tmp, f1) 
         call f(x(i-1)+h/2,y(i-1,:)+h*f1/2, parf, nb, rzerob, atx, type,
     *        tmp,f2)
         call f(x(i-1)+h/2,y(i-1,:)+h*f2/2, parf, nb, rzerob, atx, type,
     *        tmp, f3)
         call f(x(i), y(i-1,:)+h*f3, parf, nb, rzerob, atx, type,
     *        rzero(i), f4)
         y(i,:) = y(i-1,:) + h*(f1+2*f2+2*f3+f4)/6
      end do
      end

c     y = {S, E, I, R}
c     par = {c, pop, sigma, gamma, nb, rzerob(nb), atx(nb)}    
      subroutine sysseir(x, y, par, nb, rzerob, atx, type, rzero, dy)
      integer nb, type
      double precision par(4), x, y(4), dy(4)
      double precision pop, sigma, gamma, rzero, beta
      double precision rzerob(nb), atx(nb)
      external R0
      c = par(1)
      pop = par(2)
      sigma = par(3)
      gamma = par(4)     
      call R0(x,rzerob, atx, nb, type, rzero)
      beta = rzero*gamma     
      dy(1) = -beta*y(1)*y(3)/pop
      dy(2) = beta*y(1)*y(3)/pop-sigma*y(2)
      dy(3) = sigma*y(2)-gamma*y(3)+c*y(4)*y(3)/pop
      dy(4) = gamma*y(3)-c*y(4)*y(3)/pop
      end      

      subroutine solveseir(n, a, b, c, sigma, gamma, y0,
     *     nb, rzerob, atx, type, rzero, x, y, h)
      integer n, nb, type
      double precision c, sigma, gamma, y0(4), x(n+1)
      double precision y(n+1,4), par(4), h, a, b
      double precision rzerob(nb), atx(nb), rzero(n+1)
      external sysseir
      par(1) = c
      par(2) = sum(y0)
      par(3) = sigma
      par(4) = gamma
      call rk4(sysseir, par, 4, nb, n, 4, a, b, y0,
     *     rzerob, atx, type, rzero, x, y, h)
      end      
      
      subroutine R0(x, rzerob, atx, nb, type, v)
      integer nb, i, type
      double precision x, v, rzerob(nb), atx(nb), m

      if (nb == 1) then
         v = rzerob(1)
      else
         do i=2,nb
            if (x < atx(i) .and. x >= atx(i-1)) then
               if (type == 1) then
                  v=rzerob(i-1)
               else
                  m = (rzerob(i)-rzerob(i-1))/(atx(i)-atx(i-1))
                  v = rzerob(i-1)+m*(x-atx(i-1))
               end if
            end if
         end do
         if (x >= atx(nb)) then
            v = rzerob(nb)
         end if
      end if
      end
    
