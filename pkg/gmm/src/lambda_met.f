      subroutine prep(gt, lam, n, q, d2)
      integer n, q, i, info, ip(q)
      double precision gt(n,q), lam(q), d2(q), dd(q,q)
      double precision tmp(n), tmp2(n), tmp3(n,q)

      call dgemv('n', n, q, 1.0d0, gt, n, lam, 1, 0.0d0, tmp, 1)
      tmp = 1/(1+tmp)
      call dgemv('t', n, q, 1.0d0, gt, n, tmp, 1, 0.0d0, d2, 1)
      tmp2 = tmp**2
      do i=1,q
         tmp3(:,i) = -gt(:,i)*tmp2
      end do
      call dgemm('t','n',q,q,n,1.0d0, gt, n, tmp3, n, 0.0d0, dd, q)
      call dgesv(q, 1, dd, q, ip, d2, q, info)
      end
      
      
      subroutine wu(gt, tol, maxit, n, q, k, conv, obj, lam)
      integer n, q, i, conv, maxit
      double precision obj, lam(q), tol, gt(n,q), k
      double precision dif, d2(q), r, tmp(n), tmp2(q)
      
      lam = 0.0d0
      i = 1
      dif = 1.0d0
      do while (dif>tol .and. i <= maxit)
         call prep(gt, lam, n, q, d2)
         dif = maxval(abs(d2))
         r = 1.0d0
         do while (r>0)
            r = 0.0d0
            tmp2 = lam-d2
            call dgemv('n', n, q, 1.0d0, gt, n, tmp2, 1, 0.0d0,
     *                 tmp, 1)
            if (minval(tmp)<=-1) then
               r = r+1
            end if
            if (r>0) then
               d2 = d2/2
            end if
         end do
         lam = tmp2
         i = i+1
      end do
      if (i>=maxit) then
         lam = 0.0d0
         conv = 1
      else
         lam = -lam
         conv = 0
      end if
      obj = sum(log(1+tmp*k))/n
      end

      subroutine ols(x, y, n, m, lwork, nrhs, info, coef)
      integer n, m, lwork, nrhs, info
      double precision x(n,m), y(n,nrhs), work(lwork) 
      double precision coef(m,nrhs)
      double precision xtmp(n,m), ytmp(n,nrhs)

      xtmp = x
      ytmp = y
      call dgels('n', n, m, nrhs, xtmp, n, ytmp, n, work, -1, info)      
      lwork = min(m*n, int(work(1)))
      if (info == 0) then
         call dgels('n',n,m,nrhs,xtmp,n,ytmp,n,work,lwork,info)
         coef = ytmp(1:m,:)
      end if
      end


      subroutine lamcue(gt, n, q, k, lam, pt, obj)
      integer n, q, lwork, info
      double precision gt(n,q), lam(q), one(n), pt(n), obj, k
      
      lwork = q*3
      one = -1.0d0
      call ols(gt, one, n, q, lwork, 1, info, lam)
      call dgemv('n', n, q, 1.0d0, gt, n, lam, 1, 0.0d0,
     *     pt, 1)
      pt = pt*k
      obj = sum(-pt-(pt**2)/2)/n
      pt = 1+pt
      pt = pt/sum(pt)
      end


      subroutine lamcuep(gt, n, q, k, maxit, conv, lam, pt, obj)
      integer n, q, i, maxit, n0, n1, conv, ind(n)
      double precision gt(n,q), lam(q), pt(n), obj, pt0(n), k
      integer wi(n), wni(n)
      logical w(n)

      call lamcue(gt, n, q, k, lam, pt, obj)
      ind = (/ (i, i=1,n) /)
      i = 1
      conv = 0
      w = (pt < 0)
      do while (.not. all(pt >= 0))
         n0 = count(w)
         n1 = n-n0
         wi(1:n1) = pack(ind, .not. w)
         wni(1:n0) = pack(ind, w)
         if (n1 < (q+1)) then
            pt = 1.0d0/n
            conv = 2
            obj = 0.0d0
            lam = 0.0d0
            exit
         end if
         if (i > maxit) then
            pt = 1.0d0/n
            conv = 1
            obj = 0.0d0
            lam = 0.0d0
            exit
         end if
         call lamcue(gt(wi(1:n1),:), n1, q, k, lam,
     *        pt0(1:n1), obj)
         pt(wi(1:n1)) = pt0(1:n1)
         pt(wni(1:n0)) = 0.0d0
         i = i+1
         w = w .or. (pt<0)
      end do
      if (conv == 0) then
         obj = obj*n1/n + dble(n0)/(2*n)
      end if
      end


      
