        subroutine sode15(NEQ, NCTR, T0, Y0, u, atol,
     1       rtol, Fcn, Jcn)
            implicit none

c           mkl function
            double precision dlange


c           Uebergabeparameter
            external Fcn, Jcn

            integer NEQ, NCTR
            double precision T0, Y0, u, atol, rtol

c           Zaehlervariable
            integer Iter, m, k, j, klast, kopt, inc

            integer nfevals, nsolves, nfailed, npds,
     1           ndecomps, nsteps

            logical Jcurrent

c           helper Variablen integer
            integer E, maxIt, maxk, nconhk, difU, tdir
c           helper Variable double precision
            double precision y, yp, t, dfdy, threshold, anaPD, hmin,
     1           wt, rh, htspan, tfinal, hmax, absh, h, DfDt, ml, mu,
     2           abshlast, dif


c           parameter double precision
c           Gamma_k, alpha, invGa
            double precision G, alpha, invGa, erconst, eps

c            parameter declartion integer
            parameter(maxIt = 4, maxk = 5, inc = 1)

c           Uebergabeparameter dimension
            dimension T0(2), Y0(NEQ), u(NCTR)

c           helper dimension integer
            dimension E(NEQ, NEQ)
c           helper dimension double precision
            dimension Y(NEQ), dfdy(NEQ, NEQ), yp(NEQ), wt(NEQ),
     1           DfDt(NEQ), dif(NEQ, maxk+2)


c           parameter dimension
            dimension G(5), alpha(5), invGa(5), erconst(5), difU(5, 5)


c            parameter declartion double precision
c           ToDo: eps bestimmen
            parameter(eps = 2.2204D-16, ml = 1.0, mu= 0.0)

c           parameter declartion array
            parameter(G = (/ 1.0, 3.0/2.0, 11.0/6.0, 25.0/12.0,
     1           137.0/60.0 /), alpha = (/-37.0/200.0, -1.0/9.0,
     2            -0.0823, -0.0415, 0.0/),
     3            invGa = 1.0 / (G * (1.0 - alpha)))

            parameter (difU =(/
     $ -1,  0,  0,  0,  0,
     $ -2,  1,  0,  0,  0,
     $ -3,  3, -1,  0,  0,
     $ -4,  6, -4,  1,  0,
     $ -5,  10,-10, 5, -1/))

            parameter (erconst = (/ alpha(1) * G(1) + 1.0 / 2,
     1       alpha(2) * G(2) + 1.0 / 3, alpha(3) * G(3) + 1.0 / 4,
     2       alpha(4) * G(4) + 1.0 / 5, alpha(5) * G(5) + 1.0 / 6/))

            data k/1/, m/1/, klast/1/

            E = 0
            do Iter = 1, NEQ
                E(Iter, Iter) = 1
            enddo

            threshold = atol / rtol
            htspan = abs(T0(2) - T0(1))
            tfinal = T0(2)
c           By default, hmax = 1/10 of the interval
            hmax = 0.1 * htspan
            tdir = sign(1.0, (T0(2) - T0(1)))

c           Stats
            nsteps   = 0
            nfailed  = 0
            nfevals  = 0
            npds     = 0
            ndecomps = 0
            nsolves  = 0
            dif = 0

            t = T0(1)
            y = Y0

c            yp = 0
c            dfdy = 0

            call Fcn (NEQ, T0(1), Y0, u, yp)
            nfevals = nfevals + 1

            call Jcn(NEQ, T0(1), Y0, u, dfdy, NEQ)
            npds = npds + 1

            Jcurrent = .true.

            hmin = 16 * eps * abs(t)

            wt = max(abs(y),threshold)
c            rh = 1.25 * maxval(abs(yp/wt), NEQ) / sqrt(rtol)
c           ToDo: spaeter absh durch hmax ersetzen
c            absh = min(hmax, htspan)

c            if (absh * rh > 1) then
c                absh = 1 / rh
c            end if
c            absh = max(absh, hmin)

c            tdel = (t + tdir*min(sqrt(eps)*max(abs(t),abs(t+h)),absh))-t
c           ToDo timeinvariant deshalb kein dfdt nötig
c            f1 = feval(odeFcn,t+tdel,y,odeArgs{:});
c            call Fcn (NEQ, t+tdel, Y0, u, f1)
c            nfevals = nfevals + 1;
c            dfdt = (f1 - yp) / tdel;
            DfDt = 0
            call dgemv('N', NEQ, NEQ, ml, dfdy, NEQ, yp, inc,
     1           mu, DfDt,inc)
c
            rh = 1.25 * sqrt(0.5 * maxval(abs( DfDt / wt)) / rtol)
            absh = min(hmax, htspan)

            if (absh * rh .GT. 1) then
                absh = 1 / rh
            end if

            absh = max(absh, hmin)
c Matlab 417
            h = tdir * absh
            abshlast = absh

            dif(:,1) = h * yp
            hinvGak = h * invGa(k)

c Matlab 425

c            % Initialize.
ck = 1;                                  % start at order 1 with BDF1
cK = 1;                                  % K = 1:k
cklast = k;
cabshlast = absh;

cdif = zeros(neq,maxk+2);
cdif(:,1) = h * yp;

chinvGak = h * invGa(k);
cnconhk = 0;                             % steps taken with current h and k

cMiter = Mt - hinvGak * dfdy;



		end subroutine
