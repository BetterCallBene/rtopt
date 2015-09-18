		program rtopt
		use omp_lib
		use mkl_service
cDEC$ DEFINE debug=1		

c		implicit none

c 		Wichtig FEX als external deklarien sonst segment fault
c 		Externe Funktionen (Pointer)
		external FEX, JacFx, dnrm2

c 		Konstanten
		integer n_state, n_contr, n_var, 
	1		n_intervals,  NEQ, inc, n_vec,
	2		n_q1, n_q2

c 		Speichervariablen
		double precision vec, state

c 		Zaehlervariablen
		integer i, j, timepoint, IOUT

c 		current controls
		double precision u	

c 		Initial Matrizen M0, N0
		double precision M0, N0

		double precision Y
		double precision F, M, N

c 		Normen
		double precision dnrmF, nrm2v, invNrm2v

c   	Time and Stepsize
		double precision T, TOUT, h
		double precision tspan

c 		Parameter fuer solver
		integer ITASK, ISTATE, LRW, MF, 
	1		ITOL, IOPT, IWORK, LIW
		double precision RTOL, ATOL, RWORK

c 		Funktionen 
		double precision dnrm2, dlange
c 		Helper
		double precision tmpQ, alpha

c 		cpu time
		double precision start_time, end_time

		double precision anaPD, numPD, sumPD

		double precision nrmMatrix

		integer ML, MU
		

		parameter(n_state = 13, n_var = 17, n_contr = 4, 
	1		n_q1 = 4, n_q2 = 7, n_intervals = 10000, 
	2		inc = 1, n_vec = n_var * (n_intervals + 1),
	3		NEQ = 234, LRW = 22 +  9*NEQ + NEQ**2, 
	4 		LIW = 20 + NEQ, MF = 21)


		parameter(ITOL = 1, RTOL = 1e-3, ATOL = 1e-4, h = 2D-1)

c 		helper
		parameter (alpha = 1)

c 		Speichervariablen
		dimension vec(n_vec), state(n_state), u(n_contr)
		dimension Y(NEQ), RWORK(LRW), IWORK(LIW)

		dimension anaPD(NEQ, NEQ), numPD(NEQ, NEQ), sumPD(NEQ, NEQ)

c		dimension Y(NEQ), YDOT(NEQ), u(n_contr), F(n_state),
c	1		anaPD(NEQ, NEQ), numPD(NEQ, NEQ),
c	2		PD(NEQ, NEQ)

		dimension F(n_state), M(n_state, n_state), N(n_state, n_contr),
	1		J(n_state, n_var)

c 		Initial Matrizen M0, N0
		dimension M0(n_state, n_state), N0(n_state, n_contr)
		dimension tspan(2)

		dimension tmpQ(n_contr)

		common /contr/ u

c 		Init M0, N0
		N0 = 0
		do 10 i = 1, 13
			M0(i, i) = 1.0
10		continue 

		call rand(n_vec, vec)

		

		call cpu_time(start_time)
c$OMP PARALLEL DEFAULT(PRIVATE) SHARED(vec)
c$OMP DO
		do 40 IOUT = 1, n_intervals

			state = vec((IOUT - 1) * n_var +1:
	1			(IOUT - 1) * n_var + n_state)
		
			u = vec((IOUT - 1) * n_var + n_state + 1:
	1				IOUT * n_var)
						
c 		Normierung der Quaternionen
			nrm2v = dnrm2(n_q1, state(n_q1:n_q2), 
	1			inc)

			invNrm2v = 1.0 / nrm2v

			call daxpy(n_q1, invNrm2v, state(n_q1:n_q2), inc, tmpQ, inc)
			state(n_q1:n_q2) = 0
			call daxpy(n_q1, alpha, tmpQ, inc, state(n_q1:n_q2), inc)
c 		Ende
c
			call helperCreateVektor(state, M0, N0, Y)
			call diff_nD(Y, FEX, NEQ, numPD)
c
c			print *, numPD
c
			call JacFx(NEQ, T, Y, ML, MU, anaPD, NEQ)
			sumPD = numPD - anaPD

			nrmMatrix = dlange( 'F', NEQ, NEQ, sumPD, NEQ, 0 )

			print *, 'Norm Matrix=', nrmMatrix
40			continue
c$OMP END DO 		
c$OMP END PARALLEL
		call cpu_time(end_time)
		print *, end_time -start_time
		STOP
80		WRITE(6,90)  ISTATE
90			FORMAT(///' Error halt.. ISTATE =',I3)
		STOP


	end program rtopt


