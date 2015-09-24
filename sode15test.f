		program sode15test
c		use omp_lib
c		use mkl_service
cDEC$ DEFINE debug=1		

		implicit none

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
		integer i, IOUT

c 		current controls
		double precision u	

c 		Initial Matrizen M0, N0
		double precision M0, N0

		double precision Y, YP0, T

c 		Normen
		double precision nrm2v, invNrm2v

c   	Time and Stepsize
		double precision h

c 		Parameter fuer solver
		integer LRW, MF,
	1		ITOL, LIW
		double precision RTOL, ATOL

c 		Funktionen 
		double precision dnrm2, dlange
c 		Helper
		double precision tmpQ, alpha

c 		cpu time
c		double precision start_time, end_time

		double precision anaPD, numPD, sumPD

		double precision nrmMatrix

		parameter(n_state = 13, n_var = 17, n_contr = 4, 
	1		n_q1 = 4, n_q2 = 7, n_intervals = 3,
	2		inc = 1, n_vec = n_var * (n_intervals + 1),
	3		NEQ = 234, LRW = 22 +  9*NEQ + NEQ**2, 
	4 		LIW = 20 + NEQ, MF = 21)


		parameter(ITOL = 1, RTOL = 1e-7, ATOL = 1e-8, h = 5D-1)

c 		helper
		parameter (alpha = 1)

c 		Speichervariablen
		dimension vec(n_vec), state(n_state), u(n_contr)
		dimension Y(NEQ), T(2), YP0(NEQ)

		dimension anaPD(NEQ, NEQ), numPD(NEQ, NEQ), sumPD(NEQ, NEQ)

c	    Initial Matrizen M0, N0
		dimension M0(n_state, n_state), N0(n_state, n_contr)

		dimension tmpQ(n_contr)

c 		Init M0, N0
		N0 = 0
		do 10 i = 1, 13
			M0(i, i) = 1.0
10		continue 

c		call rand(n_vec, vec)
        call getTestData(vec)
		T = 0
c		call cpu_time(start_time)

		do 40 IOUT = 1, n_intervals

			state = vec((IOUT - 1) * n_var +1:
	1			(IOUT - 1) * n_var + n_state)
		
			u = vec((IOUT - 1) * n_var + n_state + 1:
	1				IOUT * n_var)

	        T(1) = (IOUT - 1) * h
	        T(2) = IOUT * h
						
c 		Normierung der Quaternionen
			nrm2v = dnrm2(n_q1, state(n_q1:n_q2), 
	1			inc)

			invNrm2v = 1.0 / nrm2v

			call daxpy(n_q1, invNrm2v, state(n_q1:n_q2), inc, tmpQ, inc)
			state(n_q1:n_q2) = 0
			call daxpy(n_q1, alpha, tmpQ, inc, state(n_q1:n_q2), inc)
			call helperCreateVektor(state, M0, N0, Y)

			call JacFx(NEQ, T(1), Y, u, anaPD, NEQ)
			call diff_nD1(NEQ, Y, u, FEX, numPD)

			sumPD = numPD - anaPD
			nrmMatrix = dlange( 'F', NEQ, NEQ, sumPD, NEQ, 0 )

			if (nrmMatrix .GE. 1D-3) then
			    goto 60
			end if

			call sode15(NEQ, n_contr, T, Y, u, ATOL,
     1       RTOL, FEX, JacFx)



40			continue
c		call cpu_time(end_time)
c		print *, end_time -start_time
		STOP
60      print *, 'Error Norm Matrix'
        STOP
	end program sode15test


