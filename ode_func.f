	SUBROUTINE  FEX (NEQ, T, Y, u, YDOT)
		INTEGER  NEQ
		
		double precision YDOT, Y, T, F

		integer n_state, n_var, n_contr
		double precision state, Min, Nin, Mout, Nout
		double precision vec, u

		parameter(n_state = 13, n_var = 17, n_contr = 4)


		dimension Y(NEQ), YDOT(NEQ), F(n_state), 
	1		state(n_state), Min(n_state, n_state), Nin(n_state, n_contr),
	2		vec(n_var), u(n_contr), 
	3		Mout(n_state, n_state), Nout(n_state, n_contr)


        call helperCreateMatrizen(Y, state, Min, Nin)

        T = 0
		vec = 0
		YDOT = 0

		call dcopy(n_state, state, 1, vec, 1)
		call dcopy(n_contr, u, 	   1, vec(n_state+1:n_var), 1)

		call func(vec, F)
		
		call kM(vec, Min, Mout)
		call kN(vec, Nin, Nout)

		call helperCreateVektor(F, Mout, Nout, YDOT)

	end SUBROUTINE



	SUBROUTINE JacFx(NEQ, T, Y, u, PD, NROWPD)

c 		Uebergabeparameter
		INTEGER  NEQ, NROWPD
		DOUBLE PRECISION  T, Y, u, PD

c 		Konstanten
		integer n_state, n_var, n_contr, inc, 
	1		n_mx, n_nx
		double precision alpha, beta
c       Zaehlvariablen
		integer i, j
c 		in-Funktionaufrufvariablen
		double precision F, M, N
c 		helper Variablen
		double precision JacM, HessM, Hx, 
	1		Msp, Nsp, Jx, MxTmp, NxTmp

c 		zentrale Speicher
		double precision vec

c 		Indezierierung
		
		integer irowMx_start, irowNx_start, start_index, end_index
c       Initialisierung Konstanten

		parameter(n_state = 13, n_var = 17, n_contr = 4, 
	1		inc = 1, n_mx = n_state * n_state,  
	2		n_nx = n_state * n_contr, alpha = 1.0, 
	3		beta = 0.0)
c 			Konstanten fuer Indezierung
		parameter(irowMx_start = n_state, 
	1		irowNx_start = n_state + n_mx)

c 		Festlegung der Dimension
c 		Uebergabeparameter
		dimension Y(NEQ), PD(NROWPD, NEQ)
c       Der common Parameter
		dimension u(n_contr)		
c 		Rueckgabewerte der in - Funktionsaufrufe
		dimension JacM(n_state, n_var), HessM(n_state, n_var, n_var),
	1		M(n_state, n_state), N(n_state, n_contr), F(n_state)

c 		Helper Variablen	
		dimension Jx(n_state, n_state), Hx(n_state, n_state, n_state),
	1		Msp(n_state), Nsp(n_state),
	2		MxTmp(n_state, n_state), NxTmp(n_state, n_state)

c 		zentrale Speicher
		dimension vec(n_var)

		T = 0
		PD = 0

c 		Generiere die Matrizen M, N		
		call helperCreateMatrizen(Y, F, M, N)

		call dcopy(n_state, F, 1, vec(1:n_state), 1)
		call dcopy(n_contr, u, 1, vec(n_state+1:n_var), 1)

c       Jacobimatrix
		call jac(vec, JacM)
c 		Hessematrix?
		call hess(vec, HessM)

		Jx = JacM(1:n_state, 1:n_state)
		Hx = HessM(:, 1:n_state, 1:n_state)

		do 20 i = 1, n_state
			Msp = M(:, i)
			Nsp = N(:, i)
			MxTmp = 0
			do 10 j = 1, n_state

				call dgemv('N', n_state, n_state, alpha, 
	1				Hx(j, :, :), n_state, Msp, inc, 
	2				beta, MxTmp(:, j), inc)

				call dgemv('N', n_state, n_state, alpha, 
	1				Hx(j, :, :), n_state, Nsp, inc, 
	2				beta, NxTmp(:, j), inc)
				
10			continue

			if (i .LE. 4) then

				start_index = irowNx_start + (i-1)*n_state + 1
				end_index = irowNx_start + i * n_state

				call mkl_domatadd('r', 't', 'n', n_state, n_state,
	1					alpha, NxTmp, n_state, alpha, 
	2					HessM(:, i + n_state, 1:n_state), 
	3					n_state, PD(start_index:end_index, 1:n_state), n_state)

			end if

c 			mkl_dimatcopy(ordering, trans, rows, cols, alpha, ab, lda, ldb)
			call mkl_domatcopy('r', 't', n_state, n_state, 
	1			alpha, MxTmp, n_state, 
	2			PD(irowMx_start + (i-1)*n_state+1:
	3			irowMx_start + i*n_state, 1:n_state), n_state)

			
			start_index = (i-1)*n_state+1
			end_index = i * n_state

			call mkl_domatcopy('r', 'n', n_state, n_state, 
	1			alpha, Jx, n_state, 
	2			PD(start_index:	end_index, 
	3				start_index: end_index), n_state)

20		continue

		do 30 i = 14, 18

			start_index = (i-1)*n_state+1
			end_index = i * n_state

			call mkl_domatcopy('r', 'n', n_state, n_state, 
	1			alpha, Jx, n_state, 
	2			PD(start_index: end_index, 
	3				start_index:end_index), n_state)

 30		continue
	
	end SUBROUTINE


*     Auxiliary routine: printing a matrix.
*
	SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
	CHARACTER*(*)    DESC
	INTEGER          M, N, LDA
	DOUBLE PRECISION A( LDA, * )
*
	INTEGER          I, J
*
	WRITE(*,*)
	WRITE(*,*) DESC
	DO I = 1, M
		WRITE(*,*) ( A( I, J ), J = 1, N )
c		WRITE(*,*) ""
	END DO
*
c 10	FORMAT( 11(:,1X,F6.4) )
	RETURN
	END


c   C:= A * B + C (mxk * k x n = nxm )
	
	SUBROUTINE multiM(m, k, n, A, B, C) 
		INTEGER          m, k, n
		double precision A, B, C
		double precision beta
		parameter (beta= 0.0)
		dimension A(m, k), B(k, n), C(m, n)

		CALL multiM2(m, k, n, A, B, C, beta)

	end SUBROUTINE

	SUBROUTINE multiM2(m, k, n, A, B, C, beta) 
		INTEGER          m, k, n
		double precision A, B, C
		double precision alpha, beta
		parameter (alpha = 1.0)
		dimension A(m, k), B(k, n), C(m, n)

		CALL DGEMM('N','N',m,n,k,alpha,A,m,B,k,beta,C,m)

	end SUBROUTINE

	SUBROUTINE kM(x, Min, Mout)
		integer n_state, n_var
		parameter(n_state = 13, n_var = 17)

		double precision x, Min, Mout
		double precision Jx, PD
		dimension x(n_var), Min(n_state, n_state),
	1		Mout(n_state, n_state)

		dimension PD(n_state, n_var), Jx(n_state, n_state)

		call jac(x, PD)

		Jx = PD(1:n_state, 1:n_state)
		call multiM(n_state, n_state, n_state, Jx, Min, Mout)
		
	end SUBROUTINE

	SUBROUTINE kN(x, Nin, Nout)
		integer n_state, n_var, n_contr
		parameter(n_state = 13, n_var = 17, n_contr = 4)

		double precision x, Nin, Nout, beta
		double precision Jx, PD

		parameter(beta = 1.0)

		dimension x(n_var), Nin(n_state, n_contr),
	2		Nout(n_state, n_contr)
		dimension PD(n_state, n_var), Jx(n_state, n_state)

		call jac(x, PD)
		
		Jx = PD(1:n_state, 1:n_state)
		Nout = PD(1:n_state, n_state+1:n_var)
c 		(m, k, n, A, B, C, beta) 
		call multiM2(n_state, n_state, n_contr, Jx, Nin, Nout, beta)
c       CALL DGEMM('N','N',	m,			n, 		k,		alpha,A,m,B,k,beta,C,m)
c		CALL DGEMM('N','N', n_state, n_state, n_contr, alpha, Jx, n_state, Nin, n_contr, beta, Nout, n_state)
	end SUBROUTINE

	SUBROUTINE helperCreateMatrizen(Y, state, M, N)
		integer n_state, n_contr, M_size, N_size,
	1       LDAM, LDBM, LDAN, LDBN

		double precision Y
		double precision state, M, N

		
		parameter(n_state = 13, n_contr = 4, M_size = 169, N_size = 52, 
	1		LDAM = n_state + 1, LDBM = n_state + M_size, 
	2		LDAN = n_state + M_size + 1, LDBN = n_state + M_size + N_size)

		dimension Y(234), state(n_state), M(n_state, n_state),
	1		N(n_state, n_contr)

		state = 0
		M  = 0
		N = 0
		call dcopy(n_state, Y, 1, state, 1)
		
		M = reshape(Y(LDAM:LDBM), (/n_state, n_state/)) 
		N = reshape(Y(LDAN:LDBN), (/n_state, n_contr/))

	end SUBROUTINE


	SUBROUTINE helperCreateVektor(state, M, N, Y)
		integer n_state, n_contr, M_size, N_size, LDAM, LDBM, LDAN,
	1	LDBN

		double precision state, M, N
		double precision Y

		parameter(n_state = 13, n_contr = 4, M_size = 169, N_size = 52, 
	1		LDAM = n_state + 1, LDBM = n_state + M_size, 
	2		LDAN = n_state + M_size + 1, LDBN = n_state + M_size + N_size)			

		dimension Y(234), state(n_state), M(n_state, n_state),
	1		N(n_state, n_contr)

		Y = 0
c       call dcopy(n, x, incx, y, incy) y:=x		
		call dcopy(n_state, state, 1, Y, 1)

		Y(LDAM:LDBM) =	reshape(M, (/M_size/))
		Y(LDAN:LDBN) = reshape(N, (/N_size/))

	end SUBROUTINE
