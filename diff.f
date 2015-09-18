	SUBROUTINE diff_nD(vec_old, F, n, PD)
		
c 		Uebergabeparameter
		external F
		double precision vec_old
		integer n
		double precision PD

c 		Konstanten
		integer n_state, n_var, n_contr, inc
		double precision alpha, beta, eps

c 		Zaehlervariablen		
		integer i
c 		in - function calls
		double precision vec, func_p, func_n

c 		sonstige Parameter
		INTEGER  NEQ
		double precision T


		parameter(n_state = 13, n_var = 17, n_contr = 4, 
	1		inc = 1, eps = 1e-2, alpha = -1.0, beta = 1/(2*eps))

		dimension vec(n), func_p(n), 
	1		func_n(n), PD(n, n), vec_old(n)

		PD = 0

		NEQ = 0
		T = 0

		func_p = 0

		do 10 i = 1, n
			vec = vec_old
			vec(i) = vec_old(i) + eps
			call F(NEQ, T, vec, func_p)
			vec(i) = vec_old(i) - eps
			call F(NEQ, T, vec, func_n)

			call daxpy(n, alpha, func_n, inc, func_p, inc)
			call daxpy(n, beta, func_p, inc, PD(:, i), inc)
10		continue		

	end SUBROUTINE