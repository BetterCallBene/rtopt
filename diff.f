	SUBROUTINE diff_nD1(NEQ, vec_old, u, diffFunc, numPD)
		
c 		Uebergabeparameter
		external diffFunc
		integer NEQ
		double precision vec_old, numPD, u

c 		Konstanten
		integer n_state, n_var, n_contr, inc
		double precision alpha, beta, eps

c 		Zaehlervariablen		
		integer i
c 		in - function calls
		double precision vec, func_p, func_n

c 		sonstige Parameter
		double precision T


		parameter(n_state = 13, n_var = 17, n_contr = 4,
	1		inc = 1, eps = 1e-2, alpha = -1.0, beta = 1/(2*eps))

		dimension vec(NEQ), func_p(NEQ), u(n_contr),
	1		func_n(NEQ), vec_old(NEQ),
	2       numPD(NEQ, NEQ)

        T = 0
       	numPD =0
		func_p = 0

		do 10 i = 1, NEQ
			vec = vec_old
			vec(i) = vec_old(i) + eps
			call diffFunc(NEQ, T, vec, u, func_p)

			vec(i) = vec_old(i) - eps
			call diffFunc(NEQ, T, vec, u, func_n)

			call daxpy(NEQ, alpha, func_n, inc, func_p, inc)
			call daxpy(NEQ, beta, func_p, inc, numPD(:, i), inc)
10		continue

	end SUBROUTINE
