c configDim := Matrix([dimX,  dimDotMatrix,     dimJ,       dimH])
c                      17     13x1              13x17       13x17x17
c                      17   13, 1        13, 17  13, 17, 17              
      SUBROUTINE func(Y, YDOT)
            DOUBLE PRECISION  Y, YDOT
            DOUBLE PRECISION r, q, v, omega, u

            double PRECISION  g, Iges, m, kT, kQ, d, IM


            dimension Y(17), YDOT(13, 1)
            dimension r(3), q(4), v(3), u(4), omega(3)

            dimension Iges(3)

c           Folgende Eigenschaften wurden mit Hilfe unteres GUI_Modeling.fig
c           erstellt und ist im Ordner KonfigurationQuadrocopter zufinden
c           parameter (g = 9.81, Iges = (/0.0093886, 0.0093886, 0.018406/), m = 1.022,  kT  = 1.5e1, kQ  = 3e-01, d = 0.22)
            parameter (Iges = (/0.0093886,0.0093886,0.018406/))
            parameter (g = 9.81, m = 1.022,  kT  = 1.5e1)
            parameter (kQ = 3e-01, d = 0.22, IM = 4.4466e-06)

            YDOT = 0

            r = Y(1:3)
            q = Y(4:7)
            v = Y(8:10)
            omega = Y(11:13)
            u = Y(14:17)  

      YDOT(1,1) = 2 * q(1) * q(3) * v(3) - 2 * q(1) * q(4) * v(2) + 2 * 
     #q(2) * q(3) * v(2) + 2 * q(2) * q(4) * v(3) - 2 * q(3) ** 2 * v(1)
     # - 2 * q(4) ** 2 * v(1) + v(1)
      YDOT(2,1) = -2 * q(1) * q(2) * v(3) + 2 * q(1) * q(4) * v(1) - 2 *
     # q(2) ** 2 * v(2) + 2 * q(2) * q(3) * v(1) + 2 * q(3) * q(4) * v(3
     #) - 2 * q(4) ** 2 * v(2) + v(2)
      YDOT(3,1) = 2 * q(1) * q(2) * v(2) - 2 * q(1) * q(3) * v(1) - 2 * 
     #q(2) ** 2 * v(3) + 2 * q(2) * q(4) * v(1) - 2 * q(3) ** 2 * v(3) +
     # 2 * q(3) * q(4) * v(2) + v(3)
      YDOT(4,1) = -dble(q(1) ** 3) - dble(q(1) * q(2) ** 2) - dble(q(1) 
     #* q(3) ** 2) - dble(q(1) * q(4) ** 2) + dble(q(1)) - dble(q(2) * o
     #mega(1)) / 0.2D1 - dble(q(3) * omega(2)) / 0.2D1 - dble(q(4) * ome
     #ga(3)) / 0.2D1
      YDOT(5,1) = -dble(q(1) ** 2 * q(2)) - dble(q(2) ** 3) - dble(q(2) 
     #* q(3) ** 2) - dble(q(2) * q(4) ** 2) + dble(q(2)) + dble(q(1) * o
     #mega(1)) / 0.2D1 - dble(q(4) * omega(2)) / 0.2D1 + dble(q(3) * ome
     #ga(3)) / 0.2D1
      YDOT(6,1) = -dble(q(1) ** 2 * q(3)) - dble(q(2) ** 2 * q(3)) - dbl
     #e(q(3) ** 3) - dble(q(3) * q(4) ** 2) + dble(q(3)) + dble(q(4) * o
     #mega(1)) / 0.2D1 + dble(q(1) * omega(2)) / 0.2D1 - dble(q(2) * ome
     #ga(3)) / 0.2D1
      YDOT(7,1) = -dble(q(1) ** 2 * q(4)) - dble(q(2) ** 2 * q(4)) - dbl
     #e(q(3) ** 2 * q(4)) - dble(q(4) ** 3) + dble(q(4)) - dble(q(3) * o
     #mega(1)) / 0.2D1 + dble(q(2) * omega(2)) / 0.2D1 + dble(q(1) * ome
     #ga(3)) / 0.2D1
      YDOT(8,1) = 2 * g * q(1) * q(3) - 2 * g * q(2) * q(4) - omega(2) *
     # v(3) + omega(3) * v(2)
      YDOT(9,1) = -2 * g * q(1) * q(2) - 2 * g * q(3) * q(4) + omega(1) 
     #* v(3) - omega(3) * v(1)
      YDOT(10,1) = (2 * g * m * q(2) ** 2 + 2 * g * m * q(3) ** 2 + kT *
     # u(1) ** 2 + kT * u(2) ** 2 + kT * u(3) ** 2 + kT * u(4) ** 2 - m 
     #* omega(1) * v(2) + m * omega(2) * v(1) - m * g) / m
      YDOT(11,1) = -(-d * kT * u(2) ** 2 + d * kT * u(4) ** 2 + IM * ome
     #ga(2) * u(1) - IM * omega(2) * u(2) + IM * omega(2) * u(3) - IM * 
     #omega(2) * u(4) - omega(3) * Iges(2) * omega(2) + omega(2) * Iges(
     #3) * omega(3)) / Iges(1)
      YDOT(12,1) = (-d * kT * u(1) ** 2 + d * kT * u(3) ** 2 + IM * omeg
     #a(1) * u(1) - IM * omega(1) * u(2) + IM * omega(1) * u(3) - IM * o
     #mega(1) * u(4) - omega(3) * Iges(1) * omega(1) + omega(1) * Iges(3
     #) * omega(3)) / Iges(2)
      YDOT(13,1) = -(kQ * u(1) ** 2 - kQ * u(2) ** 2 + kQ * u(3) ** 2 - 
     #kQ * u(4) ** 2 - omega(2) * Iges(1) * omega(1) + omega(1) * Iges(2
     #) * omega(2)) / Iges(3)


            
      end SUBROUTINE func


      SUBROUTINE jac(Y, PD)
            DOUBLE PRECISION  Y, PD
            DOUBLE PRECISION r, q, v, omega, u

            double PRECISION  g, Iges, m, kT, kQ, d, IM


            dimension Y(17), PD(13, 17)
            dimension r(3), q(4), v(3), u(4), omega(3)

            dimension Iges(3)

c           Folgende Eigenschaften wurden mit Hilfe unteres GUI_Modeling.fig
c           erstellt und ist im Ordner KonfigurationQuadrocopter zufinden
c           parameter (g = 9.81, Iges = (/0.0093886, 0.0093886, 0.018406/), m = 1.022,  kT  = 1.5e1, kQ  = 3e-01, d = 0.22)
            parameter (Iges = (/0.0093886,0.0093886,0.018406/))
            parameter (g = 9.81, m = 1.022,  kT  = 1.5e1)
            parameter (kQ = 3e-01, d = 0.22, IM = 4.4466e-06)

            PD = 0

            r = Y(1:3)
            q = Y(4:7)
            v = Y(8:10)
            omega = Y(11:13)
            u = Y(14:17)

      PD(1,4) = 2 * q(3) * v(3) - 2 * q(4) * v(2)
      PD(1,5) = 2 * q(3) * v(2) + 2 * q(4) * v(3)
      PD(1,6) = 2 * q(1) * v(3) + 2 * q(2) * v(2) - 4 * q(3) * v(1)
      PD(1,7) = -2 * q(1) * v(2) + 2 * q(2) * v(3) - 4 * q(4) * v(1)
      PD(1,8) = -2 * q(3) ** 2 - 2 * q(4) ** 2 + 1
      PD(1,9) = -2 * q(1) * q(4) + 2 * q(2) * q(3)
      PD(1,10) = 2 * q(1) * q(3) + 2 * q(2) * q(4)
      PD(2,4) = -2 * q(2) * v(3) + 2 * q(4) * v(1)
      PD(2,5) = -2 * q(1) * v(3) - 4 * q(2) * v(2) + 2 * q(3) * v(1)
      PD(2,6) = 2 * q(2) * v(1) + 2 * q(4) * v(3)
      PD(2,7) = 2 * q(1) * v(1) + 2 * q(3) * v(3) - 4 * q(4) * v(2)
      PD(2,8) = 2 * q(1) * q(4) + 2 * q(2) * q(3)
      PD(2,9) = -2 * q(2) ** 2 - 2 * q(4) ** 2 + 1
      PD(2,10) = -2 * q(1) * q(2) + 2 * q(3) * q(4)
      PD(3,4) = 2 * q(2) * v(2) - 2 * q(3) * v(1)
      PD(3,5) = 2 * q(1) * v(2) - 4 * q(2) * v(3) + 2 * q(4) * v(1)
      PD(3,6) = -2 * q(1) * v(1) - 4 * q(3) * v(3) + 2 * q(4) * v(2)
      PD(3,7) = 2 * q(2) * v(1) + 2 * q(3) * v(2)
      PD(3,8) = -2 * q(1) * q(3) + 2 * q(2) * q(4)
      PD(3,9) = 2 * q(1) * q(2) + 2 * q(3) * q(4)
      PD(3,10) = -2 * q(2) ** 2 - 2 * q(3) ** 2 + 1
      PD(4,4) = -3 * q(1) ** 2 - q(2) ** 2 - q(3) ** 2 - q(4) ** 2 + 1
      PD(4,5) = -dble(2 * q(1) * q(2)) - omega(1) / 0.2D1
      PD(4,6) = -dble(2 * q(1) * q(3)) - omega(2) / 0.2D1
      PD(4,7) = -dble(2 * q(1) * q(4)) - omega(3) / 0.2D1
      PD(4,11) = -dble(q(2)) / 0.2D1
      PD(4,12) = -dble(q(3)) / 0.2D1
      PD(4,13) = -dble(q(4)) / 0.2D1
      PD(5,4) = -dble(2 * q(1) * q(2)) + omega(1) / 0.2D1
      PD(5,5) = -q(1) ** 2 - 3 * q(2) ** 2 - q(3) ** 2 - q(4) ** 2 + 1
      PD(5,6) = -dble(2 * q(2) * q(3)) + omega(3) / 0.2D1
      PD(5,7) = -dble(2 * q(2) * q(4)) - omega(2) / 0.2D1
      PD(5,11) = dble(q(1)) / 0.2D1
      PD(5,12) = -dble(q(4)) / 0.2D1
      PD(5,13) = dble(q(3)) / 0.2D1
      PD(6,4) = -dble(2 * q(1) * q(3)) + omega(2) / 0.2D1
      PD(6,5) = -dble(2 * q(2) * q(3)) - omega(3) / 0.2D1
      PD(6,6) = -q(1) ** 2 - q(2) ** 2 - 3 * q(3) ** 2 - q(4) ** 2 + 1
      PD(6,7) = -dble(2 * q(3) * q(4)) + omega(1) / 0.2D1
      PD(6,11) = dble(q(4)) / 0.2D1
      PD(6,12) = dble(q(1)) / 0.2D1
      PD(6,13) = -dble(q(2)) / 0.2D1
      PD(7,4) = -dble(2 * q(1) * q(4)) + omega(3) / 0.2D1
      PD(7,5) = -dble(2 * q(2) * q(4)) + omega(2) / 0.2D1
      PD(7,6) = -dble(2 * q(3) * q(4)) - omega(1) / 0.2D1
      PD(7,7) = -q(1) ** 2 - q(2) ** 2 - q(3) ** 2 - 3 * q(4) ** 2 + 1
      PD(7,11) = -dble(q(3)) / 0.2D1
      PD(7,12) = dble(q(2)) / 0.2D1
      PD(7,13) = dble(q(1)) / 0.2D1
      PD(8,4) = 2 * g * q(3)
      PD(8,5) = -2 * g * q(4)
      PD(8,6) = 2 * g * q(1)
      PD(8,7) = -2 * g * q(2)
      PD(8,9) = omega(3)
      PD(8,10) = -omega(2)
      PD(8,12) = -v(3)
      PD(8,13) = v(2)
      PD(9,4) = -2 * g * q(2)
      PD(9,5) = -2 * g * q(1)
      PD(9,6) = -2 * g * q(4)
      PD(9,7) = -2 * g * q(3)
      PD(9,8) = -omega(3)
      PD(9,10) = omega(1)
      PD(9,11) = v(3)
      PD(9,13) = -v(1)
      PD(10,5) = 4 * g * q(2)
      PD(10,6) = 4 * g * q(3)
      PD(10,8) = omega(2)
      PD(10,9) = -omega(1)
      PD(10,11) = -v(2)
      PD(10,12) = v(1)
      PD(10,14) = 2 * kT * u(1) / m
      PD(10,15) = 2 * kT * u(2) / m
      PD(10,16) = 2 * kT * u(3) / m
      PD(10,17) = 2 * kT * u(4) / m
      PD(11,12) = -(dble(IM * u(1)) - dble(IM * u(2)) + dble(IM * u(3)) 
     #- dble(IM * u(4)) - Iges(2) * omega(3) + Iges(3) * omega(3)) / Ige
     #s(1)
      PD(11,13) = omega(2) * (Iges(2) - Iges(3)) / Iges(1)
      PD(11,14) = -dble(IM) * omega(2) / Iges(1)
      PD(11,15) = (dble(2 * d * kT * u(2)) + dble(IM) * omega(2)) / Iges
     #(1)
      PD(11,16) = -dble(IM) * omega(2) / Iges(1)
      PD(11,17) = (-dble(2 * d * kT * u(4)) + dble(IM) * omega(2)) / Ige
     #s(1)
      PD(12,11) = (dble(IM * u(1)) - dble(IM * u(2)) + dble(IM * u(3)) -
     # dble(IM * u(4)) - Iges(1) * omega(3) + Iges(3) * omega(3)) / Iges
     #(2)
      PD(12,13) = -omega(1) * (Iges(1) - Iges(3)) / Iges(2)
      PD(12,14) = (-dble(2 * d * kT * u(1)) + dble(IM) * omega(1)) / Ige
     #s(2)
      PD(12,15) = -dble(IM) * omega(1) / Iges(2)
      PD(12,16) = (dble(2 * d * kT * u(3)) + dble(IM) * omega(1)) / Iges
     #(2)
      PD(12,17) = -dble(IM) * omega(1) / Iges(2)
      PD(13,11) = omega(2) * (Iges(1) - Iges(2)) / Iges(3)
      PD(13,12) = omega(1) * (Iges(1) - Iges(2)) / Iges(3)
      PD(13,14) = -0.2D1 * kQ * dble(u(1)) / Iges(3)
      PD(13,15) = 0.2D1 * kQ * dble(u(2)) / Iges(3)
      PD(13,16) = -0.2D1 * kQ * dble(u(3)) / Iges(3)
      PD(13,17) = 0.2D1 * kQ * dble(u(4)) / Iges(3)

      

      end SUBROUTINE

      SUBROUTINE hess(Y, PDD)
            DOUBLE PRECISION  Y, PDD
            DOUBLE PRECISION r, q, v, omega, u

            double PRECISION  g, Iges, m, kT, kQ, d, IM

            dimension Y(17), PDD(13, 17, 17)
            dimension r(3), q(4), v(3), u(4), omega(3)

            dimension Iges(3)

c           Folgende Eigenschaften wurden mit Hilfe unteres GUI_Modeling.fig
c           erstellt und ist im Ordner KonfigurationQuadrocopter zufinden
c           parameter (g = 9.81, Iges = (/0.0093886, 0.0093886, 0.018406/), m = 1.022,  kT  = 1.5e1, kQ  = 3e-01, d = 0.22)
            parameter (Iges = (/0.0093886,0.0093886,0.018406/))
            parameter (g = 9.81, m = 1.022,  kT  = 1.5e1)
            parameter (kQ = 3e-01, d = 0.22, IM = 4.4466e-06)

            PDD = 0

            r = Y(1:3)
            q = Y(4:7)
            v = Y(8:10)
            omega = Y(11:13)
            u = Y(14:17)

      PDD(1,4,6) = 2 * v(3)
      PDD(1,4,7) = -2 * v(2)
      PDD(1,4,9) = -2 * q(4)
      PDD(1,4,10) = 2 * q(3)
      PDD(1,5,6) = 2 * v(2)
      PDD(1,5,7) = 2 * v(3)
      PDD(1,5,9) = 2 * q(3)
      PDD(1,5,10) = 2 * q(4)
      PDD(1,6,4) = 2 * v(3)
      PDD(1,6,5) = 2 * v(2)
      PDD(1,6,6) = -4 * v(1)
      PDD(1,6,8) = -4 * q(3)
      PDD(1,6,9) = 2 * q(2)
      PDD(1,6,10) = 2 * q(1)
      PDD(1,7,4) = -2 * v(2)
      PDD(1,7,5) = 2 * v(3)
      PDD(1,7,7) = -4 * v(1)
      PDD(1,7,8) = -4 * q(4)
      PDD(1,7,9) = -2 * q(1)
      PDD(1,7,10) = 2 * q(2)
      PDD(1,8,6) = -4 * q(3)
      PDD(1,8,7) = -4 * q(4)
      PDD(1,9,4) = -2 * q(4)
      PDD(1,9,5) = 2 * q(3)
      PDD(1,9,6) = 2 * q(2)
      PDD(1,9,7) = -2 * q(1)
      PDD(1,10,4) = 2 * q(3)
      PDD(1,10,5) = 2 * q(4)
      PDD(1,10,6) = 2 * q(1)
      PDD(1,10,7) = 2 * q(2)
      PDD(2,4,5) = -2 * v(3)
      PDD(2,4,7) = 2 * v(1)
      PDD(2,4,8) = 2 * q(4)
      PDD(2,4,10) = -2 * q(2)
      PDD(2,5,4) = -2 * v(3)
      PDD(2,5,5) = -4 * v(2)
      PDD(2,5,6) = 2 * v(1)
      PDD(2,5,8) = 2 * q(3)
      PDD(2,5,9) = -4 * q(2)
      PDD(2,5,10) = -2 * q(1)
      PDD(2,6,5) = 2 * v(1)
      PDD(2,6,7) = 2 * v(3)
      PDD(2,6,8) = 2 * q(2)
      PDD(2,6,10) = 2 * q(4)
      PDD(2,7,4) = 2 * v(1)
      PDD(2,7,6) = 2 * v(3)
      PDD(2,7,7) = -4 * v(2)
      PDD(2,7,8) = 2 * q(1)
      PDD(2,7,9) = -4 * q(4)
      PDD(2,7,10) = 2 * q(3)
      PDD(2,8,4) = 2 * q(4)
      PDD(2,8,5) = 2 * q(3)
      PDD(2,8,6) = 2 * q(2)
      PDD(2,8,7) = 2 * q(1)
      PDD(2,9,5) = -4 * q(2)
      PDD(2,9,7) = -4 * q(4)
      PDD(2,10,4) = -2 * q(2)
      PDD(2,10,5) = -2 * q(1)
      PDD(2,10,6) = 2 * q(4)
      PDD(2,10,7) = 2 * q(3)
      PDD(3,4,5) = 2 * v(2)
      PDD(3,4,6) = -2 * v(1)
      PDD(3,4,8) = -2 * q(3)
      PDD(3,4,9) = 2 * q(2)
      PDD(3,5,4) = 2 * v(2)
      PDD(3,5,5) = -4 * v(3)
      PDD(3,5,7) = 2 * v(1)
      PDD(3,5,8) = 2 * q(4)
      PDD(3,5,9) = 2 * q(1)
      PDD(3,5,10) = -4 * q(2)
      PDD(3,6,4) = -2 * v(1)
      PDD(3,6,6) = -4 * v(3)
      PDD(3,6,7) = 2 * v(2)
      PDD(3,6,8) = -2 * q(1)
      PDD(3,6,9) = 2 * q(4)
      PDD(3,6,10) = -4 * q(3)
      PDD(3,7,5) = 2 * v(1)
      PDD(3,7,6) = 2 * v(2)
      PDD(3,7,8) = 2 * q(2)
      PDD(3,7,9) = 2 * q(3)
      PDD(3,8,4) = -2 * q(3)
      PDD(3,8,5) = 2 * q(4)
      PDD(3,8,6) = -2 * q(1)
      PDD(3,8,7) = 2 * q(2)
      PDD(3,9,4) = 2 * q(2)
      PDD(3,9,5) = 2 * q(1)
      PDD(3,9,6) = 2 * q(4)
      PDD(3,9,7) = 2 * q(3)
      PDD(3,10,5) = -4 * q(2)
      PDD(3,10,6) = -4 * q(3)
      PDD(4,4,4) = -6 * q(1)
      PDD(4,4,5) = -2 * q(2)
      PDD(4,4,6) = -2 * q(3)
      PDD(4,4,7) = -2 * q(4)
      PDD(4,5,4) = -2 * q(2)
      PDD(4,5,5) = -2 * q(1)
      PDD(4,5,11) = -0.1D1 / 0.2D1
      PDD(4,6,4) = -2 * q(3)
      PDD(4,6,6) = -2 * q(1)
      PDD(4,6,12) = -0.1D1 / 0.2D1
      PDD(4,7,4) = -2 * q(4)
      PDD(4,7,7) = -2 * q(1)
      PDD(4,7,13) = -0.1D1 / 0.2D1
      PDD(4,11,5) = -0.1D1 / 0.2D1
      PDD(4,12,6) = -0.1D1 / 0.2D1
      PDD(4,13,7) = -0.1D1 / 0.2D1
      PDD(5,4,4) = -2 * q(2)
      PDD(5,4,5) = -2 * q(1)
      PDD(5,4,11) = 0.1D1 / 0.2D1
      PDD(5,5,4) = -2 * q(1)
      PDD(5,5,5) = -6 * q(2)
      PDD(5,5,6) = -2 * q(3)
      PDD(5,5,7) = -2 * q(4)
      PDD(5,6,5) = -2 * q(3)
      PDD(5,6,6) = -2 * q(2)
      PDD(5,6,13) = 0.1D1 / 0.2D1
      PDD(5,7,5) = -2 * q(4)
      PDD(5,7,7) = -2 * q(2)
      PDD(5,7,12) = -0.1D1 / 0.2D1
      PDD(5,11,4) = 0.1D1 / 0.2D1
      PDD(5,12,7) = -0.1D1 / 0.2D1
      PDD(5,13,6) = 0.1D1 / 0.2D1
      PDD(6,4,4) = -2 * q(3)
      PDD(6,4,6) = -2 * q(1)
      PDD(6,4,12) = 0.1D1 / 0.2D1
      PDD(6,5,5) = -2 * q(3)
      PDD(6,5,6) = -2 * q(2)
      PDD(6,5,13) = -0.1D1 / 0.2D1
      PDD(6,6,4) = -2 * q(1)
      PDD(6,6,5) = -2 * q(2)
      PDD(6,6,6) = -6 * q(3)
      PDD(6,6,7) = -2 * q(4)
      PDD(6,7,6) = -2 * q(4)
      PDD(6,7,7) = -2 * q(3)
      PDD(6,7,11) = 0.1D1 / 0.2D1
      PDD(6,11,7) = 0.1D1 / 0.2D1
      PDD(6,12,4) = 0.1D1 / 0.2D1
      PDD(6,13,5) = -0.1D1 / 0.2D1
      PDD(7,4,4) = -2 * q(4)
      PDD(7,4,7) = -2 * q(1)
      PDD(7,4,13) = 0.1D1 / 0.2D1
      PDD(7,5,5) = -2 * q(4)
      PDD(7,5,7) = -2 * q(2)
      PDD(7,5,12) = 0.1D1 / 0.2D1
      PDD(7,6,6) = -2 * q(4)
      PDD(7,6,7) = -2 * q(3)
      PDD(7,6,11) = -0.1D1 / 0.2D1
      PDD(7,7,4) = -2 * q(1)
      PDD(7,7,5) = -2 * q(2)
      PDD(7,7,6) = -2 * q(3)
      PDD(7,7,7) = -6 * q(4)
      PDD(7,11,6) = -0.1D1 / 0.2D1
      PDD(7,12,5) = 0.1D1 / 0.2D1
      PDD(7,13,4) = 0.1D1 / 0.2D1
      PDD(8,4,6) = 2 * g
      PDD(8,5,7) = -2 * g
      PDD(8,6,4) = 2 * g
      PDD(8,7,5) = -2 * g
      PDD(8,9,13) = 1
      PDD(8,10,12) = -1
      PDD(8,12,10) = -1
      PDD(8,13,9) = 1
      PDD(9,4,5) = -2 * g
      PDD(9,5,4) = -2 * g
      PDD(9,6,7) = -2 * g
      PDD(9,7,6) = -2 * g
      PDD(9,8,13) = -1
      PDD(9,10,11) = 1
      PDD(9,11,10) = 1
      PDD(9,13,8) = -1
      PDD(10,5,5) = 4 * g
      PDD(10,6,6) = 4 * g
      PDD(10,8,12) = 1
      PDD(10,9,11) = -1
      PDD(10,11,9) = -1
      PDD(10,12,8) = 1
      PDD(10,14,14) = 2 * kT / m
      PDD(10,15,15) = 2 * kT / m
      PDD(10,16,16) = 2 * kT / m
      PDD(10,17,17) = 2 * kT / m
      PDD(11,12,13) = -(-dble(Iges(2)) + dble(Iges(3))) / dble(Iges(1))
      PDD(11,12,14) = -IM / Iges(1)
      PDD(11,12,15) = IM / Iges(1)
      PDD(11,12,16) = -IM / Iges(1)
      PDD(11,12,17) = IM / Iges(1)
      PDD(11,13,12) = -(-dble(Iges(2)) + dble(Iges(3))) / dble(Iges(1))
      PDD(11,14,12) = -IM / Iges(1)
      PDD(11,15,12) = IM / Iges(1)
      PDD(11,15,15) = 2 * d * kT / Iges(1)
      PDD(11,16,12) = -IM / Iges(1)
      PDD(11,17,12) = IM / Iges(1)
      PDD(11,17,17) = -2 * d * kT / Iges(1)
      PDD(12,11,13) = (-Iges(1) + Iges(3)) / Iges(2)
      PDD(12,11,14) = IM / Iges(2)
      PDD(12,11,15) = -IM / Iges(2)
      PDD(12,11,16) = IM / Iges(2)
      PDD(12,11,17) = -IM / Iges(2)
      PDD(12,13,11) = (-Iges(1) + Iges(3)) / Iges(2)
      PDD(12,14,11) = IM / Iges(2)
      PDD(12,14,14) = -2 * d * kT / Iges(2)
      PDD(12,15,11) = -IM / Iges(2)
      PDD(12,16,11) = IM / Iges(2)
      PDD(12,16,16) = 2 * d * kT / Iges(2)
      PDD(12,17,11) = -IM / Iges(2)
      PDD(13,11,12) = -(-Iges(1) + Iges(2)) / Iges(3)
      PDD(13,12,11) = -(-Iges(1) + Iges(2)) / Iges(3)
      PDD(13,14,14) = -2 * kQ / Iges(3)
      PDD(13,15,15) = 2 * kQ / Iges(3)
      PDD(13,16,16) = -2 * kQ / Iges(3)
      PDD(13,17,17) = 2 * kQ / Iges(3)



      end SUBROUTINE

      SUBROUTINE getTestData(DMat)
            DOUBLE PRECISION DMat
            dimension DMat(51, 1)
            DMat = 0

      DMat(1,1) = 0.545886208890829838D0
      DMat(2,1) = 0.943169839655478737D0
      DMat(3,1) = 0.321473069869989025D0
      DMat(4,1) = 0.806466803793333575D0
      DMat(5,1) = 0.601398755252809170D0
      DMat(6,1) = 0.789620465123038406D0
      DMat(7,1) = 0.799185035044315595D0
      DMat(8,1) = 0.495647655782109897D-1
      DMat(9,1) = 0.283198633840945568D0
      DMat(10,1) = 0.653456823979964518D0
      DMat(11,1) = 0.489655345328471103D0
      DMat(12,1) = 0.972852237010802923D0
      DMat(13,1) = 0.748489909977081980D0
      DMat(14,1) = 0.567841149738807283D0
      DMat(15,1) = 0.298964163457197118D0
      DMat(16,1) = 0.256109781061838127D0
      DMat(17,1) = 0.886563794411546979D0
      DMat(18,1) = 0.446800862782800512D0
      DMat(19,1) = 0.815987253296207360D0
      DMat(20,1) = 0.983373017446377284D-1
      DMat(21,1) = 0.859593453475272384D0
      DMat(22,1) = 0.276290081611750837D-1
      DMat(23,1) = 0.899156434256364379D0
      DMat(24,1) = 0.899935500526175547D0
      DMat(25,1) = 0.524106009617362978D0
      DMat(26,1) = 0.120199512647134799D0
      DMat(27,1) = 0.177794091822278899D0
      DMat(28,1) = 0.706107594508614378D0
      DMat(29,1) = 0.831359756387049287D0
      DMat(30,1) = 0.348338205677289903D-1
      DMat(31,1) = 0.757838710385932535D0
      DMat(32,1) = 0.957112178001164104D0
      DMat(33,1) = 0.342870514996847620D0
      DMat(34,1) = 0.638243709424705874D0
      DMat(35,1) = 0.343005811970325092D0
      DMat(36,1) = 0.216471396182540854D0
      DMat(37,1) = 0.786200627030544497D0
      DMat(38,1) = 0.723089959545471506D0
      DMat(39,1) = 0.278838904787457209D0
      DMat(40,1) = 0.582431428704480059D0
      DMat(41,1) = 0.421005781368280152D0
      DMat(42,1) = 0.920687067543677351D-1
      DMat(43,1) = 0.240273250067485344D-1
      DMat(44,1) = 0.491145804011365428D0
      DMat(45,1) = 0.278267020151879096D0
      DMat(46,1) = 0.339757164151064717D0
      DMat(47,1) = 0.287349609515075932D0
      DMat(48,1) = 0.170903235249478369D0
      DMat(49,1) = 0.399263314281429471D0
      DMat(50,1) = 0.697649568457898273D0
      DMat(51,1) = 0.203676437129014953D0


      end SUBROUTINE


