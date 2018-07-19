import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
import numpy as np

Print = PETSc.Sys.Print

class ShellMatrix(object):

    def __init__(self, m, n, KL11, KL22, L21, L22, M11, M12):
        self.m, self.n = m, n
        scalar = PETSc.ScalarType
        self.KL11 = KL11
        self.KL22 = KL22
        self.L21 = L21
        self.L22 = L22
        self.M11 = M11
        self.M12 = M12
        self.workvec, _ = L22.createVecs()
        self.workvec.set(0)

    def mult(self, A, x, y):
        """
        Second version: isolating \Phi_2
        """
        w1 = self.M11 * x
        w2 = self.L21 * x
        self.KL22.solve(w2, self.workvec)
        w4 = w1 + self.M12 * self.workvec
        self.KL11.solve(w4, y)
       

def construct_operator(m, n, KL11, KL22, L21, L22, M11, M12):
    # Create shell matrix
    context = ShellMatrix(m,n,KL11, KL22, L21, L22, M11, M12)
    A = PETSc.Mat().createPython([m,n], context)
    A.setUp()
    return A

def solve_eigensystem(A, problem_type=SLEPc.EPS.ProblemType.NHEP):
    # Create the result vectors
    xr, xi = A.createVecs()

    # Setup the eigensolver
    E = SLEPc.EPS().create()
    E.setOperators(A,None)
    E.setDimensions(3,PETSc.DECIDE)
    E.setProblemType( problem_type )
    E.setFromOptions()

    # Solve the eigensystem
    E.solve()
    Print("")
    its = E.getIterationNumber()
    Print("Number of iterations of the method: %i" % its)
    sol_type = E.getType()
    Print("Solution method: %s" % sol_type)
    nev, ncv, mpd = E.getDimensions()
    Print("Number of requested eigenvalues: %i" % nev)
    tol, maxit = E.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))
    nconv = E.getConverged()
    Print("Number of converged eigenpairs: %d" % nconv)
    if nconv > 0:
        Print("")
        Print("        k          ||Ax-kx||/||kx|| ")
        Print("----------------- ------------------")
        for i in range(nconv):
            k = E.getEigenpair(i, xr, xi)
            error = E.computeError(i)
            if k.imag != 0.0:
              Print(" %9f%+9f j  %12g" % (k.real, k.imag, error))
            else:
              Print(" %12f       %12g" % (k.real, error))
        Print("")


def main():
    opts = PETSc.Options()
    # load from file
    viewer = PETSc.Viewer().createBinary('../matrices/ringhals1.petsc', 'r')

    L11 = PETSc.Mat().load(viewer)
    L22 = PETSc.Mat().load(viewer)

    L21 = PETSc.Vec().load(viewer)

    M11 = PETSc.Vec().load(viewer)
    M12 = PETSc.Vec().load(viewer)

    # create linear solver
    KL11 = PETSc.KSP()
    KL11.create(PETSc.COMM_WORLD)
    # use conjugate gradients
    KL11.setType('cg')
    # and incomplete Cholesky
    KL11.getPC().setType('none')
    # obtain sol & rhs vectors
    KL11.setOperators(L11)
    KL11.setFromOptions()

    ((lr,gr),(lc,gc)) = L11.getSizes()

    # create linear solver
    KL22 = PETSc.KSP()
    KL22.create(PETSc.COMM_WORLD)
    # use conjugate gradients
    KL22.setType('cg')
    # and incomplete Cholesky
    KL22.getPC().setType('none')
    # obtain sol & rhs vectors
    KL22.setOperators(L22)
    KL22.setFromOptions()

    Print("gr={:}\n".format(gr))
    Print("gc={:}\n".format(gc))
    Print("lr={:}\n".format(lr))
    Print("lc={:}\n".format(lc))
    Print("Standard Non-Symmetric Eigenvalue Problem (matrix-free)")
    A = construct_operator(gr,gc, KL11, KL22, L21, L22, M11, M12)
    solve_eigensystem(A)

if __name__ == '__main__':
    main()
