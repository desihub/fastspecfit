#import numba
#@numba.jit(nopython=True)
def nnls(A, b, maxiter=None, eps=1e-7, v=False):
    # https://en.wikipedia.org/wiki/Non-negative_least_squares#Algorithms
    # not sure what eps should be set to
    m, n = A.shape
    if b.shape[0] != m:
        raise ValueError()#f"Shape mismatch: {A.shape} {b.shape}. Expected: (m, n) (m) ")
    if maxiter is None:
        # this is what scipy does when maxiter is not specified
        maxiter = 3 * n
    # init
    P = np.zeros(n).astype(bool)
    R = np.ones(n).astype(bool)
    x = np.zeros(n)
    s = np.zeros(n)
    # compute these once ahead of time
    ATA = A.T.dot(A)
    ATb = A.T.dot(b)
    # main loop
    w = ATb - ATA.dot(x)
    j = np.argmax(R * w)
    c = 0 # iteration count
    # while R != {} and max(w[R]) > eps
    while np.any(R) and w[j] > eps:
        #if v: print(f"{c=}", f"{P=}\n {R=}\n {w=}\n {j=}")
        # add j to P, remove j from R
        P[j], R[j] = True, False
        s[P] = np.linalg.inv(ATA[P][:, P]).dot(ATb[P])
        s[R] = 0
        d = 0 # inner loop iteration count, for debugging
        #if v: print(f"{c=}", f"{P=}\n {R=}\n {s=}")
        # while P != {} and min(s[P]) < eps
        # make sure P is not empty before checking min s[P]
        while np.any(P) and np.min(s[P]) < eps:
            i = P & (s < eps)
            #if v: print(f" {d=}", f"  {P=}\n  {i=}\n  {s=}\n  {x=}")
            a = np.min(x[i] / (x[i] - s[i]))
            x = x + a * (s - x)
            j = P & (x < eps)
            R[j], P[j] = True, False
            s[P] = np.linalg.inv(ATA[P][:, P]).dot(ATb[P])
            s[R] = 0
            d += 1
        x[:] = s
        # w = A.T.dot(b - A.dot(x))
        w = ATb - ATA.dot(x)
        j = np.argmax(R * w)
        #if v: print(f"{c=}", f"{P=}\n {R=}\n {w=}\n {j=}")
        c += 1
        if c >= maxiter:
            break
    res = np.linalg.norm(A.dot(x) - b)
    return x, res

