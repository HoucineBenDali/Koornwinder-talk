# Nonsymmetric Koornwinder polynomials E_eta, n=3, via intertwiners (main.tex normalization).
# T_i: Demazure-Lusztig (eq_Ti-middle, eq_Tn, eq_T0); intertwiner for node i>=1: T_i; node 0: T_0^{-1} z_1.
# Constants fixed formula-free: w = Phi E_eta = beta E_eta + alpha E_{s_i eta}; beta found from one Y_j.
import sys, pickle
n = 3
P = PolynomialRing(QQ, 'q,t,a,b,c,d,z1,z2,z3')
q,t,a,b,c,d,z1,z2,z3 = P.gens()
Z = [z1,z2,z3]
F = P.fraction_field()
t0 = -a*c/q; tn = -b*d

def sub(f, i, val):
    return f.subs({Z[i]: val})

def s(i, f):
    if 1 <= i < n:
        return f.subs({Z[i-1]: Z[i], Z[i]: Z[i-1]})
    if i == n:
        return f.subs({Z[n-1]: 1/Z[n-1]})
    return f.subs({Z[0]: q/Z[0]})

def T(i, f):
    f = F(f)
    if 1 <= i < n:
        x, y = Z[i-1], Z[i]
        return t*f + (t*x - y)/(x - y)*(s(i,f) - f)
    if i == n:
        x = Z[n-1]
        return -b*d*f + (1-b*x)*(1-d*x)/(1-x**2)*(s(i,f) - f)
    x = Z[0]
    return -a*c/q*f + (1-a/x)*(1-c/x)/(1-q/x**2)*(s(0,f) - f)

def tpar(i):
    return t0 if i == 0 else (tn if i == n else t)

def Tinv(i, f):
    return (T(i, f) - (tpar(i)-1)*F(f))/tpar(i)

def Y(i, f):
    # Y_i = T_i ... T_{n-1} T_n T_{n-1} ... T_1 T_0 T_1^{-1} ... T_{i-1}^{-1}  (i is 1-based)
    g = F(f)
    for k in range(1, i):          # rightmost first: T_{i-1}^{-1} ... T_1^{-1}: apply T_{i-1}^{-1} first
        pass
    for k in reversed(range(1, i)):
        g = Tinv(k, g)
    g = T(0, g)
    for k in range(1, n):
        g = T(k, g)
    g = T(n, g)
    for k in reversed(range(i, n)):
        g = T(k, g)
    return g

Kp = PolynomialRing(QQ, 'q,t,a,b,c,d')
K = Kp.fraction_field()
LR = LaurentPolynomialRing(K, 'z1,z2,z3')
def tolaurent(f):
    f = F(f)
    num, den = P(f.numerator()), P(f.denominator())
    zexps = set(tuple(e[6:]) for e in den.exponents())
    assert len(zexps) == 1, "denominator not a z-monomial"
    m = list(zexps)[0]
    dpar = Kp(sum(cf*Kp.monomial(*e[:6]) for e, cf in den.dict().items()))
    res = LR(0)
    zg = LR.gens()
    for e, cf in num.dict().items():
        res += K(cf*Kp.monomial(*e[:6]))/dpar * prod(zg[k]**(e[6+k]-m[k]) for k in range(3))
    return res

def coeffs(f):
    g = tolaurent(f)
    return {tuple(e): K(cf) for e, cf in zip(g.exponents(), g.coefficients())}

def coef(f, e):
    return coeffs(f).get(tuple(e), K(0))

def act(i, eta):
    eta = list(eta)
    if 1 <= i < n:
        eta[i-1], eta[i] = eta[i], eta[i-1]
    elif i == n:
        eta[n-1] = -eta[n-1]
    else:
        eta[0] = -1 - eta[0]
    return tuple(eta)

def step(E, eta, i):
    new = act(i, eta)
    w = T(i, E) if i >= 1 else Tinv(0, F(E)*z1)
    alpha = coef(w, new)
    assert alpha != 0
    for j in range(1, n+1):
        Yw = Y(j, w); YE = Y(j, E)
        yj = coef(YE, eta)              # E monic at eta
        yjp = coef(Yw, new)/alpha
        if yj != yjp:
            beta = coef(Yw - yjp*w, eta)/(yj - yjp)
            En = (w - beta*F(E))/alpha
            return new, En
    raise Exception("no separating Y")

def check(E, eta):
    for j in range(1, n+1):
        YE = Y(j, E); lam = coef(YE, eta)
        assert tolaurent(YE - lam*F(E)) == 0, ("not eigen", j)
    assert coef(E, eta) == 1
    return True


def act_W0(i, eta):
    return act(i, eta)

# E_delta via intertwiners from E_0 = 1
Epaths = {
 (-2,-1,0): [0,1,2,3,2, 0,1,0],
}
# F_kappa = T_i^{-1} F_{s_i kappa}, climbing from delta (qKZ)
Fpaths = {
 (0,2,1): ((-2,-1,0), [2,1,3,2,3,2]),
}
Es = {}
for target, word in Epaths.items():
    eta = (0,0,0); E = F(1)
    for i in word:
        eta, E = step(E, eta, i)
    assert eta == target
    check(E, eta)
    Es[target] = E
    print("E", target, "ok", len(tolaurent(E).exponents()), "monomials"); sys.stdout.flush()

def orbit(lam):
    from itertools import permutations, product
    res = set()
    for p in permutations(lam):
        for sg in product([1,-1], repeat=n):
            res.add(tuple(x*y for x, y in zip(p, sg)))
    return res

results = {}
for kappa, (delta, word) in Fpaths.items():
    Fk = Es[delta]; eta = delta
    # qKZ sanity checks at delta: T_i F = t_i F when s_i fixes delta
    for i in range(1, n+1):
        if act(i, delta) == delta:
            assert tolaurent(T(i, Fk) - tpar(i)*Fk) == 0, ("stabilizer qKZ fails", delta, i)
    for i in word:
        new = act(i, eta)
        Fk = Tinv(i, Fk); eta = new
        # check qKZ in the forward direction: T_i F_new = F_old  (new_i > new_{i+1}, resp. new_n > 0)
    assert eta == kappa
    g = tolaurent(Fk)
    C = dict(zip([tuple(e) for e in g.exponents()], g.coefficients()))
    assert C[kappa] == 1, "normalization"
    lam = tuple(sorted([abs(x) for x in kappa], reverse=True))
    bad = [e for e in C if e in orbit(lam) and e != kappa]
    assert not bad, ("orbit monomials present", bad)
    results[kappa] = g
    print("F", kappa, "ok", len(C), "monomials"); sys.stdout.flush()
save(results, 'F021_results')
