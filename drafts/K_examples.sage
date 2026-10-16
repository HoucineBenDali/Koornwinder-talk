# Symmetric Koornwinder K_lambda, n=2, as sum of F_mu over W_0-orbit (F_delta=E_delta, F_mu = T_i^{-1} F_{s_i mu}).
import sys
n = 2
P = PolynomialRing(QQ, 'q,t,a,b,c,d,z1,z2')
q,t,a,b,c,d,z1,z2 = P.gens(); Z=[z1,z2]
F = P.fraction_field()
t0 = -a*c/q; tn = -b*d
def s(i, f):
    if 1 <= i < n: return f.subs({Z[i-1]: Z[i], Z[i]: Z[i-1]})
    if i == n: return f.subs({Z[n-1]: 1/Z[n-1]})
    return f.subs({Z[0]: q/Z[0]})
def T(i, f):
    f = F(f)
    if 1 <= i < n:
        x, y = Z[i-1], Z[i]; return t*f + (t*x - y)/(x - y)*(s(i,f) - f)
    if i == n:
        x = Z[n-1]; return -b*d*f + (1-b*x)*(1-d*x)/(1-x**2)*(s(i,f) - f)
    x = Z[0]; return -a*c/q*f + (1-a/x)*(1-c/x)/(1-q/x**2)*(s(0,f) - f)
def tpar(i): return t0 if i == 0 else (tn if i == n else t)
def Tinv(i, f): return (T(i, f) - (tpar(i)-1)*F(f))/tpar(i)
def Y(i, f):
    g = F(f)
    for k in reversed(range(1, i)): g = Tinv(k, g)
    g = T(0, g)
    for k in range(1, n): g = T(k, g)
    g = T(n, g)
    for k in reversed(range(i, n)): g = T(k, g)
    return g
Kp = PolynomialRing(QQ, 'q,t,a,b,c,d'); K = Kp.fraction_field()
LR = LaurentPolynomialRing(K, 'z1,z2')
def tolaurent(f):
    f = F(f); num, den = P(f.numerator()), P(f.denominator())
    zexps = set(tuple(e[6:]) for e in den.exponents()); assert len(zexps) == 1
    m = list(zexps)[0]
    dpar = Kp(sum(cf*Kp.monomial(*e[:6]) for e, cf in den.dict().items()))
    zg = LR.gens(); res = LR(0)
    for e, cf in num.dict().items():
        res += K(cf*Kp.monomial(*e[:6]))/dpar * prod(zg[k]**(e[6+k]-m[k]) for k in range(n))
    return res
def coef(f, e):
    g = tolaurent(f); C = {tuple(x): K(y) for x, y in zip(g.exponents(), g.coefficients())}
    return C.get(tuple(e), K(0))
def act(i, eta):
    eta = list(eta)
    if 1 <= i < n: eta[i-1], eta[i] = eta[i], eta[i-1]
    elif i == n: eta[n-1] = -eta[n-1]
    else: eta[0] = -1 - eta[0]
    return tuple(eta)
def step(E, eta, i):
    new = act(i, eta)
    w = T(i, E) if i >= 1 else Tinv(0, F(E)*z1)
    alpha = coef(w, new); assert alpha != 0
    for j in range(1, n+1):
        Yw = Y(j, w); YE = Y(j, E)
        yj = coef(YE, eta); yjp = coef(Yw, new)/alpha
        if yj != yjp:
            beta = coef(Yw - yjp*w, eta)/(yj - yjp)
            return new, (w - beta*F(E))/alpha
    raise Exception
def check(E, eta):
    for j in range(1, n+1):
        YE = Y(j, E); lam = coef(YE, eta)
        assert tolaurent(YE - lam*F(E)) == 0
cases = {
 (1,0): ([0], {(-1,0):[], (0,-1):[1], (0,1):[1,2], (1,0):[1,2,1]}),
 (1,1): ([0,1,0], {(-1,-1):[], (-1,1):[2], (1,-1):[2,1], (1,1):[2,1,2]}),
}
for lam,(Ew,Fw) in cases.items():
    delta = tuple(-x for x in lam)
    eta=(0,0); E=F(1)
    for i in Ew: eta,E = step(E,eta,i)
    assert eta==delta, eta
    check(E,eta)
    Ksum = F(0)
    for mu,word in Fw.items():
        Fm=E; e=delta
        for i in word: Fm=Tinv(i,Fm); e=act(i,e)
        assert e==mu,(e,mu)
        Ksum += Fm
    g = tolaurent(Ksum); g = g/ g.monomial_coefficient(LR.monomial(*lam))
    C = {tuple(e):cf for e,cf in zip(g.exponents(), g.coefficients())}
    for e,cf in C.items():
        for e2 in [(e[1],e[0]),(e[0],-e[1])]:
            assert C.get(e2)==cf, ("not symmetric", e, e2)
    print(lam)
    for e,cf in zip(g.exponents(), g.coefficients()):
        print("  ", e, factor(cf))
    sys.stdout.flush()
