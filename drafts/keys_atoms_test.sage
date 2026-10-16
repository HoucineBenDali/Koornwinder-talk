import sage.combinat.sf.ns_macdonald as NS
from itertools import product
n = 3
Kp.<q,t> = QQ['q,t']; K = Kp.fraction_field(); q,t = K.gens()
R = PolynomialRing(K, ['x%d'%i for i in range(n)], implementation='generic'); X = R.gens()
Ksp.<s> = QQ['s']; Ks = Ksp.fraction_field(); s = Ks.gen()

def HHL_E(alpha):
    p = NS.E(list(alpha), q, t)
    return R(sage_eval(str(p), locals=dict({'x%d'%i: X[i] for i in range(n)}, q=q, t=t)))

def swap(f, i):  # s_i, 0-based i <-> i+1
    return f.subs({X[i]: X[i+1], X[i+1]: X[i]})
def exdiv(num, i):
    qq, rr = num.quo_rem(X[i]-X[i+1]); assert rr == 0; return qq
def T(f, i):
    return t*f + exdiv((t*X[i]-X[i+1])*(swap(f,i)-f), i)
def pi(f, i):
    return exdiv(X[i]*f - X[i+1]*swap(f,i), i)

def comps(maxsum):
    return [a for a in product(range(maxsum+1), repeat=n) if 0 < sum(a) <= maxsum]

def sort_desc(a): return tuple(sorted(a, reverse=True))

def climb(alpha, base, op):
    """build object for alpha from dominant (decreasing) rearrangement via op_i when a_i > a_{i+1} -> s_i a"""
    lam = sort_desc(alpha)
    # BFS over rearrangements
    known = {lam: base(lam)}; frontier = [lam]
    while frontier:
        new = []
        for b in frontier:
            for i in range(n-1):
                if b[i] > b[i+1]:
                    c = list(b); c[i], c[i+1] = c[i+1], c[i]; c = tuple(c)
                    if c not in known:
                        known[c] = op(known[b], i); new.append(c)
        frontier = new
    return known[alpha]

key  = lambda a: climb(a, lambda l: prod(X[i]**l[i] for i in range(n)), lambda f,i: pi(f,i))
atom = lambda a: climb(a, lambda l: prod(X[i]**l[i] for i in range(n)), lambda f,i: pi(f,i)-f)
F_A  = lambda a: climb(a, HHL_E, lambda f,i: T(f,i))

def spec(f, mode):
    """mode '0': q=t=0 ; 'inf': q=t=1/s, s->0 ; returns polynomial over QQ or None if pole"""
    out = 0
    for cf, m in zip(f.coefficients(), f.monomials()):
        cf = K(cf)
        if mode == '0':
            num, den = cf.numerator(), cf.denominator()
            if den(0,0) == 0: return None
            v = num(0,0)/den(0,0)
        else:
            g = Ks(cf.numerator()(1/s,1/s)) / Ks(cf.denominator()(1/s,1/s))
            N, D = g.numerator(), g.denominator()
            vn = min(N.exponents()); vd = min(D.exponents())
            if vn < vd: return None
            v = (N[vn]/D[vd]) if vn == vd else 0
        out += v * m
    return out

Rq = PolynomialRing(QQ, ['x%d'%i for i in range(n)]); Y = Rq.gens()
def toQ(f): return Rq(f) if f is not None else None
def rev_vars(f): return f.subs({Y[i]: Y[n-1-i] for i in range(n)}) if f is not None else None

alphas = comps(3)
tests = {}
for a in alphas:
    E = HHL_E(a); Fa = F_A(a)
    kap = {b: toQ(key(b)) for b in [a, a[::-1]]}; at = {b: toQ(atom(b)) for b in [a, a[::-1]]}
    for name, P in [('E', E), ('F', Fa)]:
        for mode in ['0', 'inf']:
            val = toQ(spec(P, mode))
            for vr in [False, True]:
                v = rev_vars(val) if vr else val
                for idx in ['a', 'rev a']:
                    b = a if idx == 'a' else a[::-1]
                    for tname, tgt in [('key', kap[b]), ('atom', at[b])]:
                        lab = (name, 'q=t=0' if mode=='0' else 'q=t->inf', 'x reversed' if vr else 'x', tname + '_' + idx)
                        tests.setdefault(lab, []).append(v == tgt)
print(len(alphas), "compositions tested")
for lab, res in sorted(tests.items()):
    if all(res): print("ALWAYS:", lab)
for lab, res in sorted(tests.items()):
    if any(res) and not all(res): print("sometimes (%d/%d):" % (sum(res), len(res)), lab)
# sample
a=(0,1,2)
print("E_(0,1,2) at q=t=0:", toQ(spec(HHL_E(a),'0')), "| F_(0,1,2) at q=t=0:", toQ(spec(F_A(a),'0')))
print("key_(0,1,2):", toQ(key(a)), "| atom_(0,1,2):", toQ(atom(a)))
