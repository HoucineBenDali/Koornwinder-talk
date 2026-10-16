S.<Q,T,s,p,u,v> = PolynomialRing(QQ, order='deglex')
r2 = S(load('N020_sym'))
def mono_tex(m):
    e = dict(zip(S.gens(), m.exponents()[0]))
    out = ''
    for var, sym in [(Q,'q'),(T,'t')]:
        k = e[var]
        if k == 1: out += sym
        elif k > 1: out += '%s^{%d}' % (sym, k)
    for var, one, two in [(p,'ac','a^2c^2'),(v,'bd','b^2d^2')]:
        k = e[var]
        if k == 1: out += one
        elif k == 2: out += two
        elif k > 2: out += '(%s)^{%d}' % (one, k)
    for var, sym in [(s,'(a+c)'),(u,'(b+d)')]:
        k = e[var]
        if k == 1: out += sym
        elif k > 1: out += sym + '^{%d}' % k
    return out
def mono_sage(m):
    e = m.exponents()[0]
    names = ['q','t','(a+c)','(a*c)','(b+d)','(b*d)']
    return '*'.join('%s^%d' % (nm, k) for nm, k in zip(names, e) if k) or '1'
def poly_tex(P):
    terms = sorted(zip(P.coefficients(), P.monomials()), key=lambda x: x[1], reverse=True)
    out = ''
    for i, (cf, m) in enumerate(terms):
        mt = mono_tex(m)
        a_ = abs(cf); sg = '-' if cf < 0 else ('+' if i else '')
        coeff = '' if (a_ == 1 and mt) else str(a_)
        out += sg + coeff + mt
    return out
def poly_sage(P):
    return '+'.join('(%s)*%s' % (cf, mono_sage(m)) for cf, m in zip(P.coefficients(), P.monomials()))
from collections import defaultdict
byq = defaultdict(lambda: S(0))
for cf, m in zip(r2.coefficients(), r2.monomials()):
    k = m.degree(Q); byq[k] += cf*m
lines_tex = []; total = ''
for k in sorted(byq, reverse=True):
    P = byq[k]
    # factor out the common monomial (q^k t^j (bd)^l ...)
    g = gcd([mm for mm in P.monomials()])
    P0 = S(P / g)
    lead_sign = 1
    if P0.lc() < 0 and len(P0.monomials()) > 1:
        lead_sign = -1; P0 = -P0
    inner = poly_tex(P0)
    pre = mono_tex(g)
    if len(P0.monomials()) > 1:
        body = pre + '\\bigl(' + inner + '\\bigr)'
    else:
        body = (str(P0.lc()) if P0.lc() != 1 else '') + pre
    lines_tex.append(('-' if lead_sign < 0 else '+') + body)
    total += ' + (%d)*(%s)*(%s)' % (lead_sign, mono_sage(g), poly_sage(P0))
lines_tex[0] = lines_tex[0].lstrip('+')
# verify
Kp.<q,t,a,b,c,d> = PolynomialRing(QQ)
val = sage_eval(total, locals=dict(q=q,t=t,a=a,b=b,c=c,d=d))
R = load('F020_results'); g0 = R[(0,2,0)]; K = g0.parent().base_ring()
C = dict(zip([tuple(e) for e in g0.exponents()], g0.coefficients()))
qq,tt,aa,bb,cc,dd = K.gens()
D = (1-qq*tt)*(1-qq*tt^4*aa*bb*cc*dd)*(1-qq^2*tt^4*aa*bb*cc*dd)
assert K(-val.subs({}) if False else -K(val)/D) == K(C[(0,0,0)]) if False else K(-K(Kp(val))/D) == K(C[(0,0,0)])
print("verified")
for L in lines_tex: print(L)
open('N020_lines.tex','w').write('\n'.join(lines_tex))
