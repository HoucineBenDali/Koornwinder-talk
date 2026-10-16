R = load('F020_results'); g = R[(0,2,0)]
K = g.parent().base_ring(); Kp = K.ring()
C = dict(zip([tuple(e) for e in g.exponents()], g.coefficients()))
q,t,a,b,c,d = Kp.gens()
D = (1-q*t)*(1-q*t^4*a*b*c*d)*(1-q^2*t^4*a*b*c*d)
N = Kp(-(C[(0,0,0)]*K(D)))
# express N in s=a+c, p=ac, u=b+d, v=bd
S.<Q,T,s,p,u,v> = PolynomialRing(QQ, order='deglex')
B.<q2,t2,a2,b2,c2,d2,s2,p2,u2,v2> = PolynomialRing(QQ, order='lex(4),degrevlex(6)')
I = B.ideal([s2-a2-c2, p2-a2*c2, u2-b2-d2, v2-b2*d2])
G = I.groebner_basis()
Nb = N.subs({q:q2,t:t2,a:a2,b:b2,c:c2,d:d2})
r = Nb.reduce(G)
print("vars in reduced:", r.variables())
r2 = S(str(r).replace('q2','Q').replace('t2','T').replace('s2','s').replace('p2','p').replace('u2','u').replace('v2','v'))
print(len(r2.monomials()), "terms in (q,t,a+c,ac,b+d,bd)")
# group by powers of q
from collections import defaultdict
byq = defaultdict(lambda: S(0))
for cf, m in zip(r2.coefficients(), r2.monomials()):
    e = m.degree(Q); byq[e] += cf*m/Q**e
for e in sorted(byq, reverse=True):
    print("q^%d:" % e, factor(byq[e]))
save(r2, 'N020_sym')
