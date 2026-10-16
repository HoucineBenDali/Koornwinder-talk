R = load('F020_results'); g = R[(0,2,0)]
K = g.parent().base_ring(); Kp = K.ring(); q,t,a,b,c,d = K.gens()
C = dict(zip([tuple(e) for e in g.exponents()], g.coefficients()))
X = q*t^2*a*c*(b+d) - a - c
A = (1-t)/(1-q*t); B = q*t*(1-t)*X/((1-q*t)*(1-q^2*t^4*a*b*c*d))
chk = {(0,2,0):1,(1,1,0):A,(0,1,1):q*A,(-1,1,0):q*A,(0,1,-1):q*A,(1,0,0):B,(0,0,1):B,(0,0,-1):B,(-1,0,0):q*B,(0,-1,0):q*t*X/(1-q^2*t^4*a*b*c*d)}
for e,v in chk.items(): assert K(C[e])==K(v), e
print("checked", len(chk))
print("z2 coeff denominator:", factor(C[(0,1,0)].denominator()))
print("const denominator:", factor(C[(0,0,0)].denominator()))
D1 = (1-q*t)*(1-q^2*t^4*a*b*c*d)
N1 = Kp(-(C[(0,1,0)]*K(D1))); print("N1 (coef = -N1/D1):", N1)
