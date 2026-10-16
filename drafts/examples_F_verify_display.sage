R = load('F_results')
K = R[(0,1,0)].parent().base_ring()
q,t,a,b,c,d = K.gens()
M = (a+c)*(b+d) - t*(a*c*(b+d)^2 + b*d*(a+c)^2 - a*b*c*d) + t^2*a*b*c*d*((a+c)*(b+d)-1) - t^3*(1-t)*(a*b*c*d)^2
X = q*t^2*a*c*(b+d) - a - c
disp = {
 (0,1,0): {(0,1,0): 1, (0,0,0): t*(t^2*a*c*(b+d)-a-c)/(1-t^4*a*b*c*d)},
 (0,1,-1): {(0,1,-1): 1,
            (0,1,0): (t*b*d*(a+c)-b-d)/(1-t^2*a*b*c*d),
            (0,0,-1): t*(t*a*c*(b+d)-a-c)/(1-t^2*a*b*c*d),
            (0,0,0): t*M/((1-t^2*a*b*c*d)*(1-t^3*a*b*c*d))},
 (0,2,-1): {(0,2,-1): 1,
            (1,1,-1): (1-t)/(1-q*t^2),
            (-1,1,-1): q*(1-t)/(1-q*t^2),
            (0,2,0): (t*b*d*(a+c)-b-d)/(1-t^2*a*b*c*d),
            (1,0,-1): (1-t)*q*t^2*X/((1-q*t^2)*(1-q^2*t^4*a*b*c*d)),
            (-1,0,-1): (1-t)*q^2*t^2*X/((1-q*t^2)*(1-q^2*t^4*a*b*c*d)),
            (0,-1,-1): t*q*X/(1-q^2*t^4*a*b*c*d)},
}
for eta, D in disp.items():
    g = R[eta]
    C = dict(zip([tuple(e) for e in g.exponents()], g.coefficients()))
    if eta != (0,2,-1): assert set(C) == set(D), ("monomials differ", eta)
    for e, v in D.items():
        assert K(C[e]) == K(v), (eta, e)
    print(eta, "OK", len(C), "monomials")
