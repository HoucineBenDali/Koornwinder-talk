R = load('E_results')
K = R[(0,1,0)].parent().base_ring()
q,t,a,b,c,d = K.gens()
N = (q+a*c) - q*t*a*c - q*t^2*(q+a*b*c*d) - q*t^3*a*c*(1-q*b*d) + q^2*t^4*a*c*(1-b*d) + q^2*t^5*a*b*c*d
disp = {
 (0,1,0): {
   (0,1,0): 1,
   (0,0,1): q*t*(1-t)/(1-q*t^2),
   (0,0,-1): q*t*(1-t)/(1-q*t^2),
   (-1,0,0): -a*c*(1-t)*(1-q*t^3)*(1-q*t^2*b*d)/((1-q*t^2)*(1-q*t^4*a*b*c*d)),
   (0,-1,0): t*N/((1-q*t^2)*(1-q*t^4*a*b*c*d)),
   (0,0,0): -t*(1-q*t^3)*(a+c+q*(b+d)-q*t^2*(a*b*c+a*b*d+a*c*d+b*c*d))/((1-q*t^2)*(1-q*t^4*a*b*c*d)),
 },
 (0,1,-1): {
   (0,1,-1): 1,
   (-1,1,0): a*b*c*d*(1-t)/(1-t*a*b*c*d),
   (0,1,0): (t*b*d*(a+c)-b-d)/(1-t*a*b*c*d),
 },
 (0,2,-1): {
   (0,2,-1): 1,
   (1,1,-1): (1-t)/(1-q*t),
   (0,1,-2): q*(1-t)/(1-q*t),
 },
}
for eta, D in disp.items():
    g = R[eta]
    C = dict(zip([tuple(e) for e in g.exponents()], g.coefficients()))
    full = (eta == (0,1,0))
    if full: assert set(C) == set(D), "monomial sets differ"
    for e, v in D.items():
        assert K(C[e]) == K(v), (eta, e)
    print(eta, "OK", len(C), "monomials")
c00 = dict(zip([tuple(e) for e in R[(0,2,-1)].exponents()], R[(0,2,-1)].coefficients()))[(0,0,0)]
print("const term of E_(0,2,-1): numerator terms", len(c00.numerator().monomials()), "denominator", factor(c00.denominator()))
