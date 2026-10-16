# Check Macdonald VI (6.24)(ii) psi_{lam/mu} = branching coefficient in P_lam(x_1..x_n) = sum psi x_n^k P_mu(x_1..x_{n-1})
R = QQ['q','t'].fraction_field(); q,t = R.gens()
Sym = SymmetricFunctions(R); P = Sym.macdonald().P()
qq = QQ(1)/7; tt = QQ(2)/5; RF = RealField(200)
def fnum(a,b,N=300):
    v = RF(1)
    for k in range(N):
        v *= (1 - qq**(a+k)*tt**(b+1))/(1 - qq**(a+k+1)*tt**b)
    return v
def psi(lam, mu):
    l = len(mu); L = list(lam)+[0]*6; M = list(mu)+[0]*6
    v = RF(1)
    for i in range(1,l+1):
        for j in range(i,l+1):
            v *= fnum(M[i-1]-M[j-1], j-i)*fnum(L[i-1]-L[j], j-i)/(fnum(L[i-1]-M[j-1], j-i)*fnum(M[i-1]-L[j], j-i))
    return v
for lam,n in [([2,1],3),([3,1],3),([2,2],3),([3,2],3),([3,2,1],3),([4,2,1],4)]:
    lam = Partition(lam)
    X = PolynomialRing(R, ['x%d'%i for i in range(1,n+1)]); xs = X.gens()
    Pl = X(P(lam).expand(n, alphabet=xs))
    worst = 0; count=0
    for k in range(lam.size()+1):
        # coefficient of x_n^k as polynomial in x_1..x_{n-1}
        res = sum(c*X.monomial(*e) for e,c in Pl.dict().items() if e[n-1]==k)
        res = res/xs[n-1]**k if k>0 else res
        res = X(res)
        mus = [mu for mu in Partitions(lam.size()-k) if len(mu)<=n-1]
        mus.sort(key=lambda m: list(m), reverse=True)  # lex decreasing refines dominance
        for mu in mus:
            mono = X.monomial(*(list(mu)+[0]*(n-len(mu))))
            c = res.monomial_coefficient(mono)
            if mu.size()>0:
                Pm = X(P(mu).expand(n-1, alphabet=xs[:n-1]))
            else: Pm = X(1)
            res = res - c*Pm
            strip = lam.contains(mu) and all(mu.get_part(i) >= lam.get_part(i+1) for i in range(len(lam)))
            if strip:
                err = abs(RF(c.subs(q=qq,t=tt)) - psi(lam,mu)); worst=max(worst,err); count+=1
            else:
                assert c==0, ("non-strip coefficient", lam, mu)
        assert res==0
    print(lam, n, "strips checked:", count, "max err:", worst)
