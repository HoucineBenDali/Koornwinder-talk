# Check psi_{lam/mu} = prod_{s in R-C} b_mu(s)/b_lam(s), b_lam(s)=(1-q^a t^(l+1))/(1-q^(a+1) t^l), against the product formula
load('psi_check.sage') if False else None
R = QQ['q','t'].fraction_field(); q,t = R.gens()
Sym = SymmetricFunctions(R); P = Sym.macdonald().P()
def b(lam, i, j):  # cell (i,j) 0-based row i, col j
    a = lam[i]-j-1; l = lam.conjugate()[j]-i-1
    return (1-q**a*t**(l+1))/(1-q**(a+1)*t**l)
def psi_hook(lam, mu):
    rows = {i for i in range(len(lam)) if lam[i] > mu.get_part(i)}
    cols = {j for j in range(lam[0]) if lam.conjugate().get_part(j) > mu.conjugate().get_part(j)}
    v = R(1)
    for i in rows:
        for j in range(mu.get_part(i)):
            if j not in cols: v *= b(mu,i,j)/b(lam,i,j)
    return v
for lam,n in [([2,1],3),([3,1],3),([3,2],3),([3,2,1],3),([4,2,1],4)]:
    lam = Partition(lam)
    X = PolynomialRing(R, ['x%d'%i for i in range(1,n+1)]); xs = X.gens()
    Pl = X(P(lam).expand(n, alphabet=xs)); bad=0; cnt=0
    for k in range(lam.size()+1):
        res = X(sum(c*X.monomial(*e) for e,c in Pl.dict().items() if e[n-1]==k)/xs[n-1]**k)
        mus = sorted([mu for mu in Partitions(lam.size()-k) if len(mu)<=n-1], key=list, reverse=True)
        for mu in mus:
            c = res.monomial_coefficient(X.monomial(*(list(mu)+[0]*(n-len(mu)))))
            res -= c*(X(P(mu).expand(n-1, alphabet=xs[:n-1])) if mu.size()>0 else X(1))
            if c!=0:
                cnt+=1
                if c != psi_hook(lam,mu): bad+=1
    print(lam, "strips:", cnt, "mismatches:", bad)
