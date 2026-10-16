# n=1: check the paper's nonsymmetric Askey-Wilson recurrence (prop_AW_recurrence)
P = PolynomialRing(QQ,'q,t,a,b,c,d,z'); q,t,a,b,c,d,z = P.gens(); F=P.fraction_field()
def s1(f): return f.subs({z:1/z})
def s0(f): return f.subs({z:q/z})
def T1(f): f=F(f); return -b*d*f+(1-b*z)*(1-d*z)/(1-z**2)*(s1(f)-f)
def T0(f): f=F(f); return -a*c/q*f+(1-a/z)*(1-c/z)/(1-q/z**2)*(s0(f)-f)
t0=-a*c/q; tn=-b*d
def Tinv(T,tp,f): return (T(f)-(tp-1)*F(f))/tp
def Y(f): return T1(T0(F(f)))   # Y_1 = T_n T_0 for n=1
def coefz(f,e):
    f=F(f); num=P(f.numerator()); den=P(f.denominator())
    # den independent of z? make Laurent: multiply
    R=PolynomialRing(F.base_ring() if False else QQ['q','t','a','b','c','d'].fraction_field(),'Z'); 
    raise
K=QQ['q','t','a','b','c','d'].fraction_field()
L=LaurentPolynomialRing(K,'zz'); zz=L.gen()
def lau(f):
    f=F(f); num=P(f.numerator()); den=P(f.denominator())
    dz=set(e[6] for e in den.exponents()); assert len(dz)==1; m=dz.pop()
    dp=sum(cf*K(prod(g**e[i] for i,g in enumerate([q,t,a,b,c,d]))) for e,cf in den.dict().items())
    return sum(K(cf)*K(prod(g**e[i] for i,g in enumerate(K.gens())))/K(dp)*zz**(e[6]-m) for e,cf in num.dict().items())
def co(f,e):
    g=lau(f); return K(g[e])
def act(i,eta):
    return -1-eta if i==0 else -eta
def step(E,eta,i):
    new=act(i,eta)
    w = T1(E) if i==1 else Tinv(T0,t0,F(E)*z)
    al=co(w,new); assert al!=0
    yE=co(Y(E),eta); Yw=Y(w); yw=co(Yw,new)/al
    beta=K((lau(Yw)-yw*lau(w))[eta])/(yE-yw)
    return new,(w-beta*F(E))/al
# E_{(-k)} by path 0 -s0-> -1 -s1-> 1 -s0-> -2 ...
E=F(1); eta=0; Es={}
word=[0,1,0,1,0]
for i in word:
    eta,E=step(E,eta,i)
    if eta<0: Es[eta]=E
Fm={k:Es[-k] for k in [1,2,3]}
Fp={k:Tinv(T1,tn,Fm[k]) for k in [1,2,3]}
for k in [1,2,3]:
    assert co(Fp[k],k)==1 and co(Fm[k],-k)==1
def vk(k): return q**(k-1)*(a*c+q-a*c*q**k-a*b*c*d*q**k)/(1-a*b*c*d*q**(2*k-1))
def wk(k): return q**k*(-(a+c)+a*c*(b+d)*q**k)/(1-a*b*c*d*q**(2*k))
def vt(k): return (1+b*d-b*d*q**k-a*b*c*d*q**(k-1))/(1-a*b*c*d*q**(2*k-1))
def wt(k): return (-(b+d)+b*d*(a+c)*q**k)/(1-a*b*c*d*q**(2*k))
print("F1 base:", lau(Fp[1]-(z+wk(0)))==0, " F-1 base:", lau(Fm[1]-(1/z+wt(0)))==0)
for k in [1,2]:
    u=wk(k)*vt(k); ut=vk(k)*wt(k)
    r1=lau(Fp[k+1]-((z+u)*Fp[k]+(vk(k)*z+wk(k))*Fm[k]))
    r2=lau(Fm[k+1]-((vt(k)/z+wt(k))*Fp[k]+(1/z+ut)*Fm[k]))
    print("k=",k, r1==0, r2==0)
