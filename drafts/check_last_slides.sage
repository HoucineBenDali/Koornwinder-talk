# (a) quadratic + braid relations of the Demazure-Lusztig operators (n=3), (b) K_(2,0) slide transcription (n=2)
P = PolynomialRing(QQ,'q,t,a,b,c,d,z1,z2,z3'); q,t,a,b,c,d,z1,z2,z3 = P.gens(); Z=[z1,z2,z3]; F=P.fraction_field(); n=3
def s(i,f):
    if 1<=i<n: return f.subs({Z[i-1]:Z[i],Z[i]:Z[i-1]})
    if i==n: return f.subs({Z[n-1]:1/Z[n-1]})
    return f.subs({Z[0]:q/Z[0]})
def T(i,f):
    f=F(f)
    if 1<=i<n: x,y=Z[i-1],Z[i]; return t*f+(t*x-y)/(x-y)*(s(i,f)-f)
    if i==n: x=Z[n-1]; return -b*d*f+(1-b*x)*(1-d*x)/(1-x**2)*(s(i,f)-f)
    x=Z[0]; return -a*c/q*f+(1-a/x)*(1-c/x)/(1-q/x**2)*(s(0,f)-f)
tp={0:-a*c/q,1:t,2:t,3:-b*d}
f = F(z1**2*z2**-1 + 3*z3 + z1*z3**-2 + 7)
for i in range(4):
    g=T(i,f); assert T(i,g)-(tp[i]-1)*g-tp[i]*f==0, ("quadratic",i)
def W(word,f):
    for i in reversed(word): f=T(i,f)
    return f
assert W([0,1,0,1],f)==W([1,0,1,0],f); print("T0T1T0T1 ok")
assert W([2,3,2,3],f)==W([3,2,3,2],f); print("T2T3T2T3 ok")
assert W([1,2,1],f)==W([2,1,2],f); print("T1T2T1 ok")
assert W([0,2],f)==W([2,0],f) and W([0,3],f)==W([3,0],f) and W([1,3],f)==W([3,1],f); print("commutations ok")
print("quadratic ok")
# (b)
R = load('/dev/null') if False else None
