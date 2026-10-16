R = load('F021_results'); g = R[(0,2,1)]
K = g.parent().base_ring(); Kp = K.ring()
q,t,a,b,c,d = Kp.gens()
BIG = 12  # factors with more terms than this go outside the fraction (breakable)

def is_mono_minus_one(h):
    ms = h.monomials(); cs = h.coefficients()
    if len(ms) == 2 and 1 in ms:
        i = ms.index(1); j = 1 - i
        return cs[j] == 1 and cs[i] == -1, ms[j]
    return False, None

def split(poly):
    """returns sign, list of (latex, sage, exponent, big?)"""
    if poly == 1: return 1, []
    f = factor(poly); sign = f.unit(); out = []
    for h, e in f:
        ok, m = is_mono_minus_one(h)
        if ok:
            sign *= (-1)**e
            out.append(("(1-%s)" % latex(m), "(1-(%s))" % m, e, False))
        elif len(h.monomials()) == 1:
            out.append((latex(h), "(%s)" % h, e, False))
        else:
            big = len(h.monomials()) > BIG
            out.append(("\\bigl(%s\\bigr)" % latex(h) if big else "(%s)" % latex(h), "(%s)" % h, e, big))
    return sign, out

def prodlatex(fs):
    s = ""
    for L, S, e, big in fs:
        s += L + ("^{%d}" % e if e > 1 else "")
    return s

def prodsage(fs):
    return "*".join("%s^%d" % (S, e) for L, S, e, big in fs) if fs else "1"

def zmono(ex):
    s = ""
    for k, e in enumerate(ex):
        if e == 1: s += "z_{%d}" % (k+1)
        elif e != 0: s += "z_{%d}^{%d}" % (k+1, e)
    return s

terms = sorted(zip(g.exponents(), g.coefficients()), key=lambda x: (-sum(x[0]), x[0]))
latex_terms = []; ok = True
for ex, cf in terms:
    sn, fn = split(Kp(cf.numerator())); sd, fd = split(Kp(cf.denominator()))
    sign = sn * sd
    small = [x for x in fn if not x[3]]; big = [x for x in fn if x[3]]
    num = prodlatex(small) or "1"
    den = prodlatex(fd)
    coef = ("\\frac{%s}{%s}" % (num, den)) if fd else (num if num != "1" else "")
    coef += prodlatex(big)
    body = coef + zmono(ex)
    if body == "": body = "1"
    latex_terms.append((sign, body, sum(ex) == 3))
    # verification: rebuild in Sage from the same pieces
    val = sign * eval(preparse(prodsage(fn))) / eval(preparse(prodsage(fd)))
    if K(val) != K(cf):
        ok = False; print("MISMATCH", ex)
print("verified:", ok, len(latex_terms), "terms")
out = []
for i, (sign, body, top) in enumerate(latex_terms):
    s = ("" if i == 0 and sign > 0 else ("+" if sign > 0 else "-")) + body
    out.append("\\textcolor{uibred}{%s}" % s if top else s)
open('F021_terms.tex', 'w').write("\n".join(out))
print(sum(len(x) for x in out), "characters")
