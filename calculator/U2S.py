# U2S: molarity, mass %, dilution
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def getn(p):
  k=L.pick(p+" given as",["grams","moles"])
  if k==2:
    return L.num("moles")
  E.tipcase()
  s=E.getf("Solute formula")
  M=E.molar(s)
  m=L.M("mass")
  n=m/M
  out("n = m/M = "+g(m)+"/"+L.mm(M)+" = "+f(n)+" mol")
  return n

def molar():
  L.hdr("[A] = n / V")
  k=L.pick("Find",["[A] (molarity)","moles/grams","volume"])
  if k==1:
    out("1) grams -> moles:")
    n=getn("Solute")
    out("2) volume -> liters:")
    V=L.V("Solution volume")
    c=n/V
    out("3) [A] = n/V = "+f(n)+"/"+f(V))
    out("  = "+f(c)+" M (mol/L)")
  elif k==2:
    c=L.num("[A] (M)")
    V=L.V("Volume")
    n=c*V
    out("n = [A] x V = "+g(c)+" x "+f(V)+" = "+f(n)+" mol")
    if L.yes("Grams too"):
      E.tipcase()
      s=E.getf("Solute formula")
      M=E.molar(s)
      out("m = n x M = "+f(n)+" x "+L.mm(M)+" = "+f(n*M)+" g")
  else:
    c=L.num("[A] (M)")
    n=getn("Solute")
    V=n/c
    out("V = n/[A] = "+f(n)+"/"+g(c)+" = "+f(V)+" L = "+f(V*1000)+" mL")

def mpct():
  L.hdr("Mass % solutions")
  k=L.pick("Find",["mass %","mass % -> molarity","molarity -> mass %"])
  if k==1:
    a=L.num("mass solute (g)")
    b=L.num("mass solution (g) (blank: solute+solvent)",-1)
    if b<0:
      sv=L.num("mass solvent (g)")
      b=a+sv
      out("m solution = "+g(a)+" + "+g(sv)+" = "+g(b)+" g")
    p=a/b*100
    out("mass % = (m solute/m solution) x 100")
    out(" = "+g(a)+"/"+g(b)+" x 100 = "+f(p)+" %")
    return
  E.tipcase()
  s=E.getf("Solute formula")
  M=E.molar(s)
  d=L.num("Solution density g/mL",1.0)
  if k==2:
    p=L.num("mass %")
    out("Assume 1 L = 1000 mL solution:")
    ms=1000*d
    out(" mass soln = 1000 x "+g(d)+" = "+f(ms)+" g")
    a=ms*p/100
    out(" solute = "+f(ms)+" x "+g(p)+"/100 = "+f(a)+" g")
    n=a/M
    out(" moles = "+f(a)+"/"+L.mm(M)+" = "+f(n)+" mol")
    out(" [A] = "+f(n)+" mol / 1 L = "+f(n)+" M")
  else:
    c=L.num("[A] (M)")
    out("Assume 1 L solution:")
    a=c*M
    out(" solute = "+g(c)+" mol x "+L.mm(M)+" = "+f(a)+" g")
    ms=1000*d
    out(" mass soln = 1000 mL x "+g(d)+" = "+f(ms)+" g")
    out(" mass % = "+f(a)+"/"+f(ms)+" x 100 = "+f(a/ms*100)+" %")

def dil():
  L.hdr("Dilution M1V1 = M2V2")
  out("Blank = unknown. Same V units.")
  a=L.numu("M1")
  b=L.numu("V1")
  c=L.numu("M2")
  d=L.numu("V2")
  if not L.one([a,b,c,d]):
    return
  out("moles solute stay the same: [A]1 V1 = [A]2 V2")
  if a is None:
    a=c*d/b
    out("M1 = M2 V2/V1 = "+f(a))
  elif b is None:
    b=c*d/a
    out("V1 = M2 V2/M1 = "+f(b))
  elif c is None:
    c=a*b/d
    out("M2 = M1 V1/V2 = "+f(c))
  else:
    d=a*b/c
    out("V2 = M1 V1/M2 = "+f(d))
    out("solvent to add = "+f(d-b))
