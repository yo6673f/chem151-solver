# U2D: empirical / molecular formulas, combustion analysis, mass %
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def ratio(els,ms):
  # els: symbols, ms: grams -> returns integer subscripts with steps
  out("2) grams / molar mass -> moles:")
  mol=[]
  for i in range(len(els)):
    M=E.mass(els[i])
    n=ms[i]/M
    mol.append(n)
    out(" "+els[i]+": "+f(ms[i])+"/"+f(M)+" = "+f(n))
  sm=min(mol)
  out("3) divide by smallest ("+f(sm)+"):")
  r=[x/sm for x in mol]
  for i in range(len(els)):
    out(" "+els[i]+": "+f(r[i]))
  mult=1
  for k in range(1,9):
    ok=1
    for x in r:
      if abs(x*k-round(x*k))>0.1:
        ok=0
    if ok:
      mult=k
      break
  if mult>1:
    out("4) not whole -> x"+str(mult)+" (.5->x2, .33/.67->x3, .25/.75->x4, .2->x5)")
    for i in range(len(els)):
      out(" "+els[i]+": "+f(r[i])+" x "+str(mult)+" = "+f(r[i]*mult))
  else:
    out("4) all ~whole numbers")
  sub=[int(round(x*mult)) for x in r]
  fs="".join([els[i]+(str(sub[i]) if sub[i]>1 else "") for i in range(len(els))])
  out("5) Empirical formula: "+fs)
  return fs,sub

def molform(els,sub):
  Me=0.0
  for i in range(len(els)):
    Me+=sub[i]*E.mass(els[i])
  out("M(empirical) = "+L.mm(Me))
  M=L.numu("M of molecule (M+ peak; blank=skip)")
  if M is None:
    return
  k=int(round(M/Me))
  out("6) k = M/M(emp) = "+g(M)+"/"+L.mm(Me)+" = "+f(M/Me)+" ~ "+str(k))
  fs="".join([els[i]+(str(sub[i]*k) if sub[i]*k>1 else "") for i in range(len(els))])
  out("Molecular formula: "+fs)

def frompct():
  L.hdr("Empirical from mass %")
  E.tipcase()
  els=[]
  ms=[]
  k=L.pick("Data as",["mass %","grams"])
  out("1) Assume 100.0 g -> % = grams" if k==1 else "1) grams given")
  while 1:
    e=L.ask("Element (blank=done)?","")
    if e=="":
      break
    e=E.parse(e)[0][0][0]
    v=L.numu(" "+("%" if k==1 else "grams")+" (blank=rest)")
    if v is None:
      v=(100.0 if k==1 else L.num(" total sample g"))-sum(ms)
      out(" "+e+" = by difference = "+f(v))
    els.append(e)
    ms.append(v)
  fs,sub=ratio(els,ms)
  molform(els,sub)

def comb():
  L.hdr("Combustion analysis")
  out("C_xH_yO_z + O2 -> CO2 + H2O")
  ms=L.num("sample mass (g)")
  mc=L.num("mass CO2 (g)")
  mh=L.num("mass H2O (g)")
  hasO=L.yes("Sample may contain O")
  C=mc*12.01/44.01
  H=mh*2*1.008/18.02
  out("m_C = m_CO2 x 12.01/44.01")
  out(" = "+g(mc)+" x 12.01/44.01 = "+f(C)+" g")
  out("m_H = m_H2O x 2(1.008)/18.02")
  out(" = "+g(mh)+" x 2.016/18.02 = "+f(H)+" g")
  els=["C","H"]
  m=[C,H]
  if hasO:
    O=ms-C-H
    out("m_O = m_sample - m_C - m_H")
    out(" = "+g(ms)+" - "+f(C)+" - "+f(H)+" = "+f(O)+" g")
    if O>0.001*ms:
      els.append("O")
      m.append(O)
    else:
      out(" ~0 -> no O")
  out("Mass %: "+", ".join([els[i]+" "+f(m[i]/ms*100)+"%" for i in range(len(els))]))
  fs,sub=ratio(els,m)
  molform(els,sub)

def pct():
  L.hdr("Mass % from formula")
  E.tipcase()
  s=E.getf()
  M=E.molar(s)
  out("%X = (subscript x M_X)/M x 100")
  for e,n in E.parse(s)[0]:
    out(" "+e+": "+str(n)+"x"+g(E.mass(e))+"/"+L.mm(M)+"x100 = "+f(n*E.mass(e)/M*100)+"%")

def run():
  L.run_menu("2 Formulas",
    ["Empirical from % or g","Combustion (CO2, H2O)","Mass % from formula"],
    [frompct,comb,pct])
