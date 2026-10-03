# U2B: ideal gas law (molar form), identify gas
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def getM():
  k=L.pick("Molar mass",["from formula","type value"])
  if k==1:
    E.tipcase()
    s=E.getf()
    return E.molar(s)
  return L.num("M (g/mol)")

def pvnrt():
  L.hdr("PV = nRT")
  out("R = 0.08206 L atm/(K mol). Units -> atm, L, K. Leave ONE blank.")
  P=L.P("P",1)
  V=L.V("V",1)
  T=L.T("T",1)
  k=L.pick("Amount as",["moles n","mass m (needs M)"])
  M=None
  if k==1:
    n=L.numu("n (mol)")
  else:
    m=L.M("m",1)
    M=getM()
    n=None if m is None else m/M
    if m is not None:
      out("n = m/M = "+g(m)+"/"+L.mm(M)+" = "+f(n)+" mol")
  R=L.RL
  if not L.one([P,V,T,n]):
    return
  if P is None:
    P=n*R*T/V
    out("P = nRT/V = "+f(n)+" x 0.08206 x "+g(T)+" / "+f(V))
    out(" = "+f(P)+" atm = "+f(P*760)+" mmHg")
  elif V is None:
    V=n*R*T/P
    out("V = nRT/P = "+f(n)+" x 0.08206 x "+g(T)+" / "+f(P))
    out(" = "+f(V)+" L")
  elif T is None:
    T=P*V/(n*R)
    out("T = PV/(nR) = "+f(P)+" x "+f(V)+" / ("+f(n)+" x 0.08206)")
    out(" = "+f(T)+" K = "+f(T-273.15)+" C")
  else:
    n=P*V/(R*T)
    out("n = PV/(RT) = "+f(P)+" x "+f(V)+" / (0.08206 x "+g(T)+")")
    out(" = "+f(n)+" mol")
    if M:
      out("m = n x M = "+f(n*M)+" g")
  out("N = n x N_A = "+f(n*L.NA)+" particles")

def mdn():
  L.hdr("M, d, n/V of a gas")
  k=L.pick("Find",["M = mRT/PV","d = PM/RT","n/V = P/RT","M from density"])
  P=L.P("P")
  T=L.T("T")
  R=L.RL
  if k==1:
    m=L.M("m")
    V=L.V("V")
    M=m*R*T/(P*V)
    out("M = mRT/(PV)")
    out(" = "+g(m)+" x 0.08206 x "+g(T)+" / ("+f(P)+" x "+f(V)+")")
    out(" = "+f(M)+" g/mol")
    match(M)
  elif k==2:
    M=getM()
    d=P*M/(R*T)
    out("d = PM/(RT) = "+f(P)+" x "+L.mm(M)+" / (0.08206 x "+g(T)+")")
    out(" = "+f(d)+" g/L")
  elif k==3:
    c=P/(R*T)
    out("n/V = P/(RT) = "+f(P)+"/(0.08206 x "+g(T)+")")
    out(" = "+f(c)+" mol/L")
    out(" = "+f(c*L.NA/1000)+" molecules/mL")
  else:
    d=L.num("density (g/L)")
    M=d*R*T/P
    out("M = dRT/P = "+g(d)+" x 0.08206 x "+g(T)+"/"+f(P))
    out(" = "+f(M)+" g/mol")
    match(M)

GASES=("H2","He","CH4","NH3","H2O","HF","Ne","C2H2","HCN","CO","N2","C2H4","NO","C2H6","CH2O","O2","SiH4","PH3","H2S","HCl","F2","Ar","CO2","N2O","C3H8","C3H6","NO2","O3","C4H10","SO2","BF3","Cl2","SO3","HBr","CF4","Kr","SF6","Xe")

def match(M):
  out("3) Match M to a formula:")
  r=[]
  for s in GASES:
    r.append((abs(E.molar(s,0)-M),s))
  r.sort()
  for d,s in r[:3]:
    out(" "+s+" "+L.mm(E.molar(s,0))+" (off "+f(d,2)+")")
  if L.yes("Search from elements"):
    els=L.ask("Elements (e.g. C,H,O)?","C,H,O").replace(" ","").split(",")
    ms=[E.mass(e) for e in els]
    best=[]
    def rec(i,M0,cnt):
      if i==len(els):
        if M0>0 and abs(M0-M)<=max(0.5,0.01*M):
          best.append((abs(M0-M),cnt[:]))
        return
      for c in range(0,13):
        if M0+c*ms[i]>M*1.02+0.5:
          break
        cnt.append(c)
        rec(i+1,M0+c*ms[i],cnt)
        cnt.pop()
    rec(0,0.0,[])
    best.sort()
    for d,c in best[:6]:
      s="".join([els[i]+(str(c[i]) if c[i]>1 else "") for i in range(len(els)) if c[i]])
      out(" "+s+" = "+L.mm(E.molar(s,0)))
    if not best:
      out(" none found")

def stpv():
  L.hdr("Molar volume / gas mixes")
  out("1 mol ideal gas: STP (0 C, 1 atm) = 22.41 L; 25 C, 1 atm = 24.47 L")
  k=L.pick("Find",["V from mol (STP)","mol from V (STP)","moles of part in mix"])
  if k==1:
    n=L.num("n (mol)")
    out("V = n x 22.41 = "+f(n*22.41)+" L")
  elif k==2:
    V=L.V("V")
    out("n = V/22.41 = "+f(V/22.41)+" mol")
  else:
    out("Volume % = mole %.")
    p=L.num("vol % of gas (e.g. 78.08)")
    n=L.num("total moles")
    out("n gas = "+g(p)+"/100 x "+g(n)+" = "+f(p/100*n)+" mol")
    out("mole fraction = "+f(p/100))

def run():
  L.run_menu("2 Gases",
    ["PV = nRT solver","M, d, n/V; ID gas","Molar volume / mixes","Pollutant conversions","Ratio problems"],
    [pvnrt,mdn,stpv,"U2P.poll","U2P.ratios"])
