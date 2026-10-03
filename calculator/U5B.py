# U5B: stoichiometry, limiting reactant, yield, CO2 emission factor, AFR
import CLIB as L
import CELEM as E
import BAL
out=L.out
f=L.f
g=L.g

def geteq():
  out("Equation like C3H8+O2>CO2+H2O")
  E.tipcase()
  s=L.ask("Equation?","")
  sp=BAL.sides(s)
  co=BAL.solve(s)
  if co is None:
    out("Could not balance")
    return None
  allsp=sp[0]+sp[1]
  st=[(str(co[j])+" " if co[j]>1 else "")+allsp[j] for j in range(len(allsp))]
  out("1) Balanced: "+" + ".join(st[:len(sp[0])])+" -> "+" + ".join(st[len(sp[0]):]))
  Ms=[E.molar(x,0) for x in allsp]
  return allsp,co,Ms,len(sp[0])

def choose(t,allsp):
  return L.pick(t,allsp)-1

def moles(nm,M):
  k=L.pick(nm+" given as",["grams","moles","liquid mL+density","gas P,V,T","solution M x V"])
  if k==1:
    m=L.M(" mass")
    n=m/M
    out(" n = m/M = "+g(m)+"/"+L.mm(M)+" = "+f(n)+" mol")
  elif k==2:
    n=L.num(" moles")
  elif k==3:
    v=L.V(" volume")*1000
    d=L.num(" density g/mL")
    m=v*d
    n=m/M
    out(" m = rho x V = "+g(d)+" x "+g(v)+" = "+f(m)+" g")
    out(" n = m/M = "+f(n)+" mol")
  elif k==4:
    P=L.P(" P")
    V=L.V(" V")
    T=L.T(" T")
    n=P*V/(L.RL*T)
    out(" n = PV/RT = "+f(P)+"x"+f(V)+"/(0.08206x"+g(T)+")")
    out("  = "+f(n)+" mol")
  else:
    c=L.num(" [A] (M)")
    V=L.V(" V")
    n=c*V
    out(" n = [A] x V = "+g(c)+" x "+f(V)+" = "+f(n)+" mol")
  return n

def report(nm,n,M):
  out(" mass "+nm+" = n x M = "+f(n)+" x "+L.mm(M)+" = "+f(n*M)+" g")
  out(" particles = n x N_A = "+f(n*L.NA))
  out(" gas at STP: V = n x 22.41 = "+f(n*22.41)+" L")

def atob():
  L.hdr("Grams A -> grams B")
  r=geteq()
  if r is None:
    return
  sp,co,Ms,nr=r
  a=choose("Given A:",sp)
  b=choose("Find B:",sp)
  out("2) moles of "+sp[a]+":")
  nA=moles(sp[a],Ms[a])
  nB=nA*co[b]/co[a]
  out("3) mol B = mol A x (coef B/coef A)")
  out(" = "+f(nA)+" x "+str(co[b])+"/"+str(co[a])+" = "+f(nB)+" mol "+sp[b])
  out("4) grams B = mol B x M_B:")
  report(sp[b],nB,Ms[b])

def limit():
  L.hdr("Limiting reactant & yield")
  r=geteq()
  if r is None:
    return
  sp,co,Ms,nr=r
  out("1) Convert every reactant to moles:")
  ns=[]
  for i in range(nr):
    if sp[i]=="O2" and L.yes("O2 in excess (skip)"):
      ns.append(None)
      continue
    ns.append(moles(sp[i],Ms[i]))
  out("2-4) mol / coefficient (smallest = limiting):")
  best=None
  for i in range(nr):
    if ns[i] is None:
      continue
    v=ns[i]/co[i]
    out(" "+sp[i]+": "+f(ns[i])+"/"+str(co[i])+" = "+f(v))
    if best is None or v<best[0]:
      best=(v,i)
  lim=best[1]
  out("LIMITING: "+sp[lim])
  idx=[i for i in range(nr) if ns[i] is not None]
  if len(idx)==2:
    x,y=idx
    ra=ns[x]/ns[y]
    rs=co[x]/float(co[y])
    out(" (ratio check: actual "+sp[x]+"/"+sp[y]+" = "+f(ra)+" vs stoich "+f(rs)+(" -> "+sp[x]+" limiting" if ra<rs else (" -> "+sp[y]+" limiting" if ra>rs else " -> both used up"))+")")
  out("5) Theoretical yield (from limiting only):")
  for j in range(nr,len(sp)):
    n=ns[lim]*co[j]/co[lim]
    out(" "+sp[j]+": "+f(ns[lim])+" x "+str(co[j])+"/"+str(co[lim])+" = "+f(n)+" mol = "+f(n*Ms[j])+" g")
  out("6) Excess left over:")
  for i in idx:
    if i==lim:
      continue
    used=ns[lim]*co[i]/co[lim]
    left=ns[i]-used
    out(" "+sp[i]+": "+f(ns[i])+" - "+f(used)+" = "+f(left)+" mol = "+f(left*Ms[i])+" g")
  act=L.numu("Actual yield g (blank=skip)")
  if act is not None:
    j=nr+choose("of product",sp[nr:])
    th=ns[lim]*co[j]/co[lim]*Ms[j]
    out("% yield = actual/theoretical x 100")
    out(" = "+g(act)+"/"+f(th)+" x 100 = "+f(act/th*100)+" %")

def ef():
  L.hdr("CO2 emission factor & AFR")
  E.tipcase()
  s=E.getf("Fuel formula")
  c=E.counts(s)
  M=E.molar(s)
  pc=c.get("C",0)*12.01/M*100
  out("%C = "+str(c.get("C",0))+"x12.01/"+L.mm(M)+" x100 = "+f(pc)+" %")
  EF=pc/100*44.01/12.01
  out("EF = %C/100 x 44.01/12.01")
  out(" = "+f(EF)+" kg CO2 per kg fuel")
  x=c.get("C",0)
  y=c.get("H",0)
  z=c.get("O",0)
  nO2=x+y/4.0-z/2.0
  out("Combustion: 1 "+s+" + "+g(nO2)+" O2 -> "+str(x)+" CO2 + "+g(y/2.0)+" H2O")
  mO2=nO2*32.00
  out("m_O2 per mol fuel = "+g(nO2)+" x 32.00 = "+f(mO2)+" g (per "+L.mm(M)+" g fuel)")
  afr=mO2/0.233/M
  out("AFR = (m_O2/0.233)/m_fuel")
  out(" = ("+f(mO2)+"/0.233)/"+L.mm(M)+" = "+f(afr))
  out("(air = 23.3% O2 by mass; 21% by moles)")
  a=L.numu("Actual AFR (blank=skip)")
  if a is not None:
    if a<afr:
      out("Below stoich -> RICH: fuel excess, O2 limiting -> CO + unburned fuel")
    elif a>afr:
      out("Above stoich -> LEAN: O2 in excess")
    else:
      out("Stoichiometric")
  m=L.numu("kg fuel burned (blank=skip)")
  if m is not None:
    out("CO2 = "+g(m)+" x "+f(EF)+" = "+f(m*EF)+" kg")
  out("Higher %C (bigger molecules, more double/triple bonds) -> higher EF. Octane EF ~3.08; gasoline AFR ~14.7")

def run():
  L.run_menu("5 Stoichiometry",
    ["Grams A -> grams B","Limiting reactant/yield","CO2 EF and AFR"],
    [atob,limit,ef])
