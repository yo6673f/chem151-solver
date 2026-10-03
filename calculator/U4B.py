# U4B: ions, formula units, unit cells, Coulomb
import CLIB as L
import CELEM as E
import IONIC as I
out=L.out
f=L.f
g=L.g

def ionof():
  L.hdr("Ion of an element")
  E.tipcase()
  e=E.parse(L.ask("Element?",""))[0][0][0]
  Z=E.z(e)
  p,gp=E.pg(Z)
  out(e+": period "+str(p)+", group "+str(gp))
  out("Usual charges: grp1 +1, 2 +2, 13 +3, 15 -3, 16 -2, 17 -1")
  q=I.fixed(e) if E.metal(e) else I.anq(e)
  if gp==13:
    q=3
  if q==0:
    out("Transition metal: several possible charges -> Roman numeral (Fe3+ = iron(III), Cu+ = copper(I))")
    return
  ne=Z-q
  out("Ion: "+I.ion(e,q)+" ("+I.ionname(e,q)+")")
  out(("Metal loses " if q>0 else "Nonmetal gains ")+str(abs(q))+" e- -> "+str(ne)+" e-")
  for z,nm in ((2,"He"),(10,"Ne"),(18,"Ar"),(36,"Kr"),(54,"Xe"),(86,"Rn")):
    if ne==z:
      out("same config as noble gas "+nm)

def getion(p):
  out(p+": symbol/polyatomic (Na, Fe, NH4, SO4)")
  s=""
  while s=="":
    s=L.ask(" ion?","").replace(" ","")
  q=None
  for pp in I.POLY:
    if pp[0]==s:
      q=pp[1]
  if q is None:
    e=E.parse(s)[0][0][0]
    s=e
    q=I.fixed(e) if E.metal(e) else I.anq(e)
    if E.pg(E.z(e))[1]==13 and E.metal(e):
      q=3
    if q==0 or not E.metal(e) and q==0:
      q=int(L.num(" charge of "+e))
    else:
      v=L.numu(" charge (blank="+str(q)+")")
      if v is not None:
        q=int(v)
  return s,q

def funit():
  L.hdr("Formula unit from ions")
  c,cq=getion("Cation")
  a,aq=getion("Anion")
  out("1) Ions: "+I.ion(c,cq)+" and "+I.ion(a,aq))
  fs,nc,na=I.unit(c,cq,a,aq)
  out("2) Criss-cross charge sizes: "+c+" x"+str(abs(aq))+", "+a+" x"+str(abs(cq))+"; reduce -> "+str(nc)+":"+str(na))
  out(" check: "+str(nc)+"("+str(cq)+") + "+str(na)+"("+str(aq)+") = "+str(nc*cq+na*aq))
  out("3) Parentheses if >1 polyatomic ion")
  out("Formula: "+fs)
  out("Name: "+I.ionname(c,cq)+" "+I.ionname(a,aq))
  M=E.molar(fs,0)
  out("M = "+L.mm(M)+" g/mol")

def cell():
  L.hdr("Unit cell counting")
  out("corner 1/8, edge 1/4, face 1/2, inside 1")
  res=[]
  for nm in ("Cation","Anion"):
    s=L.ask(nm+" symbol?","")
    a=L.num(" # at corners",0)
    b=L.num(" # on edges",0)
    c=L.num(" # on faces",0)
    d=L.num(" # inside",0)
    t=a/8.0+b/4.0+c/2.0+d
    out(" "+s+": "+g(a)+"/8 + "+g(b)+"/4 + "+g(c)+"/2 + "+g(d)+" = "+g(t))
    res.append((s,t))
  x,y=res[0][1],res[1][1]
  k=1
  for k in range(1,13):
    if abs(x*k/min(x,y)-round(x*k/min(x,y)))<0.01 and abs(y*k/min(x,y)-round(y*k/min(x,y)))<0.01:
      break
  nx=int(round(x*k/min(x,y)))
  ny=int(round(y*k/min(x,y)))
  d=L.gcd(nx,ny)
  out("Ratio "+g(x)+":"+g(y)+" = "+str(nx//d)+":"+str(ny//d))
  out("Formula: "+res[0][0]+(str(nx//d) if nx//d>1 else "")+res[1][0]+(str(ny//d) if ny//d>1 else ""))
  out("NaCl: 6 neighbors (FCC); CsCl: 8 (BCC)")

def coul():
  L.hdr("F = K q1 q2 / r^2")
  z1=L.num("charge # ion 1 (e.g. 2)")
  z2=L.num("charge # ion 2 (e.g. -1)")
  k=L.pick("r from",["distance","ionic radii (pm)"])
  if k==1:
    r=L.LN("r between centers")
  else:
    a=L.num(" cation radius pm")
    b=L.num(" anion radius pm")
    r=(a+b)*1e-12
    out(" r = "+g(a)+"+"+g(b)+" = "+g(a+b)+" pm = "+f(r)+" m")
  e=L.EC
  F=L.KC*z1*z2*(e/r)*(e/r)
  out("q = charge # x 1.602E-19 C")
  out("F = 8.988E9 x ("+g(z1)+"x1.602E-19)("+g(z2)+"x1.602E-19)/("+f(r)+")^2")
  out(" = "+f(F)+" N")
  out(("negative -> ATTRACTION" if F<0 else "positive -> repulsion"))
  out("Bigger charges, smaller ions -> stronger -> higher T_m, T_b, dH fus.")

def notes():
  L.hdr("Ionic notes")
  for t in ("Cations SMALLER than their atom; anions LARGER; same charge bigger down a group & right-to-left.",
    "Brittle: a blow shifts layers so like charges line up and repel -> shatters.",
    "Conduct only molten or dissolved (ions free to move).",
    "Dissolves if ion-water beats ion-ion (lower E_p) AND mixing gives more configurations; one of two -> T-dependent; neither -> insoluble.",
    "Polyatomic: OH- NO3- CO3 2- SO4 2- PO4 3- NH4+ CN- CH3COO-"):
    out("- "+t)

def run():
  L.run_menu("4 Ionic compounds",
    ["Ion of an element","Formula unit from ions","Unit cell counting","Coulomb's law F","Rank ionic m.p.","Soluble? / dissolve eq","Precipitation reaction","Notes"],
    [ionof,funit,cell,coul,"U4S.mprank","U4S.solq","U4S.precip",notes])
