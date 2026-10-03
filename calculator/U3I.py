# U3I: geometry table, organic shortcuts, bond dipoles, IR, functional groups
import CLIB as L
import CELEM as E
import LGEO as W
out=L.out
f=L.f
g=L.g

def vsepr():
  L.hdr("Geometry of a center")
  b=int(L.num("# bonded atoms"))
  l=int(L.num("# lone pairs",0))
  d=b+l
  gm=W.GEO.get((d,l))
  out("domains = "+str(b)+" + "+str(l)+" = "+str(d))
  if gm:
    out("e- geometry: "+gm[0])
    out("molecular geometry: "+gm[1])
    out("angle: "+gm[2])
  else:
    out("not in table")
  out("Lone pairs & double bonds take more room -> squeeze angles.")

def organic():
  L.hdr("Organic shortcut")
  for t in ("C, 4 single bonds -> tetrahedral ~109.5","C, one double -> trigonal planar ~120","C, triple or two doubles -> linear 180","N, 3 single -> pyramidal ~107","O, 2 single -> bent ~105"):
    out(" "+t)
  out("Skeletal drawings: C at every corner/line end; add H so each C has 4 bonds.")
  out("Wedge = toward you; dashed = away; line = in page.")

def dipole():
  L.hdr("Bond polarity / dipole")
  out("|delta| = dchi/(chiA+chiB); mu = |delta| e d")
  E.tipcase()
  a=E.parse(L.ask("Atom A?",""))[0][0][0]
  b=E.parse(L.ask("Atom B?",""))[0][0][0]
  ca=E.chi(a)
  cb=E.chi(b)
  ov=L.numu("chi A (blank="+g(ca)+")")
  if ov is not None:
    ca=ov
  ov=L.numu("chi B (blank="+g(cb)+")")
  if ov is not None:
    cb=ov
  d=abs(ca-cb)
  dl=d/(ca+cb)
  out("dchi = |"+g(ca)+" - "+g(cb)+"| = "+f(d,3))
  out("|delta| = "+f(d,3)+"/("+g(ca)+"+"+g(cb)+") = "+f(dl))
  neg=a if ca>cb else b
  out(neg+" gets d- ; other gets d+. Arrow points d+ -> d-.")
  out("Percent ionic ~ "+f(dl*100,3)+"%"+(" -> treat as IONIC" if dl>0.5 else ""))
  dd=L.LN("bond length (blank=skip)",1)
  if dd:
    mu=dl*L.EC*dd
    out("mu = |delta| x 1.602E-19 x d")
    out(" = "+f(dl)+" x 1.602E-19 x "+f(dd))
    out(" = "+f(mu)+" C m")
    out(" = "+f(mu/3.335e-30)+" D (1 D = 3.335E-30 C m)")

IRT=((3200,3600,"O-H, N-H stretch (broad)"),(2800,3100,"C-H stretch"),(2100,2300,"triple: C#C, C#N"),(1600,1800,"double: C=O (strong ~1700), C=C"),(1350,1500,"C-H bends (CH2, CH3)"),(1000,1300,"single: C-O, C-N, C-C"))

def irpeaks():
  L.hdr("IR bands -> bonds")
  while 1:
    w=L.numu("wavenumber cm-1 (blank=done)")
    if w is None:
      break
    hit=[t for a,b,t in IRT if a-50<=w<=b+50]
    out(" "+g(w)+": "+(" / ".join(hit) if hit else "no table match"))
  out("Stronger bonds & lighter atoms -> higher wavenumber.")

def irmol():
  m,s=W.getmol()
  if m is None:
    return
  el=m["el"]
  b=set()
  for k in range(len(m["ed"])):
    i,j=m["ed"][k]
    o=m["bo"][k]
    pr=set((el[i],el[j]))
    if "H" in pr and (pr&set(("O","N"))):
      b.add(0)
    elif pr==set(("C","H")):
      b.add(1)
    elif o==3:
      b.add(2)
    elif o==2:
      b.add(3)
    elif pr in (set(("C","O")),set(("C","N")),set(("C",))):
      b.add(5)
  for i in range(len(el)):
    if el[i]=="C" and len([1 for j,o in W.nbrs(m,i) if el[j]=="H"])>=2:
      b.add(4)
  out("Expected IR bands for "+s+":")
  for k in sorted(b):
    a,c,t=IRT[k]
    out(" "+str(a)+"-"+str(c)+": "+t)
  hetero=len(set(el))>1
  act=(len(el)>2) or hetero
  out("IR active: "+("YES" if act else "NO")+" ("+("vibrations change the dipole" if act else "homonuclear diatomic: no dipole change")+")")
  out("Greenhouse gas? "+("yes, absorbs IR" if act else "no (like N2, O2)"))
  fgroups(m)

def fgroups(m):
  el=m["el"]
  fg=[]
  for i in range(len(el)):
    nb=W.nbrs(m,i)
    if el[i]=="O":
      cs=[j for j,o in nb if el[j]=="C"]
      hs=[j for j,o in nb if el[j]=="H"]
      if len(nb)==2 and cs and hs:
        c=cs[0]
        if [1 for j,o in W.nbrs(m,c) if el[j]=="O" and o==2]:
          fg.append("carboxylic acid R-C(=O)-O-H")
        else:
          fg.append("alcohol (hydroxyl) R-O-H")
      if len(cs)==2:
        dbl=[1 for c in cs for j,o in W.nbrs(m,c) if el[j]=="O" and o==2]
        fg.append("ester R-C(=O)-O-R" if dbl else "ether R1-O-R2")
    if el[i]=="C":
      o2=[j for j,o in nb if el[j]=="O" and o==2]
      if o2:
        oth=[j for j,o in nb if j not in o2]
        if not [1 for j in oth if el[j] in ("O","N")]:
          nh=len([1 for j in oth if el[j]=="H"])
          fg.append("aldehyde R-CH=O" if nh else "ketone R1-C(=O)-R2")
    if el[i]=="N":
      if [1 for j,o in nb if el[j]=="C"]:
        fg.append("amine R-N(R1)(R2)")
  if fg:
    out("Functional groups: "+"; ".join(sorted(set(fg))))

def notes():
  L.hdr("Lewis/VSEPR notes")
  for t in ("Central atom: highest valence (most bonds); tie -> larger atom. H and F never central.",
    "Total valence e-: sum over atoms; anion add charge, cation subtract.",
    "C: 4 single / 2 single+1 double / 1 single+1 triple / 2 double. O: 2 single or 1 double + 2 lone pairs. N: 3 bonds + 1 lone pair.",
    "Resonance (O3, CO3 2-, NO3-, SO2): more than one place for a double bond; real = average.",
    "Bond length: triple < double < single; strength: triple > double > single.",
    "Odd e- -> radical (.OH, NO, NO2).",
    "Polar? Symmetric shape + identical outer atoms (CO2, CH4, CCl4, SO3, BF3, SF6) -> nonpolar. Lone pairs on central or different outer atoms (H2O, NH3, SO2, O3, CHCl3, CH2O) -> polar.",
    "ESP maps: red = negative (e- rich), blue = positive (e- poor).",
    "Functional groups: alcohol R-OH, aldehyde R-CH=O, ketone R-C(=O)-R, acid R-C(=O)-OH, ether R-O-R, amine R-NR2, aromatic = 6-C ring alt. double bonds."):
    out("- "+t)
