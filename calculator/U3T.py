# U3T: periodic trends, ion sizes, valence, radicals, notes
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def why(a,b,prop):
  pa,ga=E.pg(E.z(a))
  pb,gb=E.pg(E.z(b))
  if pa==pb:
    return "same period: "+(a if ga>gb else b)+" more protons pull same shell in"
  if ga==gb:
    return "same group: "+(a if pa>pb else b)+" has a new shell farther out"
  return "different period & group: compare data"

def trends():
  L.hdr("Rank periodic trends")
  k=L.pick("Property",["atomic radius","1st ionization E","electronegativity"])
  E.tipcase()
  s=L.ask("Elements (comma sep)?","")
  els=[E.parse(x.strip())[0][0][0] for x in s.split(",") if x.strip()]
  data=[E.rad,E.ie,lambda Z:E.chi(E.sym(Z))][k-1]
  un=[" pm"," kJ/mol",""][k-1]
  r=[(data(E.z(e)),e) for e in els]
  r.sort()
  r.reverse()
  nm=["radius","IE1","chi"][k-1]
  out("Across a period: radius dec, IE inc, chi inc. Down a group: radius inc, IE dec, chi dec.")
  out("Largest -> smallest "+nm+":")
  for v,e in r:
    p,gr=E.pg(E.z(e))
    out(" "+e+" (per "+str(p)+", grp "+str(gr)+"): "+g(v)+un)
  out("Reasons:")
  for i in range(len(r)-1):
    out(" "+r[i][1]+" > "+r[i+1][1]+": "+why(r[i][1],r[i+1][1],k))
  out("Why: e- held tighter when closer & more nuclear charge. Big jump in radius / drop in IE when a new shell starts (He->Li, Ne->Na).")
  if k==2:
    out("(Small exceptions: Be>B, N>O in IE)")

RION={"Li+":76,"Na+":102,"K+":138,"Rb+":152,"Cs+":167,"Be2+":45,"Mg2+":72,"Ca2+":100,"Sr2+":118,"Ba2+":135,"Al3+":54,"F-":133,"Cl-":181,"Br-":196,"I-":220,"O2-":140,"S2-":184,"Se2-":198,"N3-":146,"P3-":212,"Ga3+":62,"Zn2+":74,"Ag+":115}

def ions():
  L.hdr("Rank sizes of atoms/ions")
  out("Enter like Na+, O2-, Ne, Mg2+")
  s=L.ask("Species (comma sep)?","")
  r=[]
  for x in s.split(","):
    x=x.strip()
    if not x:
      continue
    y=x.replace("^","")
    sg=0
    while y and y[-1] in "+-":
      sg+=1 if y[-1]=="+" else -1
      y=y[:-1]
    d=""
    while y and y[-1].isdigit():
      d=y[-1]+d
      y=y[:-1]
    q=sg*(int(d) if d else 1) if sg else 0
    e=E.parse(y)[0][0][0]
    if q==0:
      v=E.rad(E.z(e))
    else:
      v=RION.get(e+(str(abs(q)) if abs(q)>1 else "")+("+" if q>0 else "-"),0)
    Z=E.z(e)
    r.append((v,x,Z,Z-q))
  if [x for x in r if x[2]!=x[3]]:
    r=[(x[0] if x[2]!=x[3] else 0,x[1],x[2],x[3]) for x in r]
  r.sort()
  r.reverse()
  same=len(set([x[3] for x in r]))==1
  if same:
    r.sort(key=lambda x:x[2])
  out("Largest -> smallest:")
  for v,x,Z,ne in r:
    out(" "+x+": "+(g(v)+" pm" if v else "n/a")+" (Z="+str(Z)+", "+str(ne)+" e-)")
  if same:
    out("ISOELECTRONIC: same # e-; more protons pull harder -> smaller")
  out("Cations smaller than their atom (lost e-/shell); anions larger (more e- repulsion). Same charge: bigger down a group.")

def valence():
  L.hdr("Valence (bonding capacity)")
  E.tipcase()
  s=L.ask("Element?","")
  e=E.parse(s)[0][0][0]
  ve=E.nval(e)
  mx=2 if e in ("H","He") else 8
  out("valence e- = "+str(ve)+" (group 1-2: group #; 13-18: group # - 10)")
  b=mx-ve
  out("bonds = max occupancy - valence e-")
  out(" = "+str(mx)+" - "+str(ve)+" = "+str(b))
  p,gr=E.pg(E.z(e))
  if p>=3 and gr>=15:
    out("Period "+str(p)+": can exceed octet (P 5, S 6: SF6, H2SO4)")
  out("Typical: H1 C4 N3 O2 F,Cl,Br,I 1, noble gases 0")

def radical():
  L.hdr("Radical check")
  E.tipcase()
  s=E.getf()
  lst,q=E.parse(s)
  t=0
  st=""
  for e,n in lst:
    t+=E.nval(e)*n
    st+=e+":"+str(n)+"x"+str(E.nval(e))+" "
  t-=q
  out("valence e-: "+st+("charge "+str(q) if q else ""))
  out(" total = "+str(t))
  out("ODD -> RADICAL (unpaired e-, very reactive)" if t%2 else "even -> not a radical")

def vr():
  k=L.pick("Find",["bonds of element","radical? (odd e-)","light/atom notes"])
  if k==1:
    valence()
  elif k==2:
    radical()
  else:
    notes()

def notes():
  L.hdr("Light & atom notes")
  for t in ("Longer lambda <-> lower nu <-> lower E. Regions long->short: radio/micro (rotation), IR (vibration), visible 700-400 nm (e- excitation), UV 400-10 nm (bond breaking, ionization), X-ray/gamma <10 nm (core e-).",
    "Energy transferred only in whole photons; atoms/molecules have discrete energy levels.",
    "Absorption & emission lines of an element at the SAME lambda; each line = one transition.",
    "Photoelectric: e- ejected only if h nu > threshold; KE grows with nu, not intensity.",
    "Ground state = lowest levels filled; excited = e- promoted.",
    "e- as wave: smaller space -> shorter lambda -> higher KE; delocalization lowers KE.",
    "Blackbody: hotter -> peak shifts to shorter lambda (red->white->blue), more energy.",
    "Molecular bands broad: neighbors shift levels slightly.",
    "Per mole: E = N_A h nu; kJ/mol ~ 1.197E5/lambda(nm) ~ 0.01197 x wavenumber.",
    "Structure ID: % comp -> empirical; M+ -> molecular formula; IR -> bonds; draw Lewis; check MS fragments (M-15 CH3, M-17 OH)."):
    out("- "+t)
