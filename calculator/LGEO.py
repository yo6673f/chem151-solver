# LGEO: VSEPR geometry, polarity, H-bonding; getmol() asks for a formula
import CLIB as L
import CELEM as E
out=L.out

GEO={(2,0):("linear","linear","180"),(3,0):("trigonal planar","trigonal planar","120"),
(3,1):("trigonal planar","bent","<120 (~117)"),(4,0):("tetrahedral","tetrahedral","109.5"),
(4,1):("tetrahedral","trigonal pyramidal","~107"),(4,2):("tetrahedral","bent","~104.5"),
(5,0):("trigonal bipyramidal","trigonal bipyramidal","90, 120"),(5,1):("trigonal bipyramidal","seesaw","<90, <120"),
(5,2):("trigonal bipyramidal","T-shaped","<90"),(5,3):("trigonal bipyramidal","linear","180"),
(6,0):("octahedral","octahedral","90"),(6,1):("octahedral","square pyramidal","<90"),(6,2):("octahedral","square planar","90")}

SYM=((2,0),(3,0),(4,0),(5,0),(6,0),(5,3),(6,2))

def labels(el):
  c={}
  t={}
  for e in el:
    t[e]=t.get(e,0)+1
  r=[]
  for e in el:
    c[e]=c.get(e,0)+1
    r.append(e+(str(c[e]) if t[e]>1 else ""))
  return r

def nbrs(m,i):
  r=[]
  for k in m["nb"][i]:
    a,b=m["ed"][k]
    r.append((b if a==i else a,m["bo"][k]))
  return r

def dom(m,i):
  lp=m["lp"][i]
  return len(m["nb"][i])+lp//2+lp%2,lp//2+lp%2

def bpol(a,b):
  if (a,b) in (("C","H"),("H","C")):
    return 0.0
  return E.chi(b)-E.chi(a)

def centers(m):
  return [i for i in range(len(m["el"])) if len(m["nb"][i])>=2]

def polar(m):
  # returns (polar bool, reasons list)
  el=m["el"]
  why=[]
  if len(el)==2:
    d=abs(bpol(el[0],el[1]))
    return d>0.001,["diatomic, dchi = "+str(round(d,2))]
  pol=False
  lab=labels(el)
  for i in centers(m):
    d,l=dom(m,i)
    nb=nbrs(m,i)
    ds=[round(bpol(el[i],el[j]),2) for j,o in nb]
    sym=(d,l) in SYM
    if not sym and l>0:
      why.append(lab[i]+": lone pairs, asymmetric -> dipoles don't cancel")
      pol=True
    elif max([abs(x) for x in ds])<0.001:
      why.append(lab[i]+": bonds nonpolar (C-H ~ nonpolar)")
    elif sym and len(set(ds))==1 and len(set([el[j] for j,o in nb]))==1:
      why.append(lab[i]+": symmetric, identical outer atoms -> cancel")
    elif sym and len(set(ds))==1:
      why.append(lab[i]+": symmetric, equal dipoles -> cancel")
    else:
      why.append(lab[i]+": different outer atoms -> dipoles don't cancel")
      pol=True
  return pol,why

def hbond(m):
  el=m["el"]
  don=False
  acc=False
  for k in range(len(m["ed"])):
    a,b=m["ed"][k]
    if (el[a]=="H" and el[b] in ("N","O","F")) or (el[b]=="H" and el[a] in ("N","O","F")):
      don=True
  for i in range(len(el)):
    if el[i] in ("N","O","F") and m["lp"][i]>=2:
      acc=True
  return don,acc

def getmol():
  out("Formula, or condensed formula for organics (CH3CH2OH, CH3COOH). # = triple bond.")
  E.tipcase()
  s=L.ask("Formula?","")
  if s=="":
    return None,None
  q=0
  t=s.replace(" ","")
  if "^" in t:
    t,qs=t.split("^")
    q=E._qp(qs)
  else:
    while t and t[-1] in "+-":
      q+=1 if t[-1]=="+" else -1
      t=t[:-1]
  if q==0:
    q=int(L.num("Charge",0))
  import LSOLV
  m=LSOLV.build(t,q)
  if m is None:
    out("Could not build a structure (rings or unusual bonding?).")
    return None,t
  return m,t
