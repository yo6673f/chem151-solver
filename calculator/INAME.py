# INAME: name <-> formula for ionic and binary molecular compounds
from IONIC import *
import IONIC as I

PRE=("","mono","di","tri","tetra","penta","hexa","hepta","octa","nona","deca")

def name(s,show=1):
  # name any binary/ionic compound from formula
  com={"H2O":"water","NH3":"ammonia","CH4":"methane","H2O2":"hydrogen peroxide"}
  if s in com:
    return com[s]
  if isionic(s):
    r=split(s)
    if r is None:
      return None
    c,cq,cn,a,aq,an=r
    if show:
      out("Ionic (metal + nonmetal):")
      out(" anion "+ion(a,aq)+" ("+ionname(a,aq)+")")
      out(" total (-) = "+str(an)+" x ("+str(aq)+") = "+str(an*aq))
      out(" so "+str(cn)+" "+c+" carry +"+str(-an*aq)+" -> "+ion(c,cq))
      if E.z(c) and not fixed(c):
        out(" metal has several charges -> Roman numeral ("+ROM[cq]+")")
    return ionname(c,cq)+" "+ionname(a,aq)
  lst,q=E.parse(s)
  if len(lst)!=2:
    return None
  if show:
    out("Molecular (nonmetals): prefixes; no mono on 1st; 2nd gets -ide")
  (e1,n1),(e2,n2)=lst
  a=("" if n1==1 else PRE[n1])+E.name(e1).lower()
  r=IDE.get(e2,E.name(e2).lower()+"ide")
  p=PRE[n2]
  if p[-1:] in "ao" and r[0] in "ao":
    p=p[:-1]
  return a+" "+p+r

def frm(nm):
  # name -> formula
  w=nm.lower().replace("("," (").split()
  if not w:
    return None
  com={"water":"H2O","ammonia":"NH3","methane":"CH4"}
  if nm.lower() in com:
    return com[nm.lower()]
  rom=0
  if len(w)>1 and w[1].startswith("("):
    rs=w[1].strip("()").upper()
    rom=ROM.index(rs) if rs in ROM else 0
    w=[w[0]]+w[2:]
  first=w[0]
  an=" ".join(w[1:])
  cat=None
  if first=="ammonium":
    cat="NH4"
    cq=1
  else:
    for i in range(1,E.NEL+1):
      r=E.row(i)
      if r[1].lower()==first and E.metal(r[0]):
        cat=r[0]
        cq=rom or fixed(cat)
        if cat=="Hg" and cq==1:
          cat="Hg2"
          cq=2
  if cat:
    a=None
    for p in POLY:
      if p[2]==an and p[1]<0:
        a=p[0]
        aq=p[1]
        break
    if a is None:
      for k in IDE:
        if IDE[k]==an:
          a=k
          aq=anq(k)
    if a is None or not cq:
      return None
    fs,nc,na=unit(cat,cq,a,aq)
    out("Ions: "+ion(cat,cq)+" and "+ion(a,aq))
    out("Criss-cross |charges|, reduce:")
    out(" "+str(nc)+" x ("+str(cq)+") + "+str(na)+" x ("+str(aq)+") = 0")
    return fs
  # molecular
  res=""
  for part in w:
    n=1
    for k in range(10,0,-1):
      pk=PRE[k]
      rs=[]
      if part.startswith(pk):
        rs.append(part[len(pk):])
      if pk[-1] in "ao" and part.startswith(pk[:-1]):
        rs.append(part[len(pk)-1:])
      hit=0
      for rest in rs:
        if _el(rest):
          n=k
          part=rest
          hit=1
          break
      if hit:
        break
    sym=_el(part)
    if not sym:
      return None
    res+=sym+(str(n) if n>1 else "")
  return res

def _el(word):
  for i in range(1,E.NEL+1):
    r=E.row(i)
    if r[1].lower()==word:
      return r[0]
  for k in IDE:
    if IDE[k]==word:
      return k
  return None

# ---- solubility (sheet tables) ----
