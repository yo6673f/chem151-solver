# U4C: classify substances for IMFs; boiling point ranking
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def classify(s):
  # returns dict for a substance: ion?, polar, hbond donor/acceptor, e-count, branching
  r={"s":s,"ion":False,"pol":False,"don":False,"acc":False,"ne":0,"br":0,"C":0}
  try:
    if s.replace("(","").startswith("NH4") or [1 for e,n in E.parse(s)[0] if E.metal(e)]:
      r["ion"]=True
      r["ne"]=E.nelec(s)
      return r
  except:
    pass
  import LSOLV
  import LGEO as W
  q=0
  t=s
  while t and t[-1] in "+-":
    q+=1 if t[-1]=="+" else -1
    t=t[:-1]
  m=LSOLV.build(t,q)
  if m is None:
    r["ne"]=E.nelec(t)
    return r
  r["pol"],why=W.polar(m)
  r["don"],r["acc"]=W.hbond(m)
  r["ne"]=sum([E.z(e) for e in m["el"]])-q
  r["C"]=m["el"].count("C")
  if q:
    r["ion"]=True
  br=0
  for i in range(len(m["el"])):
    hv=[1 for j,o in W.nbrs(m,i) if m["el"][j]!="H"]
    if m["el"][i]=="C" and len(hv)>=3:
      br+=len(hv)-2
  r["br"]=br
  return r

def lvl(r):
  if r["ion"]:
    return 3
  if r["don"]:
    return 2
  if r["pol"]:
    return 1
  return 0

LV=("dispersion only","dipole-dipole","H-bonding","ionic")

def desc(r):
  if r["ion"]:
    return r["s"]+": ionic/ion ("+str(r["ne"])+" e-)"
  t=r["s"]+": "+str(r["ne"])+" e-, "+("polar" if r["pol"] else "nonpolar")
  if r["don"]:
    t+=", H-bond (H on N/O/F)"
  if r["br"]:
    t+=", branched"
  return t

def eff(r):
  # effective "stickiness": e- count scaled by IMF type (sheet method)
  if r["s"] in ("H2O","HOH"):
    return 64.0
  k=(1.0,1.4,2.3,20.0)[lvl(r)]
  return r["ne"]*k*(1-0.03*r["br"])

def less(a,b):
  # True if a has LOWER boiling point than b; returns reason
  ea,eb=eff(a),eff(b)
  la,lb=lvl(a),lvl(b)
  if abs(ea-eb)>0.03*max(ea,eb):
    lo,hi=(a,b) if ea<eb else (b,a)
    if lvl(hi)>lvl(lo):
      w=LV[lvl(hi)]+" beats "+LV[lvl(lo)]
      if hi["ne"]<lo["ne"]:
        w+=" (even with fewer e-)"
    elif lvl(hi)<lvl(lo):
      w="much larger molecule (more e-, more surface) -> dispersion dominates"
    elif hi["ne"]==lo["ne"]:
      w="isomers: straight chain > branched (more contact)"
    else:
      w="same IMF types: more e- -> stronger dispersion"
    return ea<eb,w
  if la!=lb:
    return la<lb,LV[max(la,lb)]+" beats "+LV[min(la,lb)]+" (similar size)"
  if a["br"]!=b["br"]:
    return a["br"]>b["br"],"isomers: straight chain > branched (more contact)"
  return a["ne"]<b["ne"],"same IMF type: more e- -> stronger dispersion"

def rank():
  L.hdr("Rank boiling points")
  out("Enter formulas one per prompt (condensed for organics).")
  E.tipcase()
  rs=[]
  while 1:
    s=L.ask("Substance (blank=done)?","")
    if s=="":
      break
    rs.append(classify(s))
  out("1-2) Polarity & IMFs:")
  for r in rs:
    out(" "+desc(r))
  # insertion sort low -> high
  srt=[]
  for r in rs:
    i=0
    while i<len(srt) and less(srt[i],r)[0]:
      i+=1
    srt.insert(i,r)
  out("3-5) Lowest -> highest T_b:")
  out(" "+" < ".join([r["s"] for r in srt]))
  out("6) Justification:")
  for i in range(len(srt)-1):
    lo,hi=srt[i],srt[i+1]
    w=less(lo,hi)[1]
    out(" "+hi["s"]+" > "+lo["s"]+": "+w+", so more energy is needed to separate its molecules -> higher T_b.")
  out("(Same order for T_m, viscosity, surface tension; reverse for vapor pressure.)")
  out("Note: water is a famous exception (very strong H-bond network).")
