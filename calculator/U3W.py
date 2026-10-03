# U3W: Lewis structure steps, VSEPR, polarity, IMFs
import CLIB as L
import CELEM as E
import LGEO as W
out=L.out
f=L.f
g=L.g

BS=("-","=","#")

def show(m,s):
  el=m["el"]
  lab=W.labels(el)
  n=len(el)
  q=m["q"]
  L.hdr("Lewis: "+s+(" ("+str(q)+")" if q else ""))
  heavy=[i for i in range(n) if el[i]!="H"]
  if m["tag"].startswith("central"):
    c=m["tag"].split()[1]
    out("1) Central atom: "+c+" (appears once/least electroneg.; H and F never central)")
  elif m["tag"]=="chain":
    out("1) Skeleton: chain of heavy atoms in formula order (isomers possible - type a condensed formula to choose)")
  elif m["tag"]=="condensed":
    out("1) Skeleton from condensed formula")
  else:
    out("1) Skeleton: "+m["tag"])
  cnt={}
  od=[]
  for e in el:
    if e not in cnt:
      od.append(e)
    cnt[e]=cnt.get(e,0)+1
  st=" + ".join([e+" "+str(cnt[e])+"x"+str(E.nval(e)) for e in od])
  out("2) Valence e-: "+st+(" - ("+str(q)+") charge" if q else "")+" = "+str(m["tot"]))
  nb=len(m["ed"])
  out("3) "+str(nb)+" single bonds x 2 = "+str(2*nb)+" e-; left "+str(m["tot"]-2*nb))
  out("4) Lone pairs on outer atoms (to 8; H gets 2), leftovers on central.")
  mb=[]
  for k in range(nb):
    if m["bo"][k]>1:
      a,b=m["ed"][k]
      mb.append(lab[a]+BS[m["bo"][k]-1]+lab[b])
  if mb:
    out("5) Short of 8 -> move lone pairs into bonds: "+", ".join(mb))
  else:
    out("5) No multiple bonds needed.")
  out("Structure:")
  for i in heavy:
    nbs=W.nbrs(m,i)
    hs=[1 for j,o in nbs if el[j]=="H"]
    parts=[BS[o-1]+lab[j] for j,o in nbs if el[j]!="H"]
    if hs:
      parts.append("-H"+("x"+str(len(hs)) if len(hs)>1 else ""))
    lp=m["lp"][i]
    tx=" "+lab[i]+": "+" ".join(parts)
    tx+="; lone pairs "+str(lp//2)+(" + 1 unpaired e-" if lp%2 else "")
    ecount=lp+2*sum([o for j,o in nbs])
    tx+=" ("+str(ecount)+" e-)"
    out(tx)
  used=2*sum(m["bo"])+sum(m["lp"])
  out("6) Check: "+str(2*sum(m["bo"]))+" bonding + "+str(sum(m["lp"]))+" lone = "+str(used)+" e- "+("OK" if used==m["tot"] else "!!"))
  fc=m["fc"]
  nz=[lab[i]+"("+("+" if fc[i]>0 else "")+str(fc[i])+")" for i in range(n) if fc[i]]
  if nz:
    out(" Formal charges: "+" ".join(nz))
  if m["nres"]>1:
    out("RESONANCE: "+str(m["nres"])+" forms (move the double bond; draw all with <->). Real molecule = average; identical bonds of intermediate length.")
  if m["rad"]>=0:
    out("Odd e- count -> RADICAL (unpaired e- on "+lab[m["rad"]]+")")
  for i in heavy:
    e2=m["lp"][i]+2*sum([o for j,o in W.nbrs(m,i)])
    if e2>8:
      out(lab[i]+" has "+str(e2)+" e-: expanded octet (period 3+)")
    elif e2<8 and m["rad"]!=i and el[i]!="H":
      out(lab[i]+" has only "+str(e2)+" e- (B/Be/Al allowed)")
  if q:
    out("Ion: draw in brackets [ ] with charge "+str(q)+" outside.")

def geom(m):
  el=m["el"]
  lab=W.labels(el)
  cs=W.centers(m)
  if not cs:
    out("Geometry: diatomic -> linear")
    return
  out("VSEPR (domain = each lone pair or bonded atom; single/double/triple = 1):")
  for i in cs:
    d,l=W.dom(m,i)
    gm=W.GEO.get((d,l),("?","?","?"))
    out(" "+lab[i]+": "+str(d-l)+" bonded + "+str(l)+" lone = "+str(d)+" domains")
    out("  e- geom "+gm[0]+"; shape "+gm[1].upper()+"; angle "+gm[2])
    nb=W.nbrs(m,i)
    if l and d<5:
      out("  lone pairs squeeze angles smaller")
    if [1 for j,o in nb if o>1] and l==0 and d>2:
      out("  double bond squeezes nearby angles a bit")

def polarity(m):
  el=m["el"]
  seen=[]
  out("Bond polarity (Pauling chi):")
  for a,b in m["ed"]:
    x,y=el[a],el[b]
    if (x,y) in seen or (y,x) in seen:
      continue
    seen.append((x,y))
    ca,cb=E.chi(x),E.chi(y)
    d=abs(ca-cb)
    if d<0.001:
      out(" "+x+"-"+y+": dchi 0 -> nonpolar")
      continue
    neg=x if ca>cb else y
    pos=y if neg==x else x
    dl=d/(ca+cb)
    tag=" (C-H treated nonpolar)" if set((x,y))==set(("C","H")) else ""
    out(" "+x+"-"+y+": dchi "+f(d,2)+", delta "+f(dl,2)+"; "+pos+" d+ -> "+neg+" d-"+tag)
  p,why=W.polar(m)
  for w in why:
    out(" "+w)
  out("Molecule is "+("POLAR" if p else "NONPOLAR (mu = 0)"))
  return p

def imf(m,p):
  don,acc=W.hbond(m)
  ne=sum([E.z(e) for e in m["el"]])-m["q"]
  out("IMFs: dispersion (all; "+str(ne)+" e-)"+(", dipole-dipole" if p else "")+(", H-bonding (H on N/O/F)" if don else ""))
  if acc and not don:
    out(" can ACCEPT H-bonds (N/O/F lone pair) from water etc.")

def lewis():
  m,s=W.getmol()
  if m is None:
    return
  show(m,s)
  L.pause()
  geom(m)
  p=polarity(m)
  imf(m,p)
