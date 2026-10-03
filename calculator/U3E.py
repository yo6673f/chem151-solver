# U3E: electron configurations, PES, identify element
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

ORD=("1s","2s","2p","3s","3p","4s","3d","4p","5s","4d","5p","6s","4f","5d","6p","7s","5f","6d","7p")

CAP={"s":2,"p":6,"d":10,"f":14}

EXC={24:"[Ar] 4s1 3d5",29:"[Ar] 4s1 3d10",41:"[Kr] 5s1 4d4",42:"[Kr] 5s1 4d5",44:"[Kr] 5s1 4d7",45:"[Kr] 5s1 4d8",46:"[Kr] 4d10",47:"[Kr] 5s1 4d10",78:"[Xe] 6s1 4f14 5d9",79:"[Xe] 6s1 4f14 5d10"}

NOB=((2,"He"),(10,"Ne"),(18,"Ar"),(36,"Kr"),(54,"Xe"),(86,"Rn"))

def fill(ne):
  c=[]
  for s in ORD:
    if ne<=0:
      break
    k=min(ne,CAP[s[1]])
    c.append([s,k])
    ne-=k
  return c

def cstr(c):
  return " ".join([s+str(k) for s,k in c])

def core(c):
  # noble gas notation
  tot=0
  best=(0,"")
  run=0
  for i in range(len(c)):
    run+=c[i][1]
    for z,s in NOB:
      if run==z and c[i][1]==CAP[c[i][0][1]] and i<len(c)-1:
        best=(i+1,s)
  if best[0]==0:
    return cstr(c)
  return "["+best[1]+"] "+cstr(c[best[0]:])

def remove(c,k):
  # remove k e- from highest n first (then highest l)
  key=lambda x:(int(x[0][0]),"spdf".index(x[0][1]))
  for _ in range(k):
    c.sort(key=key)
    c[-1][1]-=1
    if c[-1][1]==0:
      c.pop()
  c.sort(key=lambda x:ORD.index(x[0]))
  return c

def unp(s,k):
  o=CAP[s[1]]//2
  return k if k<=o else 2*o-k

def box(s,k):
  o=CAP[s[1]]//2
  b=[0]*o
  for i in range(k):
    b[i%o]+=1
  return "".join(["["+("ud" if x==2 else ("u " if x==1 else "  "))+"]" for x in b])

def config():
  L.hdr("Electron configuration")
  E.tipcase()
  s=L.ask("Element or Z?","")
  try:
    Z=int(s)
  except:
    Z=E.z(E.parse(s)[0][0][0])
  q=int(L.num("Charge",0))
  sym=E.sym(Z)
  ne=Z-q
  out("1) Z = "+str(Z)+" ("+E.nameZ(Z)+"); #e = Z - charge = "+str(ne))
  out("2) Fill order: 1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p..")
  out(" s holds 2, p 6, d 10, f 14")
  if q>0:
    c=fill(Z)
    out(" neutral: "+cstr(c))
    c=remove(c,q)
    out(" cation: remove "+str(q)+" e- from highest n first (4s before 3d)")
  else:
    c=fill(ne)
  out("Config: "+cstr(c))
  out("Noble gas: "+core(c))
  out("3) Check: "+"+".join([str(k) for s1,k in c])+" = "+str(sum([k for s1,k in c])))
  if q==0 and Z in EXC:
    out("NOTE: real config is exception: "+EXC[Z]+" (half/full d more stable)")
  for z,nm in NOB:
    if ne==z:
      out("Isoelectronic with "+nm)
  nmax=max([int(s1[0]) for s1,k in c])
  ve=sum([k for s1,k in c if int(s1[0])==nmax and s1[1] in "sp"])
  out("Valence e- (outermost n="+str(nmax)+" s+p) = "+str(ve))
  u=sum([unp(s1,k) for s1,k in c])
  out("Unpaired e- = "+str(u)+(" (paramagnetic)" if u else " (diamagnetic)"))
  last=c[-1]
  out("Last subshell "+last[0]+": "+box(last[0],last[1])+" (singly before pairing - Hund)")
  if L.yes("Show PES"):
    pes(c)

def pes(c):
  out("PES: one peak per subshell; height = # e-; 1s at highest IE:")
  for s1,k in c:
    out(" "+s1+" |"+"="*k+" "+str(k))
  out("(listed high IE -> low IE)")

def ident():
  L.hdr("Identify from config/PES")
  k=L.pick("Given",["config text","PES peak heights"])
  if k==1:
    s=L.ask("Config (e.g. [Ne]3s2 3p1)?","")
    t=0
    for z,nm in NOB:
      if "["+nm.upper()+"]" in s.upper():
        t+=z
    s2=s.split("]")[-1]
    for tok in s2.replace(","," ").split():
      tok=tok.strip()
      if len(tok)>=3 and tok[0].isdigit() and tok[1].lower() in "spdf":
        t+=int(tok[2:])
  else:
    out("Enter heights (# e-) of each peak; blank=done")
    t=0
    hs=[]
    while 1:
      v=L.numu(" peak "+str(len(hs)+1))
      if v is None:
        break
      hs.append(int(v))
      t+=int(v)
    c=fill(t)
    out(" Subshells: "+cstr(c))
  out("Total e- = "+str(t))
  if 0<t<=E.NEL:
    out("Neutral atom -> Z = "+str(t)+": "+E.sym(t)+" ("+E.nameZ(t)+")")

# ---------- trends ----------
