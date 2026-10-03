# IONIC: ions, formula units, solubility rules
import CLIB as L
import CELEM as E
out=L.out

POLY=(("NH4",1,"ammonium"),("H3O",1,"hydronium"),("Hg2",2,"mercury(I)"),
("OH",-1,"hydroxide"),("NO3",-1,"nitrate"),("NO2",-1,"nitrite"),
("HCO3",-1,"hydrogen carbonate"),("CO3",-2,"carbonate"),
("HSO4",-1,"hydrogen sulfate"),("SO4",-2,"sulfate"),("SO3",-2,"sulfite"),
("H2PO4",-1,"dihydrogen phosphate"),("HPO4",-2,"hydrogen phosphate"),("PO4",-3,"phosphate"),
("CN",-1,"cyanide"),("CH3COO",-1,"acetate"),("C2H3O2",-1,"acetate"),
("ClO4",-1,"perchlorate"),("ClO3",-1,"chlorate"),("ClO2",-1,"chlorite"),("ClO",-1,"hypochlorite"),
("MnO4",-1,"permanganate"),("Cr2O7",-2,"dichromate"),("CrO4",-2,"chromate"),
("SCN",-1,"thiocyanate"),("C2O4",-2,"oxalate"),("O2",-2,"peroxide"))

IDE={"H":"hydride","N":"nitride","O":"oxide","F":"fluoride","P":"phosphide","S":"sulfide","Cl":"chloride","Se":"selenide","Br":"bromide","I":"iodide","C":"carbide","As":"arsenide","Te":"telluride","Si":"silicide","B":"boride"}

FIX={"Ag":1,"Zn":2,"Cd":2,"Al":3,"Ga":3,"Sc":3}

ROM=("","I","II","III","IV","V","VI","VII","VIII")

def fixed(sym):
  # usual cation charge or 0 if variable
  if sym in FIX:
    return FIX[sym]
  p,g=E.pg(E.z(sym))
  if g in (1,2) and sym!="H":
    return g
  return 0

def anq(sym):
  p,g=E.pg(E.z(sym))
  if sym=="H":
    return -1
  if g>=15 and g<=17:
    return g-18
  return 0

def ionname(ion,q):
  for p in POLY:
    if p[0]==ion and (p[1]>0)==(q>0):
      return p[2]
  if q>0:
    nm=E.name(ion).lower()
    if not fixed(ion):
      nm+="("+ROM[q]+")"
    return nm
  return IDE.get(ion,E.name(ion).lower()+"ide")

def _tok(s,i):
  # read a unit at s[i]: (group) or element; returns (text,end)
  if s[i]=="(":
    j=s.index(")",i)
    return s[i+1:j],j+1
  cs=False
  for ch in s:
    if ch.islower():
      cs=True
  e=E._sym(s,i,cs)
  if e is None:
    raise ValueError("bad formula")
  return e,i+len(e)

def _cnt(s,i):
  j=i
  while j<len(s) and s[j].isdigit():
    j+=1
  return (int(s[i:j]) if j>i else 1),j

def split(s):
  # ionic formula -> (cation,cq,ccount,anion,aq,acount) or None
  s=s.replace(" ","").split(".")[0]
  while s and s[-1] in "+-":
    s=s[:-1]
  cat=None
  i=0
  for p in ("NH4","Hg2","H3O"):
    if s.startswith("("+p+")"):
      cat=p
      i=len(p)+2
      break
    if s.startswith(p):
      cat=p
      i=len(p)
      break
  if cat is None:
    cat,i=_tok(s,0)
  cn,i=_cnt(s,i)
  rest=s[i:]
  if rest=="":
    return None
  an=None
  if rest[0]=="(":
    an,j=_tok(rest,0)
    acnt,j=_cnt(rest,j)
  else:
    for p in POLY:
      if p[1]<0 and p[0]!="O2" and rest.startswith(p[0]):
        k,j=_cnt(rest,len(p[0]))
        if j==len(rest):
          an=p[0]
          acnt=k
          break
    if an is None:
      an,j=_tok(rest,0)
      acnt,j=_cnt(rest,j)
  if j!=len(rest):
    return None
  aq=0
  if E.z(an):
    aq=anq(an)
  else:
    for p in POLY:
      if p[0]==an and p[1]<0:
        aq=p[1]
  if aq==0:
    return None
  cq=-aq*acnt
  if cq%cn:
    return None
  cq=cq//cn
  fx=fixed(cat) if E.z(cat) else 0
  if cat in ("NH4","H3O"):
    fx=1
  if cat=="Hg2":
    fx=2
  if fx and fx!=cq:
    if an=="O" and cq==2*fx and acnt%2==0:
      an="O2"
      aq=-2
      acnt=acnt//2
      cq=fx
    else:
      return None
  return (cat,cq,cn,an,aq,acnt)

def isionic(s):
  lst,q=E.parse(s)
  for e,n in lst:
    if E.metal(e):
      return True
  return s.replace("(","").startswith("NH4")

def unit(c,cq,a,aq):
  # criss-cross and reduce; returns formula string and counts
  nc=abs(aq)
  na=abs(cq)
  d=L.gcd(nc,na)
  nc//=d
  na//=d
  return fmt(c,nc,a,na),nc,na

def fmt(c,nc,a,na):
  def part(x,n):
    poly=E.z(x)==0 or x in ("Hg2",)
    if x=="Hg2":
      return "Hg2" if n==1 else "(Hg2)"+str(n)
    if poly and n>1:
      return "("+x+")"+str(n)
    return x+(str(n) if n>1 else "")
  return part(c,nc)+part(a,na)

def ion(x,q):
  # Na+, Cl-, Mg^2+, SO4^2-, NH4^+
  hat=abs(q)>1 or x[-1:].isdigit()
  return x+("^" if hat else "")+qs(q)

def qs(q):
  return (str(abs(q)) if abs(q)!=1 else "")+("+" if q>0 else "-")

ALK=("Li","Na","K","Rb","Cs")

AE=("Be","Mg","Ca","Sr","Ba")

def sol(c,a,show=1):
  # returns (soluble bool, reason)
  if a in ("F",):
    bad=("Mg","Ca","Sr","Ba","Pb")
    r=c not in bad
    why="F- mostly soluble; insoluble with alkaline earths (Mg,Ca..), Pb2+"
  elif a in ("Cl","Br","I"):
    r=c not in ("Ag","Hg2","Pb")
    why="Cl-,Br-,I- soluble except with Ag+, Hg2 2+, Pb2+"
  elif a=="NO3":
    r=True
    why="NO3- always soluble"
  elif a=="SO4":
    r=c not in ("Sr","Ba","Hg2","Pb")
    why="SO4 2- soluble except Sr2+, Ba2+, Hg2 2+, Pb2+"
  elif a=="S":
    r=c=="NH4" or c in ALK or c in AE
    why="S2- insoluble except NH4+, alkali, alkaline earths"
  elif a=="OH":
    r=c=="NH4" or c in ALK or c in ("Ca","Sr","Ba")
    why="OH- insoluble except alkali, NH4+, Ca2+, Sr2+, Ba2+"
  elif a in ("CO3","PO4"):
    r=c=="NH4" or c in ALK
    why=a+" insoluble except NH4+, alkali metals"
  elif c in ALK or c=="NH4":
    r=True
    why="alkali/NH4+ salts soluble (general rule)"
  elif a in ("CH3COO","C2H3O2","ClO4","ClO3"):
    r=True
    why=a+"- soluble (general rule, not on sheet)"
  else:
    r=False
    why="not in sheet tables: assume insoluble (check book)"
  return r,why
