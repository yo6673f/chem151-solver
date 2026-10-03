# CELEM: periodic table data, formula parser, molar mass
import CLIB as L
out=L.out
f=L.f
# sym name mass EN(0=none) IE1(kJ/mol) radius(pm)
_D="""H Hydrogen 1.008 2.20 1312 31
He Helium 4.003 0 2372 28
Li Lithium 6.941 0.98 520 128
Be Beryllium 9.012 1.57 899 96
B Boron 10.81 2.04 801 84
C Carbon 12.01 2.55 1086 76
N Nitrogen 14.01 3.04 1402 71
O Oxygen 16.00 3.44 1314 66
F Fluorine 19.00 3.98 1681 57
Ne Neon 20.18 0 2081 58
Na Sodium 22.99 0.93 496 166
Mg Magnesium 24.31 1.31 738 141
Al Aluminum 26.98 1.61 578 121
Si Silicon 28.09 1.90 787 111
P Phosphorus 30.97 2.19 1012 107
S Sulfur 32.07 2.58 1000 105
Cl Chlorine 35.45 3.16 1251 102
Ar Argon 39.95 0 1521 106
K Potassium 39.10 0.82 419 203
Ca Calcium 40.08 1.00 590 176
Sc Scandium 44.96 1.36 633 170
Ti Titanium 47.87 1.54 659 160
V Vanadium 50.94 1.63 651 153
Cr Chromium 52.00 1.66 653 139
Mn Manganese 54.94 1.55 717 139
Fe Iron 55.85 1.83 762 132
Co Cobalt 58.93 1.88 760 126
Ni Nickel 58.69 1.91 737 124
Cu Copper 63.55 1.90 745 132
Zn Zinc 65.38 1.65 906 122
Ga Gallium 69.72 1.81 579 122
Ge Germanium 72.63 2.01 762 120
As Arsenic 74.92 2.18 947 119
Se Selenium 78.97 2.55 941 120
Br Bromine 79.90 2.96 1140 120
Kr Krypton 83.80 3.00 1351 116
Rb Rubidium 85.47 0.82 403 220
Sr Strontium 87.62 0.95 550 195
Y Yttrium 88.91 1.22 600 190
Zr Zirconium 91.22 1.33 640 175
Nb Niobium 92.91 1.6 652 164
Mo Molybdenum 95.95 2.16 684 154
Tc Technetium 98 1.9 702 147
Ru Ruthenium 101.1 2.2 710 146
Rh Rhodium 102.9 2.28 720 142
Pd Palladium 106.4 2.20 804 139
Ag Silver 107.9 1.93 731 145
Cd Cadmium 112.4 1.69 868 144
In Indium 114.8 1.78 558 142
Sn Tin 118.7 1.96 709 139
Sb Antimony 121.8 2.05 834 139
Te Tellurium 127.6 2.1 869 138
I Iodine 126.9 2.66 1008 139
Xe Xenon 131.3 2.60 1170 140
Cs Cesium 132.9 0.79 376 244
Ba Barium 137.3 0.89 503 215
La Lanthanum 138.9 1.10 538 207
Ce Cerium 140.1 1.12 534 204
Pr Praseodymium 140.9 1.13 527 203
Nd Neodymium 144.2 1.14 533 201
Pm Promethium 145 1.13 540 199
Sm Samarium 150.4 1.17 545 198
Eu Europium 152.0 1.2 547 198
Gd Gadolinium 157.3 1.2 593 196
Tb Terbium 158.9 1.1 566 194
Dy Dysprosium 162.5 1.22 573 192
Ho Holmium 164.9 1.23 581 192
Er Erbium 167.3 1.24 589 189
Tm Thulium 168.9 1.25 597 190
Yb Ytterbium 173.0 1.1 603 187
Lu Lutetium 175.0 1.27 524 187
Hf Hafnium 178.5 1.3 659 175
Ta Tantalum 180.9 1.5 761 170
W Tungsten 183.8 2.36 770 162
Re Rhenium 186.2 1.9 760 151
Os Osmium 190.2 2.2 840 144
Ir Iridium 192.2 2.20 880 141
Pt Platinum 195.1 2.28 870 136
Au Gold 197.0 2.54 890 136
Hg Mercury 200.6 2.00 1007 132
Tl Thallium 204.4 1.62 589 145
Pb Lead 207.2 2.33 716 146
Bi Bismuth 209.0 2.02 703 148
Po Polonium 209 2.0 812 140
At Astatine 210 2.2 899 150
Rn Radon 222 0 1037 150
Fr Francium 223 0.7 393 260
Ra Radium 226 0.9 509 221
Ac Actinium 227 1.1 499 215
Th Thorium 232.0 1.3 587 206
Pa Protactinium 231.0 1.5 568 200
U Uranium 238.0 1.38 598 196"""
NEL=92
_S=" "
_i=0
while _i<len(_D):
  _j=_D.find(" ",_i)
  _S+=_D[_i:_j]+" "
  _i=_D.find("\n",_i)
  if _i<0:
    break
  _i+=1

def row(Z):
  # data line for element Z: [sym,name,mass,chi,IE,radius]
  i=0
  for k in range(Z-1):
    i=_D.find("\n",i)+1
  j=_D.find("\n",i)
  return _D[i:j if j>=0 else len(_D)].split()

def z(s):
  k=_S.find(" "+s+" ")
  if k<0 or not s:
    return 0
  return _S[:k+1].count(" ")

def sym(Z):
  return row(Z)[0]
def nameZ(Z):
  return row(Z)[1]
def mass(s):
  return float(row(z(s))[2])
def chi(s):
  return float(row(z(s))[3])
def ie(Z):
  return int(row(Z)[4])
def rad(Z):
  return int(row(Z)[5])
def name(s):
  return row(z(s))[1]

def pg(Z):
  # (period, group); f-block -> group 3
  for p,a,b in ((1,1,2),(2,3,10),(3,11,18),(4,19,36),(5,37,54),(6,55,86),(7,87,118)):
    if a<=Z<=b:
      k=Z-a+1
      if p==1:
        return (1,1 if Z==1 else 18)
      if p<=3:
        return (p,k if k<=2 else k+10)
      if p<=5:
        return (p,k)
      if k<=2:
        return (p,k)
      if k<=16:
        return (p,3)
      return (p,k-14)
  return (0,0)

NONM=("H","He","C","N","O","F","Ne","P","S","Cl","Ar","Se","Br","Kr","I","Xe","Rn")
SEMI=("B","Si","Ge","As","Sb","Te","Po","At")
def kind(s):
  if s in NONM:
    return "nonmetal"
  if s in SEMI:
    return "metalloid"
  return "metal"
def metal(s):
  return kind(s)=="metal"

def nval(s):
  # valence electrons (outer s + p) for main group
  Z=z(s)
  p,g=pg(Z)
  if s=="He":
    return 2
  if g<=2:
    return g
  if g>=13:
    return g-10
  return 2

def _sym(s,i,cs):
  # read element symbol at s[i]; cs: case-sensitive mode
  c=s[i]
  if cs:
    if i+1<len(s) and s[i+1].islower() and z(c+s[i+1]):
      return c+s[i+1]
    if z(c):
      return c
    return None
  c=c.upper()
  if i+1<len(s) and s[i+1].isalpha():
    two=c+s[i+1].lower()
    one=z(c) and (z(s[i+1].upper()) and s[i+1].upper()!="U")
    if z(two) and (not one or two in ("Si","Pb","Sn")):
      return two
  if z(c):
    return c
  return None

def parse(s):
  # returns (list of (sym,count) in first-seen order, charge)
  s=s.strip().replace(" ","")
  q=0
  for k in ("^","{"):
    if k in s:
      s,qs=s.split(k,1)
      qs=qs.replace("}","")
      q=_qp(qs)
  while s and s[-1] in "+-":
    q+=1 if s[-1]=="+" else -1
    s=s[:-1]
  cs=False
  for ch in s:
    if ch.islower():
      cs=True
  parts=s.replace("*",".").replace("\xb7",".").split(".")
  tot={}
  order=[]
  for part in parts:
    mult=1
    j=0
    while j<len(part) and part[j].isdigit():
      j+=1
    if j:
      mult=int(part[:j])
      part=part[j:]
    d=_grp(part,cs)
    for e in d[1]:
      if e not in tot:
        tot[e]=0
        order.append(e)
      tot[e]+=d[0][e]*mult
  return ([(e,tot[e]) for e in order],q)

def _qp(qs):
  sg=1
  if "-" in qs:
    sg=-1
  d=qs.replace("+","").replace("-","")
  return sg*(int(d) if d else 1)

def _grp(s,cs):
  st=[({},[])]
  i=0
  while i<len(s):
    c=s[i]
    if c in "([":
      st.append(({},[]))
      i+=1
    elif c in ")]":
      i+=1
      j=i
      while j<len(s) and s[j].isdigit():
        j+=1
      m=int(s[i:j]) if j>i else 1
      i=j
      d,o=st.pop()
      for e in o:
        _add(st[-1],e,d[e]*m)
    elif c.isalpha():
      e=_sym(s,i,cs)
      if e is None:
        raise ValueError("bad symbol near "+s[i:])
      i+=len(e)
      j=i
      while j<len(s) and s[j].isdigit():
        j+=1
      n=int(s[i:j]) if j>i else 1
      i=j
      _add(st[-1],e,n)
    else:
      raise ValueError("bad char "+c)
  return st[0]

def _add(t,e,n):
  if e not in t[0]:
    t[0][e]=0
    t[1].append(e)
  t[0][e]+=n

def fstr(lst,q=0):
  s=""
  for e,n in lst:
    s+=e+(str(n) if n!=1 else "")
  if q:
    s+=" "+(str(abs(q)) if abs(q)!=1 else "")+("+" if q>0 else "-")
  return s

def molar(s,show=1):
  lst,q=parse(s)
  tot=0.0
  if show:
    out("M("+s.strip()+"):")
  for e,n in lst:
    m=mass(e)
    tot+=n*m
    if show:
      ms=f(m)
      dm=len(ms.split(".")[1]) if "." in ms else 0
      out(" "+e+": "+str(n)+" x "+ms+" = "+("%."+str(dm)+"f")%(n*m))
  if show:
    out(" M = "+"%.2f"%tot+" g/mol")
  return tot

def getf(p="Formula"):
  while 1:
    s=L.ask(p+"?")
    if s=="":
      continue
    try:
      lst,q=parse(s)
      if not lst:
        continue
      cs=False
      for ch in s:
        if ch.islower():
          cs=True
      if not cs and fstr(lst,q).replace(" ","")!=s.replace(" ","").replace("^",""):
        out(" read as "+fstr(lst,q))
      return s
    except Exception as e:
      out(str(e))

def counts(s):
  d={}
  for e,n in parse(s)[0]:
    d[e]=n
  return d

def nelec(s):
  lst,q=parse(s)
  t=0
  for e,n in lst:
    t+=z(e)*n
  return t-q

def tipcase():
  out("Tip: type 2nd letter lowercase")
  out("(alpha alpha = lowercase).")
