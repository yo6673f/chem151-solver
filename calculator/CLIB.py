# CLIB: shared screen, input, number and unit helpers
import math
W=32
H=10
SF=4
_n=0
_o=0
NA=6.022e23
KB=1.380e-23
R=8.314
RL=0.08206
RC=1.986
HP=6.626e-34
C=3.00e8
EC=1.602e-19
VM=22.41
KC=8.988e9

inp=input

def cls():
  global _n
  _n=0
  try:
    import ti_system
    ti_system.disp_clr()
  except:
    print("")

def wrap(s):
  r=[]
  for ln in str(s).split("\n"):
    if len(ln)<=W:
      r.append(ln)
      continue
    k=0
    while k<len(ln) and ln[k]==" ":
      k+=1
    ind=ln[:k]+" "
    cur=ln[:k]
    for w in ln[k:].split(" "):
      if w=="":
        continue
      t=cur+(" " if cur.strip() else "")+w
      if len(t)<=W:
        cur=t
      else:
        if cur.strip():
          r.append(cur)
        cur=ind+w
        while len(cur)>W:
          r.append(cur[:W])
          cur=ind+cur[W:]
    r.append(cur)
  return r

def out(s=""):
  global _n,_o
  for ln in wrap(s):
    if _n>=H:
      pause()
    print(ln)
    _n+=1
    _o+=1

def pause():
  global _n
  try:
    inp("  [enter]")
  except:
    pass
  _n=0

def done():
  pause()

def hdr(s):
  out("-- "+s+" --")

def ask(p,d=None):
  global _n,_o
  if _n>=H and _o:
    pause()
  s=inp(p+" ").strip()
  _n+=1
  _o=0
  if s=="" and d is not None:
    return d
  return s

ENV={"pi":math.pi,"sqrt":math.sqrt,"log":math.log10,"ln":math.log,"exp":math.exp}

def val(s):
  s=s.strip().replace(" ","")
  for a in ("x10^","*10^","X10^"):
    s=s.replace(a,"e")
  s=s.replace("E","e")
  try:
    return float(s)
  except:
    pass
  return float(eval(s.replace("^","**"),ENV))

def num(p,d=None,blank=0):
  while 1:
    s=ask(p+"?")
    if s=="":
      if d is not None:
        return d
      if blank:
        return None
      continue
    try:
      return val(s)
    except:
      out("Bad number, retry")

def numu(p):
  return num(p,None,1)

def yes(p):
  s=ask(p+" (1=y 0=n)?","0")
  return s[:1] in ("1","y","Y")

def menu(t,items):
  global _n
  cls()
  print(t[:W])
  i=1
  for it in items:
    print((str(i)+" "+it)[:W])
    i+=1
  _n=0
  s=ask("#? (0=back)","0")
  try:
    k=int(s)
  except:
    return 0
  if k<0 or k>len(items):
    return 0
  return k

def pick(t,items,sep=":"):
  s=""
  i=1
  for it in items:
    s+=str(i)+sep+it+" "
    i+=1
  if _n+len(wrap(t+" "+s))>=H and _o:
    pause()
  out(t+" "+s)
  while 1:
    v=ask("#?","1")
    try:
      k=int(v)
      if 1<=k<=len(items):
        return k
    except:
      pass

def f(x,sf=0):
  sf=sf or SF
  try:
    x=float(x)
  except:
    return str(x)
  if x!=x:
    return "?"
  if x==0:
    return "0"
  if abs(x)>1e37:
    return "inf"
  a=abs(x)
  e=int(math.floor(math.log10(a)))
  m=x/(10.0**e)
  m=round(m,sf-1)
  if abs(m)>=10:
    m=m/10
    e+=1
  if -3<=e<=5:
    d=sf-1-e
    if d>0:
      return ("%."+str(d)+"f")%(m*10.0**e)
    q=10.0**(-d)
    return "%d"%int(round(x/q)*q)
  return ("%."+str(sf-1)+"f")%m+"E"+str(e)

def g(x):
  # compact: up to 6 sig figs, no trailing zeros
  s=f(x,6)
  if "E" in s:
    m,e=s.split("E")
  else:
    m,e=s,None
  if "." in m:
    m=m.rstrip("0").rstrip(".")
  return m+("E"+e if e else "")

def mm(M):
  return "%.2f"%M

def fi(x):
  if abs(x-round(x))<1e-6:
    return str(int(round(x)))
  return f(x)

def gcd(a,b):
  a=abs(a)
  b=abs(b)
  while b:
    a,b=b,a%b
  return a

def lpad(s,n):
  s=str(s)
  return s+" "*(n-len(s))

# unit tables: (name, factor to base)
UP=(("atm",1.0),("mmHg",1/760.0),("torr",1/760.0),("Pa",1/101325.0),("kPa",1/101.325),("bar",1/1.01325))
UV=(("L",1.0),("mL",0.001),("cm3",0.001),("m3",1000.0),("dm3",1.0),("uL",1e-6))
UM=(("g",1.0),("mg",0.001),("kg",1000.0),("ug",1e-6))
UE=(("J",1.0),("kJ",1000.0),("cal",4.184),("kcal",4184.0))
UL=(("m",1.0),("nm",1e-9),("pm",1e-12),("um",1e-6),("cm",0.01),("mm",0.001),("A",1e-10))
UT=(("C",0),("K",1),("F",2))

def unitv(p,tab,blank=0,show=1):
  # read a value with unit, convert to base unit (tab[0]); show steps
  while 1:
    s=ask(p+"?")
    if s=="" and blank:
      return None
    if s=="":
      continue
    s2=s.replace(" ","")
    u=None
    for nm,fc in tab:
      if s2.lower().endswith(nm.lower()) and len(s2)>len(nm):
        cand=s2[:-len(nm)]
        try:
          v=val(cand)
          u=(nm,fc)
          break
        except:
          pass
    if u is None:
      try:
        v=val(s2)
      except:
        out("Bad number, retry")
        continue
      if len(tab)>1:
        u=tab[pick(" unit",[t[0] for t in tab],"")-1]
      else:
        u=tab[0]
    if tab is UT:
      return tconv(v,u[1],show)
    b=v*u[1]
    if show and u[0]!=tab[0][0]:
      out(" "+g(v)+" "+u[0]+" = "+f(b)+" "+tab[0][0])
    return b

def tconv(v,k,show=1):
  if k==1:
    return v
  if k==0:
    t=v+273.15
    if show:
      out(" T = "+g(v)+"+273.15 = "+g(t)+" K")
    return t
  t=(v-32)*5/9.0+273.15
  if show:
    out(" T = ("+g(v)+"-32)x5/9+273.15")
    out("   = "+g(t)+" K")
  return t

def P(p,blank=0):
  return unitv(p+" (atm)",UP,blank)
def V(p,blank=0):
  return unitv(p+" (L)",UV,blank)
def T(p,blank=0):
  return unitv(p+" (C/K)",UT,blank)
def M(p,blank=0):
  return unitv(p+" (g)",UM,blank)
def EN(p,blank=0):
  return unitv(p+" (J)",UE,blank)
def LN(p,blank=0):
  return unitv(p+" (m)",UL,blank)

def free(keep=("CLIB","CHEM")):
  try:
    import sys
    for k in list(sys.modules.keys()):
      if k.upper()==k and k not in keep and k[:2]!="TI":
        del sys.modules[k]
  except:
    pass
  try:
    import gc
    gc.collect()
  except:
    pass

def one(vals):
  # exactly one unknown (None) allowed
  k=0
  for v in vals:
    if v is None:
      k+=1
  if k!=1:
    out("Leave exactly ONE value blank (you left "+str(k)+").")
    return False
  return True

def call(fn):
  # fn: function, or "MODULE.func" loaded on demand to save memory
  if type(fn)==str:
    m,n=fn.split(".")
    fn=getattr(__import__(m),n)
  fn()

def run_menu(t,items,funcs):
  while 1:
    k=menu(t,items)
    if k==0:
      return
    cls()
    try:
      call(funcs[k-1])
    except KeyboardInterrupt:
      out("(stopped)")
    except EOFError:
      raise
    except MemoryError:
      free(("CLIB","CHEM","CELEM"))
      out("Out of memory - go back to the main menu and retry.")
    except Exception as e:
      out("Error: "+repr(e))
    free(("CLIB","CHEM","CELEM"))
    done()
