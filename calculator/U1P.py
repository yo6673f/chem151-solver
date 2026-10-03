# U1P: phase diagrams
import CLIB as L
import math
out=L.out
f=L.f
g=L.g

WATER=(273.16,0.00604,647.1,218.0,[(373.15,1.0)],(263.15,0.002566),-0.0101)

CO2=(216.15,5.2,304.15,73.0,[],(195.15,1.0),0.0125)

def lnp(T,a,b):
  # line through a,b in (1/T, lnP)
  x=1/T
  x1=1/a[0]
  x2=1/b[0]
  y1=math.log(a[1])
  y2=math.log(b[1])
  return math.exp(y1+(y2-y1)*(x-x1)/(x2-x1))

class PD:
  def __init__(s,d):
    s.Tt,s.Pt,s.Tc,s.Pc,ex,s.sp,s.ms=d
    s.vp=[(s.Tt,s.Pt)]+ex+[(s.Tc,s.Pc)]
  def psat(s,T):
    if T<s.Tt:
      return lnp(T,(s.Tt,s.Pt),s.sp)
    v=s.vp
    for i in range(len(v)-1):
      if T<=v[i+1][0] or i==len(v)-2:
        return lnp(T,v[i],v[i+1])
  def tmelt(s,P):
    return s.Tt+s.ms*(P-s.Pt)
  def phase(s,T,P):
    if T>=s.Tc and P>=s.Pc:
      return "supercritical"
    if T>=s.Tc or P<s.psat(T):
      if T<s.tmelt(P) and P>s.Pt:
        return "solid"
      return "gas"
    return "solid" if T<s.tmelt(P) else "liquid"

TR={"solid>liquid":"MELTING","liquid>solid":"FREEZING","liquid>gas":"BOILING (vaporization)","gas>liquid":"CONDENSATION","solid>gas":"SUBLIMATION","gas>solid":"DEPOSITION"}

def tr(a,b):
  k=a+">"+b
  if k in TR:
    return TR[k]
  return "becomes "+b+" (no sharp boundary)"

def getpd():
  k=L.pick("Substance",["water","CO2","custom"])
  if k==1:
    return PD(WATER)
  if k==2:
    return PD(CO2)
  out("Enter diagram values:")
  Tt=L.T("Triple pt T")
  Pt=L.P("Triple pt P")
  Tc=L.T("Critical T")
  Pc=L.P("Critical P")
  ex=[]
  sp=None
  if Pt<1:
    tb=L.T("Normal b.p. (blank=skip)",1)
    if tb:
      ex=[(tb,1.0)]
  else:
    ts=L.T("Normal sublimation T (1 atm)")
    sp=(ts,1.0)
  lf=L.pick("Solid-liquid line tilts",["right (normal)","left (like water)"])
  ms=0.0125 if lf==1 else -0.0101
  if sp is None:
    # estimate sublimation curve: 1.2x steeper than vapor curve
    a=(Tt,Pt)
    b=ex[0] if ex else (Tc,Pc)
    sl=(math.log(b[1])-math.log(a[1]))/(1/b[0]-1/a[0])*1.2
    T2=Tt-10
    sp=(T2,math.exp(math.log(Pt)+sl*(1/T2-1/Tt)))
  return PD((Tt,Pt,Tc,Pc,ex,sp,ms))

def tp(T,P):
  return g(round(T-273.15,2))+" C ("+g(round(T,2))+" K), "+f(P)+" atm"

def phasept():
  L.hdr("Phase at T,P")
  d=getpd()
  out("Method: convert T,P to diagram units, plot the point; region = stable phase.")
  while 1:
    T=L.T("T",1)
    if T is None:
      break
    P=L.P("P")
    ph=d.phase(T,P)
    out("Point: "+tp(T,P))
    if T<d.Tc:
      out(" boundary P here = "+f(d.psat(T))+" atm")
    out(" -> stable phase: "+ph.upper())
    if abs(T-d.Tt)<0.5 and abs(P-d.Pt)/d.Pt<0.05:
      out(" (at triple point: all 3 coexist)")

def path():
  L.hdr("What changes if...")
  d=getpd()
  k=L.pick("Change",["T at const P (horiz)","P at const T (vert)"])
  if k==1:
    P=L.P("Constant P")
    a=L.T("Start T")
    b=L.T("End T")
    fx=lambda x:d.phase(x,P)
    xs=lambda i:a+(b-a)*i/400.0
  else:
    T=L.T("Constant T")
    a=L.P("Start P")
    b=L.P("End P")
    fx=lambda x:d.phase(T,x)
    la=math.log(a)
    lb=math.log(b)
    xs=lambda i:math.exp(la+(lb-la)*i/400.0)
  out("Move "+("horizontally" if k==1 else "vertically")+"; each line crossed = a phase change.")
  p0=fx(xs(0))
  out("Start: "+p0.upper())
  n=0
  for i in range(1,401):
    p1=fx(xs(i))
    if p1!=p0:
      lo=xs(i-1)
      hi=xs(i)
      for j in range(30):
        mi=(lo+hi)/2
        if fx(mi)==p0:
          lo=mi
        else:
          hi=mi
      n+=1
      if k==1:
        out(str(n)+") "+tr(p0,p1)+" at "+g(round(hi-273.15,2))+" C ("+g(round(hi,2))+" K)")
      else:
        out(str(n)+") "+tr(p0,p1)+" at "+f(hi)+" atm = "+f(hi*760)+" mmHg")
      p0=p1
  if n==0:
    out("No line crossed: stays "+p0)
  out("End: "+p0.upper())

def pdinfo():
  L.hdr("Phase diagram facts")
  out("Region = only stable phase; line = 2 phases coexist; triple pt = all 3; critical pt = end of liquid-gas line -> supercritical fluid beyond.")
  out("Liquid-gas line = vapor pressure curve.")
  out("Normal m.p./b.p. = where line crosses P = 1 atm.")
  out("Water: s-l line tilts LEFT (ice less dense) -> pressure can melt ice. TP 0.01 C, 0.006 atm (4.58 mmHg); CP 374 C, 218 atm.")
  out("CO2: TP -57 C, 5.2 atm -> at 1 atm it SUBLIMES at -78 C (never liquid at 1 atm). CP 31 C, 73 atm.")
  d=getpd()
  if d.Pt<1:
    for nm,P in (("normal m.p.",1.0),):
      out(nm+" = "+g(round(d.tmelt(P)-273.15,2))+" C")
    lo=d.Tt
    hi=d.Tc
    for i in range(40):
      mi=(lo+hi)/2
      if d.psat(mi)<1:
        lo=mi
      else:
        hi=mi
    out("normal b.p. = "+g(round(lo-273.15,2))+" C")
  else:
    lo=100.0
    hi=d.Tt
    for i in range(40):
      mi=(lo+hi)/2
      if d.psat(mi)<1:
        lo=mi
      else:
        hi=mi
    out("At 1 atm: sublimes at "+g(round(lo-273.15,1))+" C (no liquid)")

# ---------- vapor pressure ----------
