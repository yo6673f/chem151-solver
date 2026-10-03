# U1C: intensive/extensive, heating & cooling curves
import CLIB as L
import math
out=L.out
f=L.f
g=L.g

INT=("melting point","boiling point","density","conductivity","solubility","color","viscosity","specific heat","temperature","pressure","molar mass","concentration","hardness","refractive index","vapor pressure","odor","freezing point","surface tension","malleability")

EXT=("mass","volume","energy","heat","moles","amount","number of particles","length","heat capacity","weight","size","area")

def props():
  L.hdr("Intensive vs extensive")
  out("Intensive: no amount dependence (good for ID).")
  out("Extensive: depends on amount (useless for ID).")
  while 1:
    s=L.ask("Property (blank=done)?","")
    if s=="":
      break
    s=s.lower()
    hit=0
    for ex in (1,0):
      for lst,kd in ((INT,1),(EXT,0)):
        for p in lst:
          if hit or not (s==p if ex else (s in p or p in s)):
            continue
          hit=1
          if not kd:
            out(" "+p+": EXTENSIVE -> not useful to identify")
          elif p in ("temperature","pressure"):
            out(" "+p+": INTENSIVE but T,P describe the whole system, not one substance -> NOT a differentiating characteristic")
          else:
            out(" "+p+": INTENSIVE -> can be a differentiating characteristic")
    if not hit:
      out(" Not in list. Ask: does it change if you double the sample? yes=extensive")
  out("Good characteristic = unique per substance + independent of sample size + measurable selectively.")

def curve():
  L.hdr("Heating/cooling curve")
  tm=L.num("T_m (C)")
  tb=L.num("T_b (C)")
  t1=L.num("Start T (C)")
  t2=L.num("End T (C)")
  m=None
  if L.yes("Calc energy (q=mc dT)"):
    m=L.num("mass (g)")
    cs=L.num("c solid J/(g C)",2.09)
    cl=L.num("c liquid J/(g C)",4.184)
    cg=L.num("c gas J/(g C)",2.01)
    hf=L.num("dH fusion (J/g)",334)
    hv=L.num("dH vap (J/g)",2260)
  heat=t2>t1
  out("Steps ("+("heating" if heat else "cooling")+"):")
  def ph(t):
    return "solid" if t<tm else ("liquid" if t<tb else "gas")
  pts=[t1]
  for b in ((tm,tb) if heat else (tb,tm)):
    if min(t1,t2)<b<max(t1,t2):
      pts.append(b)
  pts.append(t2)
  k=1
  tot=0.0
  for i in range(len(pts)-1):
    a=pts[i]
    b=pts[i+1]
    mid=(a+b)/2.0
    p=ph(mid)
    out(str(k)+") "+p+" "+("warms" if heat else "cools")+" "+g(a)+" -> "+g(b)+" C: SLOPED (T changes)")
    if m is not None:
      c={"solid":cs,"liquid":cl,"gas":cg}[p]
      q=m*c*(b-a)
      tot+=q
      out("  q=mc dT="+g(m)+"x"+g(c)+"x("+g(b)+"-"+(g(a) if a>=0 else "("+g(a)+")")+")")
      out("   ="+f(q)+" J")
    k+=1
    if i<len(pts)-2:
      t=b
      if t==tm:
        nm="melting" if heat else "freezing"
        dh=hf if m is not None else 0
      else:
        nm="boiling" if heat else "condensing"
        dh=hv if m is not None else 0
      out(str(k)+") "+nm+" at "+g(t)+" C: FLAT, T constant until all converted; dE "+(">0 (absorbed)" if heat else "<0 (released)"))
      if m is not None:
        q=m*dh*(1 if heat else -1)
        tot+=q
        out("  q=m x dH="+g(m)+"x"+g(dh)+"="+f(q)+" J")
      k+=1
  if m is not None:
    out("Total q = "+f(tot)+" J = "+f(tot/1000)+" kJ")
  out("Boiling plateau is LONGER than melting plateau.")
  out("dE vs T plot: vertical jumps at T_m and T_b.")
  out("Sketch: mark T_m,T_b on T axis; slopes for each phase; flat lines at T_m,T_b.")

# ---------- phase diagrams ----------
