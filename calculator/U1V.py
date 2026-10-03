# U1V: vapor pressure, separations
import CLIB as L
import math
out=L.out
f=L.f
g=L.g

def vpboil():
  L.hdr("b.p. at a given pressure")
  out("Enter 2 points (T,VP) read off each VP curve.")
  Pe=L.P("External P")
  n=int(L.num("How many substances",1))
  res=[]
  for i in range(n):
    nm=L.ask("Name "+str(i+1)+"?","S"+str(i+1))
    T1=L.T(" T1")
    P1=L.P(" VP at T1")
    T2=L.T(" T2")
    P2=L.P(" VP at T2")
    # solve lnP linear in 1/T
    x=1/T1+(math.log(Pe)-math.log(P1))*(1/T2-1/T1)/(math.log(P2)-math.log(P1))
    Tb=1/x
    res.append((Tb,nm))
    out(" "+nm+" boils at "+g(round(Tb-273.15,1))+" C ("+g(round(Tb,1))+" K)")
  res.sort()
  out("Horizontal line at P="+f(Pe)+" atm; first curve hit = lowest b.p. = most volatile:")
  i=1
  for Tb,nm in res:
    out(" "+str(i)+". "+nm+" "+g(round(Tb-273.15,1))+" C")
    i+=1
  out("(interpolated as lnP linear in 1/T)")

def vprank():
  L.hdr("Rank volatility (same T)")
  n=int(L.num("How many",2))
  r=[]
  for i in range(n):
    nm=L.ask("Name "+str(i+1)+"?","S"+str(i+1))
    p=L.num(" VP (any same unit)")
    r.append((p,nm))
  r.sort()
  r.reverse()
  out("Higher VP = more volatile = weaker attractions = lower b.p.")
  out("Most -> least volatile:")
  for p,nm in r:
    out(" "+nm+" (VP "+g(p)+")")
  out("b.p. low -> high: "+", ".join([x[1] for x in r]))
  out("Boiling: VP = external P. Lower external P (altitude, Tucson ~700 mmHg) -> lower b.p.")

# ---------- separations ----------

def tech():
  L.hdr("Separation techniques")
  for s in ("Filtration: phase (solid vs fluid) -> solid from liquid/gas",
    "Crystallization: solubility (changes w/ T or conc.) -> pure solid from solution",
    "Distillation: boiling point -> liquids; most volatile out first",
    "Fractional distillation: close b.p.s -> column hot at bottom, cool at top; most volatile condenses highest",
    "Chromatography: attraction to stationary phase -> weakly attracted exit first (shorter retention time)"):
    out("- "+s)
  out("Pick the property that differs MOST.")

def getsubs():
  n=int(L.num("How many substances",3))
  r=[]
  for i in range(n):
    nm=L.ask("Name "+str(i+1)+"?","S"+str(i+1))
    tm=L.num(" T_m (C)")
    tb=L.num(" T_b (C)")
    r.append((tb,tm,nm))
  r.sort()
  return r

def design():
  L.hdr("Design a separation")
  s=getsubs()
  k=L.pick("Start as",["gas mixture (cool)","liquid mixture (heat)"])
  out("1) Components (sorted by T_b):")
  for tb,tm,nm in s:
    out("  "+nm+": T_m "+g(tm)+", T_b "+g(tb)+" C")
  tbs=[x[0] for x in s]
  tms=[x[1] for x in s]
  dtb=max(tbs)-min(tbs)
  dtm=max(tms)-min(tms)
  out("2) T_b spread "+g(dtb)+" C, T_m spread "+g(dtm)+" C -> differentiating characteristic: BOILING point (condense/boil off)")
  if dtm>dtb:
    out("  (T_m differs more: freezing out / crystallizing could also work)")
  if k==1:
    out("3) Cooling: highest T_b condenses first.")
    st=list(range(len(s)-1,0,-1))
    rem=[x[2] for x in s]
    j=1
    for i in st:
      T=(s[i][0]+s[i-1][0])/2
      out(" Step "+str(j)+": cool to ~"+g(round(T,1))+" C (between "+g(s[i-1][0])+" and "+g(s[i][0])+")")
      liq=s[i][2]
      sol=" (solid!)" if T<s[i][1] else ""
      out("  collect "+liq+" as liquid"+sol+"; gas left: "+", ".join([x[2] for x in s[:i]]))
      j+=1
    out(" Last gas remaining: "+s[0][2])
  else:
    out("3) Heating: lowest T_b boils off first.")
    j=1
    for i in range(len(s)-1):
      T=(s[i][0]+s[i+1][0])/2
      out(" Step "+str(j)+": heat to ~"+g(round(T,1))+" C (between "+g(s[i][0])+" and "+g(s[i+1][0])+")")
      out("  "+s[i][2]+" boils off (collect by condensing); remaining: "+", ".join([x[2] for x in s[i+1:]]))
      j+=1
    out(" Left in flask: "+s[-1][2])

def liqrange():
  L.hdr("Liquid over a T range?")
  s=getsubs()
  a=L.num("Low T (C)")
  b=L.num("High T (C)")
  out("Liquid from "+g(a)+" to "+g(b)+" C needs T_m < "+g(a)+" and T_b > "+g(b))
  for tb,tm,nm in s:
    ok=tm<a and tb>b
    out(" "+nm+": T_m "+g(tm)+(" <" if tm<a else " >=")+" "+g(a)+", T_b "+g(tb)+(" >" if tb>b else " <=")+" "+g(b)+" -> "+("YES" if ok else "no"))

def column():
  L.hdr("Fractional distill column")
  s=getsubs()
  out("1) Sorted by T_b:")
  for tb,tm,nm in s:
    out("  "+nm+" "+g(tb)+" C")
  k=L.pick("Trays",["one fraction per subst.","I enter tray temps"])
  if k==1:
    trays=[(s[i][0]+s[i+1][0])/2.0 for i in range(len(s)-1)]
  else:
    n=int(L.num("Number of trays",1))
    trays=[]
    for i in range(n):
      trays.append(L.num(" Tray "+str(i+1)+" T (C)"))
  trays.sort()
  out("2-3) Trays placed between consecutive b.p.s:")
  for t in trays:
    out("  tray at "+g(round(t,1))+" C")
  out("4) "+str(len(trays))+" trays (boundaries). Top exhaust: T_b below top T; bottom: T_b above bottom T.")
  out("5) Fractions:")
  lo=None
  bands=[None]+trays+[None]
  bands.reverse()
  for i in range(len(bands)-1):
    hi=bands[i]
    lo=bands[i+1]
    got=[nm for tb,tm,nm in s if (hi is None or tb<hi) and (lo is None or tb>=lo)]
    if hi is None:
      w="Bottom (liquid, T_b>"+g(round(lo,1))+")" if lo is not None else "All"
    elif lo is None:
      w="Top exhaust gas (T_b<"+g(round(hi,1))+")"
    else:
      w="Tray "+g(round(lo,1))+" C"
    out(" "+w+": "+(", ".join(got) if got else "(none)"))
  out("Hot at bottom, cool at top: most volatile rises highest.")

def sep2():
  k=L.pick("Separation",["design steps","fractional column","liquid in T range?"])
  [design,column,liqrange][k-1]()
