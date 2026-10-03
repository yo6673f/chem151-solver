# U1B: particulate model, ideal gases (particle form), PEC
import CLIB as L
import math
out=L.out
f=L.f
g=L.g

def pmass():
  k=L.pick("Particle mass from",["molar mass g/mol","mass in kg"])
  if k==1:
    M=L.num("M (g/mol)")
    m=M/L.NA/1000
    out(" m = M/N_A/1000 = "+g(M)+"/6.022E23/1000")
    out("   = "+f(m)+" kg per particle")
    return m
  return L.num("m (kg)")

def ek():
  L.hdr("<Ek> = 1/2 m <v>^2")
  out("Leave the unknown blank.")
  k=L.pick("Given",["Ek and m","m and v","Ek and v","T and M (gas)"])
  if k==4:
    T=L.T("T")
    M=L.num("M (g/mol)")
    E=1.5*L.KB*T
    out("<Ek> = 3/2 k_B T")
    out(" = 1.5 x 1.380E-23 x "+g(T))
    out(" = "+f(E)+" J per particle")
    m=M/L.NA/1000
    out("m = M/N_A/1000 = "+f(m)+" kg")
    v=math.sqrt(2*E/m)
    out("v = sqrt(2<Ek>/m) = "+f(v)+" m/s")
    out("(= sqrt(3RT/M), M in kg/mol)")
    out("Same T -> same <Ek> for every gas/phase.")
    return
  if k==1:
    E=L.num("Ek (J)")
    m=pmass()
    v=math.sqrt(2*E/m)
    out("v = sqrt(2 Ek/m)")
    out(" = sqrt(2x"+f(E)+"/"+f(m)+")")
    out(" = "+f(v)+" m/s")
  elif k==2:
    m=pmass()
    v=L.num("v (m/s)")
    E=0.5*m*v*v
    out("Ek = 1/2 m v^2 = 0.5x"+f(m)+"x"+g(v)+"^2")
    out(" = "+f(E)+" J")
  else:
    E=L.num("Ek (J)")
    v=L.num("v (m/s)")
    m=2*E/(v*v)
    out("m = 2 Ek / v^2 = "+f(m)+" kg")
    out("M = m x N_A x 1000 = "+f(m*L.NA*1000)+" g/mol")
  out("T measures <Ek>: T=(2/3)<Ek>/k_B = "+f(2*E/3/L.KB)+" K")

def speeds():
  L.hdr("Compare speeds, same T")
  M1=L.num("M gas 1 (g/mol)")
  M2=L.num("M gas 2 (g/mol)")
  out("Same T -> same <Ek>:")
  out(" 1/2 m1 v1^2 = 1/2 m2 v2^2")
  r=math.sqrt(M2/M1)
  out(" v1/v2 = sqrt(m2/m1) = sqrt("+g(M2)+"/"+g(M1)+")")
  out("  = "+f(r))
  out("Gas "+("1" if r>1 else "2")+" (lighter) moves faster.")

def igp():
  L.hdr("P = k_B N T / V")
  out("SI: P Pa, V m3, T K. Blank = unknown.")
  P=L.P("P",1)
  N=L.numu("N (particles)")
  T=L.T("T",1)
  V=L.V("V",1)
  if not L.one([P,N,T,V]):
    return
  kb=L.KB
  if P is not None:
    Pp=P*101325
    out(" P = "+f(Pp)+" Pa")
  if V is not None:
    Vm=V/1000
    out(" V = "+f(Vm)+" m3")
  if P is None:
    Pp=kb*N*T/Vm
    out("P = k_B N T/V")
    out(" = 1.380E-23 x "+f(N)+" x "+g(T)+" / "+f(Vm))
    out(" = "+f(Pp)+" Pa = "+f(Pp/101325)+" atm")
  elif N is None:
    N=Pp*Vm/(kb*T)
    out("N = P V/(k_B T)")
    out(" = "+f(Pp)+" x "+f(Vm)+" / (1.380E-23 x "+g(T)+")")
    out(" = "+f(N)+" particles")
    out(" = "+f(N/L.NA)+" mol")
  elif T is None:
    T=Pp*Vm/(kb*N)
    out("T = P V/(k_B N) = "+f(T)+" K")
  else:
    Vm=kb*N*T/Pp
    out("V = k_B N T/P = "+f(Vm)+" m3 = "+f(Vm*1000)+" L")
  out("P ~ T, P ~ N, P ~ 1/V. P does NOT depend on particle mass.")

def comb():
  L.hdr("P1V1/T1 = P2V2/T2")
  out("N fixed. Leave ONE blank. P,V any matching units; T converted to K.")
  P1=L.numu("P1")
  V1=L.numu("V1")
  T1=L.T("T1",1)
  P2=L.numu("P2")
  V2=L.numu("V2")
  T2=L.T("T2",1)
  if not L.one([P1,V1,T1,P2,V2,T2]):
    return
  out("List: P1="+s(P1)+" V1="+s(V1)+" T1="+s(T1))
  out("      P2="+s(P2)+" V2="+s(V2)+" T2="+s(T2))
  if V2 is None:
    V2=V1*(P1/P2)*(T2/T1)
    out("V2 = V1 x (P1/P2) x (T2/T1)")
    out(" = "+g(V1)+" x ("+g(P1)+"/"+g(P2)+") x ("+g(T2)+"/"+g(T1)+")")
    out(" = "+f(V2))
  elif P2 is None:
    P2=P1*(V1/V2)*(T2/T1)
    out("P2 = P1 x (V1/V2) x (T2/T1)")
    out(" = "+g(P1)+" x ("+g(V1)+"/"+g(V2)+") x ("+g(T2)+"/"+g(T1)+")")
    out(" = "+f(P2))
  elif T2 is None:
    T2=T1*(P2*V2)/(P1*V1)
    out("T2 = T1 x (P2 V2)/(P1 V1)")
    out(" = "+f(T2)+" K = "+f(T2-273.15)+" C")
  elif V1 is None:
    V1=V2*(P2/P1)*(T1/T2)
    out("V1 = V2 x (P2/P1) x (T1/T2) = "+f(V1))
  elif P1 is None:
    P1=P2*(V2/V1)*(T1/T2)
    out("P1 = P2 x (V2/V1) x (T1/T2) = "+f(P1))
  else:
    T1=T2*(P1*V1)/(P2*V2)
    out("T1 = T2 x (P1 V1)/(P2 V2) = "+f(T1)+" K")
  out("Sanity: lower P -> bigger V; lower T -> smaller V.")
  if V2>V1:
    out(" V grew x"+f(V2/V1))
  else:
    out(" V shrank x"+f(V2/V1))

def s(x):
  return "?" if x is None else g(x)

def ratio():
  L.hdr("Proportional reasoning")
  out("P = k_B N T / V. Enter factor each changes by (1 = same).")
  N=L.num("N x",1)
  T=L.num("T(K) x",1)
  V=L.num("V x",1)
  P=N*T/V
  out("P_new/P_old = (N x T)/V")
  out(" = ("+g(N)+" x "+g(T)+")/"+g(V)+" = "+f(P))
  out("Pressure is multiplied by "+f(P))
  out("Same T,P,V -> same N for any gas (Avogadro).")

def real():
  L.hdr("Real vs ideal gases")
  k=L.pick("Conditions",["high T, low P","low T","high P","low T and high P"])
  if k==1:
    out("Ideal: particles far apart, interactions negligible.")
  else:
    out("Deviations grow near condensation (low T, high P).")
    if k in (3,4):
      out("High P: particle SIZE matters -> repulsion -> less free volume -> P HIGHER than ideal.")
    if k in (2,4):
      out("Low T: ATTRACTIONS matter -> particles hit walls less -> P LOWER than ideal.")

def pec():
  L.hdr("PEC justification")
  out("y = E_p; x = # configurations. Low T favors low E_p; high T favors many configs; high P favors fewer configs; low P favors more configs (gas).")
  a=L.ask("State w/ lower Ep (e.g. liquid)?","liquid")
  b=L.ask("State w/ more configs (e.g. gas)?","gas")
  k=L.pick("Condition",["low T","high T","high P","low P"])
  win=a if k in (1,3) else b
  fac=["E_p","configurations","fewer configurations","configurations"][k-1]
  c=["low T","high T","high P","low P"][k-1]
  out("Justification:")
  out(" "+a+" has lower E_p because particles are closer and attract more strongly. "+b+" has more configurations because particles can occupy more space. At "+c+", "+fac+" dominate(s), so "+win+" is favored.")
  out("Comparing 2 liquids: one placed LOWER on E_p has stronger attractions -> higher b.p.")

def notes():
  L.hdr("Particulate model notes")
  for t in ("1) Matter = tiny (~1 nm) particles. 2) They move constantly & randomly through empty space. 3) Attract at long range, repel at short range.",
    "<Ek> = 1/2 m <v>^2. T measures <Ek>: same T -> same <Ek>, any substance or phase.",
    "Same T: lighter particles move faster. Triple point: all 3 phases same average speed.",
    "Higher T: speed distribution shifts right and flattens.",
    "E_p = 0 at infinite separation. Attraction: E_p more negative as particles approach, E_k increases. Repulsion: E_p rises, E_k falls.",
    "Phase change at constant T: added energy goes into E_p (separating particles), not E_k -> T flat.",
    "Condensing/freezing: E_p drops, extra E_k passed to surroundings as heat (latent heat).",
    "Stronger attractions -> more energy to separate -> higher T_m and T_b.",
    "Emergent (many particles only): T_m, T_b, density, viscosity, color, hardness, malleability, conductivity.",
    "TRAPS: particles do NOT melt, expand, soften, get heavier, or change size; only spacing & energy change.",
    "Ice floats: particles FARTHER APART, not lighter/smaller.",
    "Solid particles still vibrate; can't move past one another.",
    "Bubbles in boiling water = water particles, not H2 and O2."):
    out("- "+t)

def run():
  L.run_menu("1 Particles/Gases/PEC",
    ["<Ek>=1/2mv^2 solver","Compare gas speeds","P=k_B N T/V solver","P1V1/T1=P2V2/T2","What-if ratio (P~NT/V)","Real gas deviations","PEC justification","Notes: particles, Ep"],
    [ek,speeds,igp,comb,ratio,real,pec,notes])
