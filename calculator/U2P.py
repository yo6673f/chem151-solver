# U2P: pollutant concentrations, ratio problems
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g
from U2B import getM

PU=("mg/m3","g/L","mol/L","molecules/mL","ppmv")

def poll():
  L.hdr("Pollutant concentrations")
  out("mg/m3 <-> g/L <-> mol/L <-> molecules/mL <-> ppmv")
  M=getM()
  T=L.T("T",)
  P=L.P("P")
  k=L.pick("Given unit",list(PU))
  v=L.num("Value")
  air=P/(L.RL*T)
  # convert to mol/L
  if k==1:
    gl=v/1e6
    out("mg/m3 -> g/L: x(1 g/1000 mg)x(1 m3/1000 L)")
    out(" = "+g(v)+"/1E6 = "+f(gl)+" g/L")
    ml=gl/M
    out("g/L -> mol/L: /M = "+f(gl)+"/"+L.mm(M)+" = "+f(ml))
  elif k==2:
    gl=v
    ml=gl/M
    out("mol/L = g/L / M = "+f(ml))
  elif k==3:
    ml=v
  elif k==4:
    ml=v*1000/L.NA
    out("mol/L = molec/mL x 1000/N_A = "+f(ml))
  else:
    x=v*1e-6
    out("mole fraction = ppmv x 1E-6 = "+f(x))
    out("mol air/L = P/RT = "+f(P)+"/(0.08206x"+g(T)+")")
    out(" = "+f(air)+" mol/L")
    ml=x*air
    out("mol/L = x x air = "+f(ml))
  gl=ml*M
  out("Results:")
  out(" g/L = mol/L x M = "+f(gl))
  out(" mg/m3 = g/L x 1E6 = "+f(gl*1e6))
  out(" mol/L = "+f(ml))
  N=ml*L.NA/1000
  out(" molec/mL = mol/L x N_A /1000")
  out("   = "+f(N))
  out(" air: P/RT = "+f(air)+" mol/L")
  x=ml/air
  out(" mole fraction = "+f(ml)+"/"+f(air)+" = "+f(x))
  out(" ppmv = x x 1E6 = "+f(x*1e6))
  if L.yes("Change T (fixed mole frac)"):
    T2=L.T("New T")
    out("mol/L scales by T1/T2 = "+g(T)+"/"+g(T2))
    out(" new mol/L = "+f(ml*T/T2)+"; mg/m3 = "+f(gl*1e6*T/T2))
    out(" (colder -> more concentrated; ppmv same)")

def ratios():
  L.hdr("Unit 1 ratio problems")
  k=L.pick("Type",["% increase","molecules B per A","GWP per molecule"])
  if k==1:
    a=L.num("old")
    b=L.num("new")
    out("% inc = (new-old)/old x 100")
    out(" = ("+g(b)+"-"+g(a)+")/"+g(a)+" x 100 = "+f((b-a)/a*100)+" %")
  elif k==2:
    out("Gas A:")
    MA=getM()
    mA=L.num("mass A (g)")
    out("Gas B:")
    MB=getM()
    mB=L.num("mass B (g)")
    nA=mA/MA
    nB=mB/MB
    out("n_B/n_A = (m_B/M_B)/(m_A/M_A)")
    out(" = ("+g(mB)+"/"+L.mm(MB)+")/("+g(mA)+"/"+L.mm(MA)+")")
    out(" = "+f(nB)+"/"+f(nA)+" = "+f(nB/nA))
    out(f(nB/nA)+" molecules B per A")
  else:
    out("Gas:")
    M=getM()
    gw=L.num("GWP per gram")
    r=gw*M/44.01
    out("GWP/molecule = GWP/g x M/M(CO2)")
    out(" = "+g(gw)+" x "+L.mm(M)+"/44.01 = "+f(r))
    out("1 molecule = "+f(r)+" CO2 molecules")
