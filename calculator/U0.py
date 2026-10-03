# U0: constants and conversions
import CLIB as L
out=L.out
f=L.f
g=L.g

def consts():
  L.hdr("Constants")
  for s in ("N_A = 6.022E23 /mol","k_B = 1.380E-23 J/K",
    "R = 8.314 J/(K mol)","R = 0.08206 L atm/(K mol)",
    "R = 1.986 cal/(K mol)","h = 6.626E-34 J s",
    "c = 3.00E8 m/s","e = 1.602E-19 C",
    "V_m(STP) = 22.41 L/mol","STP = 273.15 K, 1 atm",
    "K = 8.988E9 N m2/C2",
    "1 atm=101325 Pa=760 mmHg","  =760 torr; 1 Pa=1 N/m2",
    "1 cal = 4.184 J; 1 kJ=1000 J","1 L=1000 mL=1000 cm3",
    "1 m3 = 1000 L","1 g=1000 mg; 1 mg=1000 ug",
    "1 kg = 1000 g","T(K) = T(C) + 273.15",
    "ppmv x 1E-6 = mole fraction"):
    out(s)

def fac(fc,u,b):
  # conversion factor written as fraction: (? base / ? unit)
  if fc>=1:
    return "("+g(fc)+" "+b+"/1 "+u+")"
  return "(1 "+b+"/"+g(1/fc)+" "+u+")"

def conv(tab,nm):
  L.hdr(nm+" conversion")
  units=[t[0] for t in tab]
  v=L.num("Value")
  a=tab[L.pick("From",units)-1]
  b=tab[L.pick("To",units)-1]
  base=tab[0][0]
  out("Factor-label (old unit on bottom):")
  x=v*a[1]
  if a[0]!=base:
    out(" "+g(v)+" "+a[0]+" x "+fac(a[1],a[0],base))
    out("  = "+f(x)+" "+base)
  y=x/b[1]
  if b[0]!=base:
    out(" "+f(x)+" "+base+" x "+fac(1/b[1],base,b[0]))
    out("  = "+f(y)+" "+b[0])
  out("Answer: "+f(y)+" "+b[0])

def temp():
  L.hdr("Temperature")
  v=L.num("Value")
  k=L.pick("Unit of value",["C","K","F"])
  if k==1:
    K=v+273.15
    out("K = C + 273.15")
    out("  = "+g(v)+" + 273.15 = "+g(K)+" K")
  elif k==2:
    K=v
  else:
    K=(v-32)*5/9.0+273.15
    out("C = (F-32) x 5/9 = "+f(K-273.15,5))
    out("K = C + 273.15 = "+f(K,5))
  C=K-273.15
  out("T = "+g(C)+" C")
  out("T = "+g(K)+" K")
  out("T = "+g(C*9/5.0+32)+" F")
  if K<0:
    out("Warning: below absolute zero!")

PRE=(("G",9),("M",6),("k",3),("h",2),("(none)",0),("c",-2),("m",-3),("u",-6),("n",-9),("p",-12),("f",-15))
def prefix():
  L.hdr("Metric prefixes")
  out("G=1E9 M=1E6 k=1E3 h=1E2")
  out("c=1E-2 m=1E-3 u=1E-6")
  out("n=1E-9 p=1E-12 f=1E-15")
  v=L.num("Value")
  nm=[p[0] for p in PRE]
  a=PRE[L.pick("From",nm)-1]
  b=PRE[L.pick("To",nm)-1]
  u=L.ask("Base unit (g,m,L,J..)?","unit")
  pa=a[0] if a[1] else ""
  pb=b[0] if b[1] else ""
  x=v*10.0**a[1]
  out(" "+g(v)+" "+pa+u+" x (1E"+str(a[1])+" "+u+"/1 "+pa+u+")")
  out("  = "+f(x)+" "+u)
  y=x/10.0**b[1]
  out(" "+f(x)+" "+u+" x (1 "+pb+u+"/1E"+str(b[1])+" "+u+")")
  out("  = "+f(y)+" "+pb+u)

def ppm():
  L.hdr("ppmv <-> mole fraction")
  k=L.pick("Given",["ppmv","mole fraction","%"])
  v=L.num("Value")
  if k==1:
    x=v*1e-6
    out("x = ppmv x 1E-6 = "+f(x))
  elif k==2:
    x=v
  else:
    x=v/100
    out("x = %/100 = "+f(x))
  out("mole fraction = "+f(x))
  out("ppmv = x x 1E6 = "+f(x*1e6))
  out("percent = x x 100 = "+f(x*100)+" %")
  out("ppbv = x x 1E9 = "+f(x*1e9))

def chain():
  L.hdr("Factor-label chain")
  out("Enter start value, then each")
  out("factor as top/bottom. Blank")
  out("top value = finish.")
  v=L.num("Start value")
  u=L.ask("Start unit?","")
  while 1:
    t=L.numu("Factor top value")
    if t is None:
      break
    tu=L.ask(" top unit?","")
    b=L.num("Factor bottom value",1)
    bu=L.ask(" bottom unit?",u)
    if bu!=u:
      out(" Note: bottom unit should be "+u+" so it cancels!")
    nv=v*t/b
    out(" "+f(v)+" "+u+" x ("+f(t)+" "+tu+"/"+f(b)+" "+bu+")")
    out("  = "+f(nv)+" "+tu)
    v=nv
    u=tu
  out("Answer: "+f(v)+" "+u)
  out("Check units left = wanted.")

def ml():
  if L.pick("Convert",["mass","length"])==1:
    conv(L.UM,"Mass")
  else:
    conv(L.UL,"Length")

def run():
  L.run_menu("0 Constants/Conversions",
    ["Constants table","Temperature C/K/F","Pressure","Energy","Volume","Mass / length","Metric prefix","ppmv/mole fraction","Factor-label chain"],
    [consts,temp,lambda:conv(L.UP,"Pressure"),lambda:conv(L.UE,"Energy"),lambda:conv(L.UV,"Volume"),ml,prefix,ppm,chain])
