# U5A: chemical vs physical, energy diagrams, extent, balancing
import CLIB as L
out=L.out
f=L.f
g=L.g

CHEM=("burn","combust","rust","cook","bake","rot","digest","ferment","explo","tarnish","corrod","react","decompos","electroly","spoil","photosynth","respir","fire","neutraliz","precipitat")
PHYS=("boil","melt","freez","dissolv","condens","evapor","cut","crush","mix","sublim","deposit","grind","bend","break","tear","filter","distill","magnet","stretch","fold")

def change():
  L.hdr("Chemical or physical?")
  while 1:
    s=L.ask("Process (blank=done)?","").lower()
    if s=="":
      break
    if [1 for k in CHEM if k in s]:
      out(" CHEMICAL: new substance, different composition (atoms rearranged, bonds broken/formed)")
    elif [1 for k in PHYS if k in s]:
      out(" PHYSICAL: same substance (only spacing/arrangement changes)")
    else:
      out(" Ask: is a NEW substance formed? yes = chemical")
  out("Reactions: mass conserved; fixed ratios; energy always transferred; rates vary; some complete, some equilibrium.")

def assume():
  L.hdr("Six reaction assumptions")
  for t in ("1) Atoms rearranged (bonds broken & formed); atoms conserved.",
    "2) Electrons redistribute -> E_p changes -> energy absorbed or released.",
    "3) Particles must collide; more collisions -> faster.",
    "4) Right orientation needed; simpler molecules / fewer particles -> faster.",
    "5) Need at least E_a to reach the transition state.",
    "6) Reactions run forward & backward; equilibrium when rates equal.",
    "Faster: higher T (more above E_a), higher conc./P (more collisions), lower E_a (catalyst = alternate path), simpler particles."):
    out(t)

def diag(Er,Et,Ep):
  # ASCII energy diagram, top = high E
  hi=max(Er,Et,Ep)
  lo=min(Er,Et,Ep)
  rows=7
  wd=21
  sp=(hi-lo) or 1
  def rw(E):
    return int(round((hi-E)/sp*(rows-1)))
  pts=[]
  for c in range(wd):
    if c<5:
      E=Er
    elif c<10:
      E=Er+(Et-Er)*(c-4)/5.0
    elif c==10:
      E=Et
    elif c<16:
      E=Et+(Ep-Et)*(c-10)/5.0
    else:
      E=Ep
    pts.append(rw(E))
  g2=[[" "]*wd for i in range(rows)]
  for c in range(wd):
    r=pts[c]
    if c<5 or c>15:
      ch="_"
    elif c==10:
      ch="^"
    elif c<10:
      ch="/" if Et>Er else "\\"
    else:
      ch="\\" if Et>Ep else "/"
    g2[r][c]=ch
    if 4<c<16 and c!=10:
      p=pts[c-1]
      for rr in range(min(p,r)+1,max(p,r)):
        g2[rr][c]=ch
  lab={rw(Er):"R",rw(Et):"TS",rw(Ep):"P"}
  out("Ep")
  for i in range(rows):
    out("|"+"".join(g2[i])+" "+lab.get(i,""))
  out("+"+"-"*wd+" path")

def energy():
  L.hdr("Energy diagram")
  out("dH = Ep(prod) - Ep(react) = Ea(fwd) - Ea(rev). Units kJ/mol. Blank = unknown.")
  k=L.pick("Given",["dH, Ea fwd/rev","Ep values (R, TS, P)"])
  if k==1:
    dH=L.numu("dH rxn")
    ef=L.numu("Ea forward")
    er=L.numu("Ea reverse")
    if not L.one([dH,ef,er]):
      return
    if dH is None:
      dH=ef-er
      out("dH = Ea(f) - Ea(r) = "+g(ef)+" - "+g(er)+" = "+g(dH))
    elif ef is None:
      ef=dH+er
      out("Ea(f) = dH + Ea(r) = "+g(dH)+" + "+g(er)+" = "+g(ef))
    else:
      er=ef-dH
      out("Ea(r) = Ea(f) - dH = "+g(ef)+" - ("+g(dH)+") = "+g(er))
    Er=0.0
  else:
    Er=L.num("Ep reactants")
    Et=L.num("Ep transition state")
    Ep=L.num("Ep products")
    dH=Ep-Er
    ef=Et-Er
    er=Et-Ep
    out("dH = Ep(P) - Ep(R) = "+g(Ep)+" - "+g(Er)+" = "+g(dH))
    out("Ea(f) = TS - R = "+g(ef))
    out("Ea(r) = TS - P = "+g(er))
  out("dH = "+g(dH)+", Ea(f) = "+g(ef)+", Ea(r) = "+g(er))
  if dH<0:
    out("EXOTHERMIC: products lower, surroundings warm, Ea(f) < Ea(r) -> tends to be product-favored")
  elif dH>0:
    out("ENDOTHERMIC: products higher, surroundings cool, Ea(f) > Ea(r) -> tends to be reactant-favored")
  out("Draw: x = reaction path, y = Ep; R left, P right, shifted by |dH|; hump height above R = Ea(f), above P = Ea(r).")
  diag(Er,Er+ef,Er+dH)

def multi():
  L.hdr("Multi-step diagram")
  n=int(L.num("# steps",2))
  E=L.num("Ep reactants",0)
  out("Enter Ea(f) and dH for each step:")
  pts=[("R",E)]
  big=(0,0)
  for i in range(n):
    ea=L.num(" step "+str(i+1)+" Ea")
    dh=L.num(" step "+str(i+1)+" dH")
    ts=E+ea
    E=E+dh
    pts.append(("TS"+str(i+1),ts))
    pts.append(("I"+str(i+1) if i<n-1 else "P",E))
    if ea>big[0]:
      big=(ea,i+1)
  out("One hump per step; intermediates in valleys:")
  for nm,v in pts:
    out(" "+nm+": "+g(v))
  out("Overall dH = "+g(E-pts[0][1]))
  hi=max([v for nm,v in pts])
  out("Highest point: "+[nm for nm,v in pts if v==hi][0]+" ("+g(hi)+")")
  out("Largest single barrier: step "+str(big[1])+" (Ea "+g(big[0])+")")

def extent():
  L.hdr("Equilibrium extent")
  out("[products]/[reactants] = P(forward)/P(backward)")
  a=L.num("% A->B per step")
  b=L.num("% B->A per step")
  n=L.num("starting mol A",1)
  r=a/b
  out("B/A = "+g(a)+"/"+g(b)+" = "+f(r))
  nA=n/(1+r)
  out("A = n/(1+B/A) = "+f(nA)+" mol")
  out("B = n - A = "+f(n-nA)+" mol")

def photo():
  L.hdr("Photo-activated reaction")
  Ea=L.num("E_a (kJ/mol)")
  lam=L.HP*L.C*L.NA/(Ea*1000)
  out("lambda_max = h c N_A / E_a")
  out(" = 6.626E-34x3.00E8x6.022E23/"+g(Ea*1000))
  out(" = "+f(lam)+" m = "+f(lam*1e9)+" nm")
  out("Longer lambda: not enough energy. (O2 split: 498.5 -> ~240 nm)")

def bal():
  L.hdr("Balance equation")
  out("Type like C3H8+O2>CO2+H2O")
  import CELEM as E
  E.tipcase()
  s=L.ask("Equation?","")
  import BAL
  r=BAL.show(s)
  out("Coefficients = mole (particle) ratios, NOT mass ratios.")

def run():
  L.run_menu("5 Reactions/Energy",
    ["Chemical vs physical","Six assumptions/rates","Energy diagram","Multi-step diagram","Equilibrium extent","lambda_max photo-react","Balance equation"],
    [change,assume,energy,multi,extent,photo,bal])
