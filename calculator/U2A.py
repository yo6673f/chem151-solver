# U2A: classify & name matter, molar mass, mass-moles-particles
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

MOLEL=("H2","N2","O2","F2","Cl2","Br2","I2","P4","S8","O3")

def classify():
  L.hdr("Classify matter")
  E.tipcase()
  s=L.ask("Particles (comma sep)?","")
  sp=[x.strip() for x in s.split(",") if x.strip()]
  if len(sp)>1:
    out("2+ different kinds of particles NOT bonded together -> MIXTURE of:")
  for x in sp:
    lst,q=E.parse(x)
    fs=E.fstr(lst)
    if len(lst)==1:
      e=lst[0][0]
      out(fs+": ELEMENTARY SUBSTANCE (one type of atom; can't be decomposed)")
      if lst[0][1]>1:
        out(" molecular element: use molecule mass M("+fs+") = "+f(E.molar(fs,0),5))
    else:
      ion=0
      for e,n in lst:
        if E.metal(e):
          ion=1
      if fs[:3]=="NH4" or x.startswith("(NH4"):
        ion=1
      out(fs+": COMPOUND (2+ atom types bonded, fixed ratio)")
      if ion:
        out(" metal + nonmetal -> IONIC: network of cations/anions, high T_m, conducts molten/dissolved; formula unit = lowest ion ratio")
      else:
        out(" nonmetals only -> MOLECULAR: discrete molecules, poor conductor")
      try:
        import INAME
        nm=INAME.name(x,1)
        if nm:
          out(" name: "+nm)
      except Exception as ex:
        pass

def naming():
  L.hdr("Name <-> formula")
  import INAME as I
  k=L.pick("Given",["formula","name"])
  if k==1:
    E.tipcase()
    s=E.getf()
    nm=I.name(s,1)
    out("Name: "+(nm if nm else "(not a binary/ionic compound I know)"))
  else:
    out("e.g. iron(III) oxide, dinitrogen pentoxide")
    nm=L.ask("Name?","")
    fs=I.frm(nm)
    out("Formula: "+(fs if fs else "(not recognized - check spelling)"))
  out("Ionic: cation + anion root-ide; no prefixes; Roman numeral for multi-charge metals.")
  out("Molecular: left/lower element first; prefixes mono..deca; no mono on 1st.")

def molmass():
  L.hdr("Molar mass")
  E.tipcase()
  s=E.getf()
  if s:
    M=E.molar(s)
    out("Mass of one particle = M/N_A")
    out(" = "+f(M/L.NA)+" g")

def mmp():
  L.hdr("grams <-> moles <-> particles")
  E.tipcase()
  s=E.getf()
  if not s:
    return
  out("1) Molar mass (count every atom):")
  M=E.molar(s)
  k=L.pick("Given",["grams","moles","particles","liquid volume","atoms of 1 element"])
  out("2)")
  if k==1:
    m=L.M("mass")
    n=m/M
    out("n = m/M = "+g(m)+"/"+L.mm(M)+" = "+f(n)+" mol")
  elif k==2:
    n=L.num("moles")
  elif k==3:
    N=L.num("particles")
    n=N/L.NA
    out("n = N/N_A = "+f(N)+"/6.022E23")
    out("  = "+f(n)+" mol")
  elif k==4:
    v=L.V("volume")*1000
    d=L.num("density (g/mL)")
    m=v*d
    out("m = V x d = "+g(v)+" mL x "+g(d)+" g/mL")
    out("  = "+f(m)+" g")
    n=m/M
    out("n = m/M = "+f(n)+" mol")
  else:
    el=L.ask("Element?","")
    na=L.num("number of "+el+" atoms")
    sub=E.counts(s).get(el,1)
    N=na/sub
    out("particles = atoms/subscript = "+f(na)+"/"+str(sub)+" = "+f(N))
    n=N/L.NA
    out("n = N/N_A = "+f(n)+" mol")
  m=n*M
  N=n*L.NA
  out("Results:")
  out(" n = "+f(n)+" mol")
  out(" m = n x M = "+f(n)+" x "+L.mm(M)+" = "+f(m)+" g")
  out(" N = n x N_A = "+f(n)+" x 6.022E23 = "+f(N)+" particles")
  out("3) atoms of each element = N x subscript:")
  for e,c in E.parse(s)[0]:
    out(" "+e+": "+f(N)+" x "+str(c)+" = "+f(N*c))
  out("Mass of 1 particle = M/N_A = "+f(M/L.NA)+" g")

def run():
  L.run_menu("2 Matter/Moles/Solutions",
    ["Classify matter","Name <-> formula","Molar mass","g <-> mol <-> particles","Molarity [A]=n/V","Mass % / convert","Dilution M1V1=M2V2"],
    [classify,naming,molmass,mmp,"U2S.molar","U2S.mpct","U2S.dil"])
