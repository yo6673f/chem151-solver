# U2C: subatomic particles, isotopes, mass spectra
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def sub():
  L.hdr("Protons/neutrons/electrons")
  out("A = Z + N_n; charge = Z - #e")
  k=L.pick("Given",["symbol, A, charge","p, n, e"])
  if k==1:
    E.tipcase()
    s=L.ask("Symbol or Z?","")
    try:
      Z=int(s)
    except:
      lst,q=E.parse(s)
      Z=E.z(lst[0][0])
    A=int(L.num("Mass number A"))
    q=int(L.num("Charge",0))
    sym=E.sym(Z)
    out("Z = "+str(Z)+" ("+E.nameZ(Z)+") -> protons = "+str(Z))
    out("N_n = A - Z = "+str(A)+" - "+str(Z)+" = "+str(A-Z))
    e=Z-q
    out("#e = Z - charge = "+str(Z)+" - ("+str(q)+") = "+str(e))
  else:
    Z=int(L.num("protons"))
    n=int(L.num("neutrons"))
    e=int(L.num("electrons"))
    sym=E.sym(Z)
    A=Z+n
    q=Z-e
    out("Z = #p = "+str(Z)+" -> "+E.nameZ(Z)+" (Z defines element)")
    out("A = Z + N_n = "+str(Z)+" + "+str(n)+" = "+str(A))
    out("charge = Z - #e = "+str(Z)+" - "+str(e)+" = "+str(q))
  qs=""
  if q:
    qs=(str(abs(q)) if abs(q)!=1 else "")+("+" if q>0 else "-")
  out("Symbol: "+sym+"-"+str(A)+" (A sup "+str(A)+", Z sub "+str(Z)+") "+qs)
  if q>0:
    out("lost e- -> CATION")
  elif q<0:
    out("gained e- -> ANION")
  out("Isotopes: same Z, different N_n.")
  out("Mass ~ A amu (p,n ~1 amu; e- ~0)")

def parts():
  L.hdr("Subatomic particles")
  out("p+: +1 (1.602E-19 C), 1.673E-27 kg, nucleus")
  out("n0: 0, 1.675E-27 kg, nucleus")
  out("e-: -1, 9.109E-31 kg, around nucleus")
  out("Nucleus ~1E-14 m vs atom ~1E-10 m: mostly empty space (Rutherford gold foil).")
  out("amu scale: C-12 = exactly 12 amu.")

def avg():
  L.hdr("Average atomic mass")
  k=L.pick("Find",["average from isotopes","% abundance (2 isotopes)"])
  if k==1:
    n=int(L.num("# isotopes",2))
    t=0.0
    tp=0.0
    big=(0,0)
    st=""
    for i in range(n):
      m=L.num(" mass "+str(i+1))
      p=L.num(" % abund "+str(i+1))
      t+=m*p
      tp+=p
      st+=("+" if i else "")+"("+g(m)+"x"+g(p)+")"
      if p>big[1]:
        big=(m,p)
    out("r.a.m. = sum(mass x %)/100")
    out(" = ["+st+"]/100")
    out(" = "+g(round(t,4))+"/100 = "+f(t/100,5)+" amu")
    if abs(tp-100)>0.5:
      out(" (% sum = "+g(tp)+", not 100!)")
    out("Check: closest to "+g(big[0])+" (most abundant).")
  else:
    m1=L.num("mass isotope 1")
    m2=L.num("mass isotope 2")
    av=L.num("average mass")
    x=(av-m2)/(m1-m2)
    out("x(m1) + (1-x)(m2) = average")
    out(" x("+g(m1)+") + (1-x)("+g(m2)+") = "+g(av))
    out(" x = (avg - m2)/(m1 - m2)")
    out("   = ("+g(av)+"-"+g(m2)+")/("+g(m1)+"-"+g(m2)+")")
    out("   = "+f(x))
    out("Isotope 1: "+f(x*100)+" %")
    out("Isotope 2: "+f((1-x)*100)+" %")

NOM={"H":1,"C":12,"N":14,"O":16,"F":19,"Si":28,"P":31,"S":32,"Cl":35,"Br":79,"I":127,"B":11,"Na":23}
LOSS=((1,"H"),(15,"CH3"),(16,"NH2 or O"),(17,"OH"),(18,"H2O"),(28,"CO or C2H4"),(29,"CHO or C2H5"),(31,"OCH3"),(35,"35Cl"),(37,"37Cl"),(43,"C3H7 or CH3CO"),(45,"COOH or OC2H5"),(79,"79Br"),(81,"81Br"),(127,"I"))

def ms():
  L.hdr("Identify MS peaks")
  out("x = m/q; highest major peak = M+ (molecular ion).")
  fm=None
  if L.yes("Know formula"):
    E.tipcase()
    fm=E.getf()
    c=E.counts(fm)
    M=0
    st=""
    for e,n0 in E.parse(fm)[0]:
      M+=NOM.get(e,int(round(E.mass(e))))*c[e]
      st+=e+":"+str(c[e])+"x"+str(NOM.get(e,int(round(E.mass(e)))))+" "
    out("M+ with common isotopes (1H,12C,14N,16O,35Cl,79Br):")
    out(" "+st)
    out(" M+ = "+str(M))
  else:
    M=L.num("M+ peak m/q")
  while 1:
    p=L.numu("Peak m/q (blank=done)")
    if p is None:
      break
    d=M-p
    if d==0:
      out(" "+g(p)+": molecular ion M+ ("+(fm or "")+"+)")
      continue
    if d<0:
      out(" "+g(p)+": M+"+g(-d)+" -> heavier isotope version (13C, 37Cl, 81Br)")
      continue
    hit=[n for m,n in LOSS if m==round(d)]
    out(" "+g(p)+": M - "+g(d)+" = lost "+(hit[0] if hit else "? (not common)"))
    if hit and fm:
      lp=hit[0].split(" ")[0].replace("35","").replace("37","").replace("79","").replace("81","")
      try:
        c=E.counts(fm)
        for e,n in E.parse(lp)[0]:
          c[e]=c.get(e,0)-n
        s="".join([e+(str(c[e]) if c[e]>1 else "") for e,n0 in E.parse(fm)[0] if c[e]>0])
        out("   ion: ["+s+"]+")
      except:
        pass
  out("Peaks at M+1, M+2: isotope versions of the same ion.")

def pat():
  L.hdr("Isotope peak pattern")
  k=L.pick("Atom",["Cl","Br","custom"])
  if k==1:
    mL,mH,fL=35,37,0.7576
  elif k==2:
    mL,mH,fL=79,81,0.5069
  else:
    mL=L.num("light mass")
    mH=L.num("heavy mass")
    fL=L.num("light % abund")/100
  fH=1-fL
  n=int(L.num("# of these atoms (X2=2)",2))
  base=L.num("mass of rest of ion",0)
  out("f_L = "+f(fL)+", f_H = "+f(fH))
  # binomial
  pr=[]
  c=1
  for j in range(n+1):
    p=c*(fL**(n-j))*(fH**j)
    pr.append(p)
    c=c*(n-j)//(j+1)
  mx=max(pr)
  c=1
  for j in range(n+1):
    m=base+(n-j)*mL+j*mH
    lab="L"*(n-j)+"H"*j
    out(" m/q "+g(m)+" ("+lab+"): "+str(c)+" x fL^"+str(n-j)+" x fH^"+str(j))
    out("   = "+f(pr[j])+" -> rel "+f(pr[j]/mx*100,3)+"%")
    c=c*(n-j)//(j+1)
  if n==2:
    out("P(LL)=fL^2, P(LH)=2 fL fH, P(HH)=fH^2")
  out("Cl: peaks 2 apart ~3:1; Br: ~1:1")

def run():
  L.run_menu("2 Atoms/Isotopes/MS",
    ["p, n, e from symbol","Particle table/notes","Average atomic mass","Identify MS peaks","Isotope peak pattern"],
    [sub,parts,avg,ms,pat])
