# U3L: photons, photoelectric effect, energy levels
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def region(lam):
  nm=lam*1e9
  if lam>1e-3:
    r="Radio/microwave: rotational transitions"
  elif nm>700:
    r="Infrared: vibrational transitions (stretch, bend)"
  elif nm>=400:
    c="red"
    for a,b in ((400,"violet"),(450,"blue"),(495,"green"),(570,"yellow"),(590,"orange"),(620,"red")):
      if nm>=a:
        c=b
    r="Visible ("+c+"): electron excitation"
  elif nm>=10:
    r="Ultraviolet: e- excitation, bond breaking, ionization"
  else:
    r="X-ray/gamma: removes core electrons"
  return r

def photon():
  L.hdr("lambda <-> nu <-> E")
  out("c = lambda nu; E = h nu = hc/lambda")
  k=L.pick("Given",["wavelength","frequency (Hz)","wavenumber cm-1","E per photon (J)","E per mole (kJ/mol)"])
  h=L.HP
  c=L.C
  if k==1:
    lam=L.LN("wavelength")
    out("1) lambda = "+f(lam)+" m")
    nu=c/lam
    out("2) nu = c/lambda = 3.00E8/"+f(lam))
    out("  = "+f(nu)+" s-1")
  elif k==2:
    nu=L.num("nu (Hz)")
    lam=c/nu
    out("lambda = c/nu = "+f(lam)+" m")
  elif k==3:
    wn=L.num("wavenumber (cm-1)")
    out("1) cm-1 x 100 = "+f(wn*100)+" m-1")
    lam=1/(wn*100)
    out(" lambda = 1/(wavenumber) = "+f(lam)+" m")
    nu=c/lam
    out(" nu = c x wavenumber = "+f(nu)+" s-1")
  elif k==4:
    Ep=L.num("E (J/photon)")
    nu=Ep/h
    out("nu = E/h = "+f(Ep)+"/6.626E-34")
    out(" = "+f(nu)+" s-1")
    lam=c/nu
    out("lambda = c/nu = hc/E = "+f(lam)+" m")
  else:
    Em=L.num("E (kJ/mol)")
    Ep=Em*1000/L.NA
    out("4) per photon E = E(J/mol)/N_A")
    out(" = "+g(Em)+"x1000/6.022E23 = "+f(Ep)+" J")
    nu=Ep/h
    lam=c/nu
    out(" lambda = hc/E = 6.626E-34x3.00E8/"+f(Ep))
    out("  = "+f(lam)+" m")
  Ep=h*nu
  out("E = h nu = 6.626E-34 x "+f(nu))
  out(" = "+f(Ep)+" J per photon")
  Em=Ep*L.NA/1000
  out("3) E_molar = N_A x h nu /1000")
  out(" = "+f(Em)+" kJ/mol")
  out("Summary:")
  out(" lambda = "+f(lam)+" m = "+f(lam*1e9)+" nm")
  out(" nu = "+f(nu)+" Hz")
  out(" wavenumber = "+f(1/(lam*100))+" cm-1")
  out(" E = "+f(Ep)+" J/photon")
  out(" E = "+f(Em)+" kJ/mol")
  out(" = "+f(Ep/L.EC)+" eV")
  out("Region: "+region(lam))
  out("Shortcut: E(kJ/mol) ~ 1.197E5/lambda(nm) = "+f(1.197e5/(lam*1e9)))

def pe():
  L.hdr("Photoelectric E_I = h nu - E_k")
  out("Leave ONE of the three blank.")
  k=L.pick("Photon given as",["wavelength","frequency","E (J)","unknown"])
  hv=None
  if k==1:
    lam=L.LN("wavelength")
    hv=L.HP*L.C/lam
    out("h nu = hc/lambda = "+f(hv)+" J")
  elif k==2:
    nu=L.num("nu (Hz)")
    hv=L.HP*nu
    out("h nu = "+f(hv)+" J")
  elif k==3:
    hv=L.num("E photon (J)")
  u=L.pick("E_I, E_k units",["J per electron","kJ/mol"])
  sc=1.0 if u==1 else 1000/L.NA
  EI=L.numu("E_I (ionization)")
  Ek=L.numu("E_k (ejected e-)")
  if not L.one([hv,EI,Ek]):
    return
  if EI is not None:
    EI*=sc
  if Ek is not None:
    Ek*=sc
  if u==2:
    out("(kJ/mol x 1000/N_A -> J)")
  if hv is None:
    hv=EI+Ek
    out("h nu = E_I + E_k = "+f(hv)+" J")
    out("lambda = hc/E = "+f(L.HP*L.C/hv)+" m")
  elif EI is None:
    EI=hv-Ek
    out("E_I = h nu - E_k = "+f(hv)+" - "+f(Ek))
    out(" = "+f(EI)+" J = "+f(EI*L.NA/1000)+" kJ/mol")
  else:
    Ek=hv-EI
    out("E_k = h nu - E_I = "+f(hv)+" - "+f(EI))
    out(" = "+f(Ek)+" J")
    if Ek<0:
      out("Negative -> h nu < threshold: NO electron ejected")
    else:
      v=(2*Ek/9.109e-31)**0.5
      out("e- speed = sqrt(2Ek/m_e) = "+f(v)+" m/s")
  if EI>0:
    out("Threshold nu0 = E_I/h = "+f(EI/L.HP)+" Hz; lambda_max = "+f(L.HP*L.C/EI*1e9)+" nm")
  out("KE grows with nu, NOT intensity (intensity = # photons).")

def levels():
  L.hdr("|dE| = E_upper - E_lower = h nu")
  u=L.pick("Level energies in",["J","kJ/mol","eV"])
  sc=[1.0,1000/L.NA,L.EC][u-1]
  n=int(L.num("# levels",2))
  E1=[]
  for i in range(n):
    E1.append(L.num(" E level "+str(i+1))*sc)
  out("Each line = one transition (absorb up, emit down; same lambda):")
  for i in range(n):
    for j in range(i+1,n):
      d=abs(E1[j]-E1[i])
      if d==0:
        continue
      lam=L.HP*L.C/d
      out(" "+str(i+1)+"<->"+str(j+1)+": dE "+f(d)+" J, lambda "+f(lam*1e9)+" nm")
      out("   "+region(lam))

def ebind():
  L.hdr("Photo-activated reaction")
  out("Longest lambda that works: photon E = E_a per molecule")
  Ea=L.num("E_a or bond E (kJ/mol)")
  lam=L.HP*L.C*L.NA/(Ea*1000)
  out("lambda_max = h c N_A / E_a")
  out(" = 6.626E-34 x 3.00E8 x 6.022E23 / "+g(Ea*1000))
  out(" = "+f(lam)+" m = "+f(lam*1e9)+" nm")
  out(region(lam))
  out("Shorter lambda (higher E) also works.")

# ---------- electron configurations ----------
