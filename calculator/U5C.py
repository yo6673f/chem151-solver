# U5C: bond energies, heat of reaction, q = m c dT
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

# (a,b,order): BE kJ/mol; * = not on formula sheet
BE={("H","H",1):436,("C","H",1):414,("C","C",1):347,("C","C",2):611,("C","C",3):737,("C","Cl",1):339,
("C","N",1):305,("C","N",2):615,("C","O",1):360,("C","O",2):736,("C","S",1):259,("H","N",1):389,
("N","N",1):163,("N","N",3):946,("N","O",1):201,("N","O",2):607,("H","O",1):464,("O","O",1):142,
("O","O",2):498,("H","S",1):368}
XB={("C","N",3):891,("C","O",3):1072,("F","H",1):565,("Cl","H",1):431,("Br","H",1):366,("H","I",1):299,
("F","F",1):159,("Cl","Cl",1):243,("Br","Br",1):193,("I","I",1):151,("C","F",1):485,("Br","C",1):276,
("C","I",1):240,("N","N",2):418,("O","S",2):523,("S","S",1):266,("F","O",1):190}
BS=("-","=","#")

def key(a,b,o):
  if a>b:
    a,b=b,a
  return (a,b,o)

def nm(k):
  a,b,o=k
  if (a,b) in (("H","O"),("H","N"),("H","S")):
    a,b=b,a
  return a+BS[o-1]+b

def look(k,co2=False):
  if co2 and k==("C","O",2):
    return 799,""
  if k in BE:
    return BE[k],""
  if k in XB:
    return XB[k],"*"
  return None,""

def table():
  L.hdr("Bond energies (kJ/mol)")
  ks=list(BE.keys())
  ks.sort()
  for k in ks:
    out(" "+nm(k)+" "+str(BE[k])+("  (799 in CO2)" if k==("C","O",2) else ""))
  out("Not on sheet (*): C#N 891, C#O 1072, H-F 565, H-Cl 431, H-Br 366, Cl-Cl 243, Br-Br 193, F-F 159, C-F 485")

def getbonds(t):
  out(t+": enter bonds like C-H, C=O, N#N (blank=done)")
  tot=0
  lst=[]
  while 1:
    s=L.ask(" bond?","").replace(" ","")
    if s=="":
      break
    o=1
    for c,v in (("-",1),("=",2),("#",3),("~",3)):
      if c in s:
        a,b=s.split(c)
        o=v
        break
    else:
      out(" use - = or #")
      continue
    k=key(E.parse(a)[0][0][0],E.parse(b)[0][0][0],o)
    be,st=look(k,s.upper()=="C=O" and L.yes(" in CO2"))
    if be is None:
      be=L.num(" BE kJ/mol (not in table)")
    n=L.num(" how many (x coef)",1)
    lst.append((nm(k),n,be))
    tot+=n*be
  return lst,tot

def showsum(t,lst):
  tot=0
  s=[]
  for b,n,be in lst:
    tot+=n*be
    s.append(g(n)+"("+g(be)+")")
  out(t+" "+" + ".join(s)+" = "+g(tot))
  return tot

def manual():
  L.hdr("dH from bond energies")
  out("dH = sum BE(broken, reactants) - sum BE(formed, products)")
  a,ta=getbonds("Broken (reactants)")
  b,tb=getbonds("Formed (products)")
  finish(a,b)

def finish(a,b):
  ta=showsum("Broken:",a)
  tb=showsum("Formed:",b)
  dH=ta-tb
  out("dH = "+g(ta)+" - "+g(tb)+" = "+g(dH)+" kJ")
  out(("Negative -> EXOTHERMIC (releases energy)" if dH<0 else "Positive -> ENDOTHERMIC (absorbs energy)"))
  out("Breaking absorbs (+); forming releases (-). Bonds don't 'store' energy.")
  return dH

def auto():
  L.hdr("dH from equation (auto)")
  out("e.g. CH4+O2>CO2+H2O")
  import BAL
  import LSOLV as W
  E.tipcase()
  s=L.ask("Equation?","")
  r=BAL.show(s)
  if r is None:
    return
  sp,co,nr=r
  L.pause()
  out("1) Lewis structures -> bonds:")
  br=[]
  fo=[]
  for j in range(len(sp)):
    m=W.build(sp[j],0)
    if m is None:
      out(" can't draw "+sp[j])
      return
    cnt={}
    od=[]
    for k in range(len(m["ed"])):
      a,b=m["ed"][k]
      kk=key(m["el"][a],m["el"][b],m["bo"][k])
      if kk not in cnt:
        od.append(kk)
      cnt[kk]=cnt.get(kk,0)+1
    t=[]
    for kk in od:
      be,st=look(kk,sp[j]=="CO2")
      if be is None:
        be=L.num(" BE of "+nm(kk)+" kJ/mol")
      t.append(str(cnt[kk])+" "+nm(kk)+st)
      (br if j<nr else fo).append((nm(kk)+" in "+sp[j],cnt[kk]*co[j],be))
    out(" "+(str(co[j])+" x " if co[j]>1 else "")+sp[j]+": "+", ".join(t))
  out("2-4) bonds x coefficient:")
  dH=finish(br,fo)
  if L.yes("Per gram / amount of a reactant"):
    i=L.pick("Fuel",sp[:nr])-1
    M=E.molar(sp[i],0)
    pg=dH/(co[i]*M)
    out("dH per gram = dH/(coef x M) = "+g(dH)+"/("+str(co[i])+" x "+L.mm(M)+")")
    out(" = "+f(pg)+" kJ/g")
    out(" per mol "+sp[i]+" = "+f(dH/co[i])+" kJ/mol")
    m=L.numu("grams burned (blank=skip)")
    if m is not None:
      n=m/M
      out("mol = "+f(n)+"; q = dH x mol/coef = "+f(dH*n/co[i])+" kJ")

def epmol():
  L.hdr("Ep of molecule vs atoms")
  out("Ep = - sum BE of its bonds")
  import LSOLV as W
  E.tipcase()
  s=L.ask("Formula?","")
  m=W.build(s,0)
  if m is None:
    out("can't draw")
    return
  tot=0
  cnt={}
  od=[]
  for k in range(len(m["ed"])):
    a,b=m["ed"][k]
    kk=key(m["el"][a],m["el"][b],m["bo"][k])
    if kk not in cnt:
      od.append(kk)
    cnt[kk]=cnt.get(kk,0)+1
  for kk in od:
    be,st=look(kk,s=="CO2")
    if be is None:
      be=L.num(" BE of "+nm(kk))
    out(" "+str(cnt[kk])+" x "+nm(kk)+" ("+str(be)+") = "+str(cnt[kk]*be))
    tot+=cnt[kk]*be
  out("Ep = -"+g(tot)+" kJ/mol")
  out("Stronger bond = deeper E_p well.")

def heat():
  L.hdr("q = m c dT")
  out("Leave ONE blank. water c = 4.184 J/(g C)")
  q=L.numu("q (J)")
  m=L.numu("m (g)")
  c=L.numu("c J/(g C) (water 4.184; blank=unknown)")
  k=L.pick("dT as",["dT","T_i and T_f"])
  if k==1:
    dT=L.numu("dT")
  else:
    ti=L.numu("T_initial")
    tf=L.numu("T_final")
    dT=None if ti is None or tf is None else tf-ti
    if dT is not None:
      out("dT = T_f - T_i = "+g(tf)+" - "+g(ti)+" = "+g(dT))
  if not L.one([q,m,c,dT]):
    return
  if q is None:
    q=m*c*dT
    out("q = m c dT = "+g(m)+" x "+g(c)+" x "+g(dT))
    out(" = "+f(q)+" J = "+f(q/1000)+" kJ")
  elif m is None:
    m=q/(c*dT)
    out("m = q/(c dT) = "+f(m)+" g")
  elif dT is None:
    dT=q/(m*c)
    out("dT = q/(m c) = "+g(q)+"/("+g(m)+" x "+g(c)+") = "+f(dT))
    if k==2:
      if ti is None:
        out("T_i = T_f - dT = "+f(tf-dT))
      else:
        out("T_f = T_i + dT = "+f(ti+dT))
  else:
    c=q/(m*dT)
    out("c = q/(m dT) = "+f(c)+" J/(g C)")
  out("q>0 absorbed by sample; q<0 released.")

def fuelwater():
  L.hdr("Fuel heats water")
  E.tipcase()
  s=E.getf("Fuel formula")
  M=E.molar(s)
  k=L.pick("Heat of combustion in",["kJ/mol","kJ/g"])
  dH=abs(L.num("|dH|"))
  mf=L.M("mass fuel")
  if k==1:
    n=mf/M
    q=n*dH*1000
    out("n = "+g(mf)+"/"+L.mm(M)+" = "+f(n)+" mol")
    out("q released = n x |dH| = "+f(q/1000)+" kJ")
  else:
    q=mf*dH*1000
    out("q released = m x |dH| = "+f(q/1000)+" kJ")
  mw=L.M("mass water")
  dT=q/(mw*4.184)
  out("q water = -q rxn (all heat to water)")
  out("dT = q/(m c) = "+f(q)+"/("+g(mw)+" x 4.184)")
  out(" = "+f(dT)+" C rise")
  ti=L.numu("T initial (blank=skip)")
  if ti is not None:
    out("T final = "+f(ti+dT)+" C")

def run():
  L.run_menu("5 Bond energy / Heat",
    ["BE table","dH from equation (auto)","dH from bonds (manual)","Ep of a molecule","q = m c dT","Fuel heats water"],
    [table,auto,manual,epmol,heat,fuelwater])
