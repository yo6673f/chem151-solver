# U4S: ionic melting points, solubility, precipitation
import CLIB as L
import CELEM as E
import IONIC as I
out=L.out
f=L.f
g=L.g

RION={"Li":76,"Na":102,"K":138,"Rb":152,"Cs":167,"Be":45,"Mg":72,"Ca":100,"Sr":118,"Ba":135,"Al":54,"Ga":62,"Zn":74,"Ag":115,"Cd":95,"Pb":119,"Sn":118,"Ni":69,"Co":75,"Mn":83,"Cr":62,"NH4":148,"Hg2":119,
"F":133,"Cl":181,"Br":196,"I":220,"O":140,"S":184,"Se":198,"Te":221,"N":146,"P":212,"OH":137,"NO3":179,"CO3":178,"SO4":230,"PO4":238,"CN":191}

RCU={1:77,2:73}

RFE={2:78,3:65}

def rad(ion,q):
  if ion=="Cu":
    return RCU.get(q,0)
  if ion=="Fe":
    return RFE.get(q,0)
  return RION.get(ion,0)

def mprank():
  L.hdr("Rank ionic melting points")
  out("Enter formulas (NaCl, MgO..)")
  rs=[]
  while 1:
    s=L.ask("Compound (blank=done)?","")
    if s=="":
      break
    p=I.split(s)
    if p is None:
      out(" can't split into ions")
      continue
    c,cq,cn,a,aq,an=p
    r=rad(c,cq)+rad(a,aq)
    rs.append((abs(cq*aq),-r if r else 0,s,c,cq,a,aq,r))
  out("1) Ions:")
  for t in rs:
    out(" "+t[2]+": "+I.ion(t[3],t[4])+" "+I.ion(t[5],t[6])+"; |q1 q2| = "+str(t[0])+"; r ~ "+(g(t[7])+" pm" if t[7] else "?"))
  rs.sort()
  rs.reverse()
  out("2-3) Highest -> lowest T_m:")
  for t in rs:
    out(" "+t[2])
  out("Reasons:")
  for i in range(len(rs)-1):
    a,b=rs[i],rs[i+1]
    if a[0]!=b[0]:
      out(" "+a[2]+" > "+b[2]+": larger charge product ("+str(a[0])+" vs "+str(b[0])+")")
    else:
      out(" "+a[2]+" > "+b[2]+": same charges, smaller ions -> smaller r")
  out("4) Coulomb: F ~ q1 q2 / r^2")

def solq():
  L.hdr("Soluble in water?")
  s=L.ask("Ionic compound?","")
  p=I.split(s)
  if p is None:
    out("Can't split into ions")
    return
  c,cq,cn,a,aq,an=p
  ok,why=I.sol(c,a)
  out("Ions: "+I.ion(c,cq)+", "+I.ion(a,aq))
  out("Rule: "+why)
  out(s+" is "+("SOLUBLE -> strong electrolyte" if ok else "INSOLUBLE (precipitate; slightly soluble -> weak electrolyte)"))
  if ok:
    out("Dissolving:")
    out(" "+s+"(s) -> "+(str(cn)+" " if cn>1 else "")+I.ion(c,cq)+"(aq) + "+(str(an)+" " if an>1 else "")+I.ion(a,aq)+"(aq)")
  out("Small charges (+1,-1) usually soluble; both >1 usually insoluble.")

def precip():
  L.hdr("Precipitate when mixed?")
  s1=L.ask("Solution 1 compound?","")
  s2=L.ask("Solution 2 compound?","")
  p1=I.split(s1)
  p2=I.split(s2)
  if p1 is None or p2 is None:
    out("Can't split into ions")
    return
  c1,q1,n1,a1,r1,m1=p1
  c2,q2,n2,a2,r2,m2=p2
  out("1) Ions: "+I.ion(c1,q1)+", "+I.ion(a1,r1)+", "+I.ion(c2,q2)+", "+I.ion(a2,r2))
  out("2) Swap partners:")
  found=0
  for c,cq,a,aq in ((c1,q1,a2,r2),(c2,q2,a1,r1)):
    fs,nc,na=I.unit(c,cq,a,aq)
    ok,why=I.sol(c,a)
    out(" "+fs+": "+("soluble (aq)" if ok else "INSOLUBLE (s)"))
    out("  "+why)
    if not ok:
      found=1
      out("3-4) Precipitate: "+fs+"(s)")
      out(" Net ionic: "+(str(nc)+" " if nc>1 else "")+I.ion(c,cq)+"(aq) + "+(str(na)+" " if na>1 else "")+I.ion(a,aq)+"(aq) -> "+fs+"(s)")
  if not found:
    out("No precipitate: all ions stay (aq)")
  else:
    try:
      import BAL
      f1=I.unit(c1,q1,a2,r2)[0]
      f2=I.unit(c2,q2,a1,r1)[0]
      ph=["(aq)","(aq)","(s)" if not I.sol(c1,a2)[0] else "(aq)","(s)" if not I.sol(c2,a1)[0] else "(aq)"]
      co=BAL.solve(s1+"+"+s2+">"+f1+"+"+f2)
      if co:
        sp=[s1,s2,f1,f2]
        st=[(str(co[i])+" " if co[i]>1 else "")+sp[i]+ph[i] for i in range(4)]
        out("Molecular eq: "+st[0]+" + "+st[1]+" -> "+st[2]+" + "+st[3])
    except Exception as ex:
      pass
