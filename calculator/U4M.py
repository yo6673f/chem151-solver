# U4M: IMFs in a substance / between two substances
import CLIB as L
import CELEM as E
from U4C import classify,desc,lvl
out=L.out
f=L.f
g=L.g

def imfs():
  L.hdr("IMFs in a substance")
  out("Type formula; condensed for organics (CH3CH2OH).")
  E.tipcase()
  s=L.ask("Formula?","")
  r=classify(s)
  out(desc(r))
  if r["ion"]:
    out("Ionic: ion-ion attractions (400-4000 kJ/mol to separate)")
    return
  out("IMFs present:")
  out(" dispersion (all molecules) 0.05-40 kJ/mol; grows with # e-, size, elongated shape")
  if r["pol"]:
    out(" dipole-dipole (polar) 5-25")
  if r["don"]:
    out(" hydrogen bonding (H directly on N/O/F + N/O/F lone pair) 10-40")
  elif r["acc"]:
    out(" (can accept H-bonds from other molecules)")
  out("Dispersion usually the BIGGEST contributor.")
  out("Stronger IMFs -> higher T_b, T_m, dH vap, viscosity, surface tension; lower VP.")

def between():
  L.hdr("IMFs between 2 substances")
  E.tipcase()
  a=classify(L.ask("Substance A?",""))
  b=classify(L.ask("Substance B?",""))
  out(desc(a))
  out(desc(b))
  out("Interactions A-B:")
  out(" dispersion (always) 0.05-40")
  ia,ib=a["ion"],b["ion"]
  pa,pb=a["pol"] or ia,b["pol"] or ib
  if ia and ib:
    out(" ion-ion")
  elif (ia and b["pol"]) or (ib and a["pol"]):
    out(" ion-dipole 4-600 (e.g. salt in water)")
  elif ia or ib:
    out(" ion-induced dipole 3-15")
  elif a["pol"] and b["pol"]:
    out(" dipole-dipole 5-25")
  elif a["pol"] or b["pol"]:
    out(" dipole-induced dipole 2-10")
  if (a["don"] and b["acc"]) or (b["don"] and a["acc"]):
    out(" hydrogen bonding 10-40")
  # mixing
  out("Mixing: 'like dissolves like'")
  la,lb=lvl(a),lvl(b)
  if (la>=1)==(lb>=1) and abs(la-lb)<=1:
    out(" A-B interactions ~ as strong as A-A and B-B -> likely MIX (miscible/soluble)")
    for r in (a,b):
      if r["don"] and r["C"]>=5:
        out(" but "+r["s"]+" has a big hydrocarbon part -> limited solubility in water")
  else:
    out(" one substance's own interactions much stronger (e.g. water H-bonds would be lost) -> IMMISCIBLE / insoluble")
  out("Mixing always increases configurations; fails if E_p rises too much.")
