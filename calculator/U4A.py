# U4A: molecular vs ionic, polymers, menu
import CLIB as L
import CELEM as E
out=L.out
f=L.f
g=L.g

def kind():
  L.hdr("Molecular or ionic?")
  E.tipcase()
  s=E.getf()
  lst,q=E.parse(s)
  mets=[e for e,n in lst if E.metal(e)]
  non=[e for e,n in lst if not E.metal(e)]
  if not mets and not s.replace("(","").startswith("NH4"):
    out("Only nonmetals -> MOLECULAR (discrete molecules held by IMFs)")
  else:
    out("Metal ("+",".join(mets)+") + nonmetal -> IONIC (network of ions, e- transferred)")
  if len(lst)==2:
    a,b=lst[0][0],lst[1][0]
    ca,cb=E.chi(a),E.chi(b)
    d=abs(ca-cb)
    dl=d/(ca+cb)
    out("dchi = |"+g(ca)+"-"+g(cb)+"| = "+f(d,3)+(" (>~2 -> ionic)" if d>2 else ""))
    out("% ionic ~ delta x 100 = "+f(dl*100,3)+"%"+(" >50% -> ionic" if dl>0.5 else ""))
  out("(NaCl ~80% ionic, delta ~ +-0.8)")

def poly():
  L.hdr("Polymers")
  for t in ("Longer chains: more contact + tangling -> higher T_m, harder, more viscous",
    "Linear chains: pack closer -> higher density, T_m (HDPE)",
    "Short branches: chains slide -> lower density (LDPE), lower viscosity",
    "Long branches: entangle -> higher viscosity",
    "Plasticizer: spaces chains apart -> softer, flexible",
    "Cross-linking: strong links between chains -> rigid, harder (vulcanized rubber, S links)",
    "Many H-bonds between chains: strong, heat resistant (nylon, Kevlar)",
    "Polymer = many monomers joined covalently; repeat unit in brackets with n",
    "Polar/H-bonding groups (-OH, -COOH, -NH2, C=O) -> hydrophilic; hydrocarbon groups -> hydrophobic",
    "Proteins: amino acids (-NH2 + -COOH + R) joined by peptide bonds (water released)",
    "1 deg = sequence; 2 deg = alpha-helix/beta-sheet (H-bonds N-H...O=C); 3/4 deg = 3D fold/assembly"):
    out("- "+t)

def run():
  L.run_menu("4 IMFs/Molecular",
    ["Molecular or ionic?","IMFs in a substance","IMFs/mixing of 2 subst.","Rank boiling points","Polymers & proteins"],
    [kind,"U4M.imfs","U4M.between","U4C.rank",poly])
