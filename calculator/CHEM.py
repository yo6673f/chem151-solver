# CHEM: Chem 151 step-by-step solver - run this file
import CLIB as L

def go(m):
  L.free()
  try:
    mod=__import__(m)
    mod.run()
  except MemoryError:
    L.cls()
    L.out("Out of memory loading "+m+". Restart Python and try again.")
    L.pause()
  except ImportError:
    L.cls()
    L.out("Missing file "+m+".py - send all files to the calculator.")
    L.pause()
  L.free()

def sub(t,items,mods):
  while 1:
    k=L.menu(t,items)
    if k==0:
      return
    go(mods[k-1])

def settings():
  while 1:
    k=L.menu("Settings",["Sig figs ("+str(L.SF)+")","Screen width ("+str(L.W)+")","Lines per page ("+str(L.H)+")"])
    if k==0:
      return
    v=L.num("New value")
    if k==1:
      L.SF=max(2,min(8,int(v)))
    elif k==2:
      L.W=max(16,int(v))
    else:
      L.H=max(4,int(v))

def hlp():
  L.cls()
  for t in ("HOW TO USE",
    "- Type the number of a menu item; 0 or blank = back.",
    "- Formulas: 2nd letter lowercase (NaCl, Cl2). On TI: alpha alpha = lowercase.",
    "- Numbers: 6.022E23 or 6.022*10^23; math ok (2*16.00).",
    "- Units: type '750mmHg' or just 750 then pick unit (blank = default unit shown).",
    "- Blank answer = the unknown in 'solve for any' tools.",
    "- [enter] at pause shows more steps.",
    "- Organic molecules: use condensed formulas like CH3CH2OH, CH3COOH, CH2CHCH3.",
    "- Equations: C3H8+O2>CO2+H2O",
    "For studying & class activities."):
    L.out(t)
  L.pause()

def main():
  while 1:
    k=L.menu("CHEM 151 SOLVER",["Sec0 Constants/units","Sec1 Phases/particles","Sec2 Moles/composition","Sec3 Light/Lewis/shape","Sec4 IMFs/ionic","Sec5 Reactions/energy","Settings","Help"])
    if k==0:
      L.cls()
      return
    if k==1:
      go("U0")
    elif k==2:
      sub("Sec1",["1.1-1.5 Phases/separation","1.6-1.10 Particles/gas/PEC"],["U1A","U1B"])
    elif k==3:
      sub("Sec2",["2.1-2.3 Matter/moles/soln","2.4 Gas law/pollutants","2.5-2.7 Atoms/isotopes/MS","2.8 Empirical formulas"],["U2A","U2B","U2C","U2D"])
    elif k==4:
      sub("Sec3",["3.1-3.4 Light/e- config","3.5-3.8 Lewis/VSEPR/polar"],["U3A","U3B"])
    elif k==5:
      sub("Sec4",["4.1-4.3 IMFs/molecular","4.4-4.6 Ionic/solubility"],["U4A","U4B"])
    elif k==6:
      sub("Sec5",["5.1-5.3 Rxns/diagrams/balance","5.4 Stoichiometry","5.5 Bond energy/heat"],["U5A","U5B","U5C"])
    elif k==7:
      settings()
    else:
      hlp()

main()
