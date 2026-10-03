"""Find the smallest MicroPython heap that runs heavy menu paths of CHEM.py.
usage: python3 tests/heap_paths.py path/to/micropython"""
import subprocess, sys, os
CALC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "calculator")
PATHS = {
 "rank b.p. (U4C+Lewis)": (["5","1","4","CH3CH2CH2CH2CH3","(CH3)2CHCH2CH3","CH3OCH3","CH3CH2OH","C3H8","H2O",""], "Lowest"),
 "Lewis CH3COOH": (["4","2","1","CH3COOH",""], "TRIGONAL PLANAR"),
 "bond energy auto": (["6","3","2","CH4+O2>CO2+H2O"], "dH = 2652"),
 "precipitation": (["5","2","7","Pb(NO3)2","KI"], "Molecular eq"),
 "classify+name": (["3","1","1","O2, N2, CO2, NaCl"], "name: sodium"),
 "phase path": (["2","1","4","1","2","-5C","0.001","","500",""], "MELTING"),
 "e- config": (["4","1","5","Fe","3","1"], "Unpaired"),
 "IR of molecule": (["4","2","6","CH3COCH3",""], "Expected IR"),
 "pollutant": (["3","2","4","2","48","25C","1","","5","0.075","0"], "ppmv = x"),
 "limiting reactant": (["6","2","2","N2+H2>NH3","1","28","","1","10","","34","1"], "% yield"),
}
mp = sys.argv[1]
for name, (keys, want) in PATHS.items():
    data = ("\n".join(keys) + "\n" * 80 + "0\n0\n0\n0\n").encode()
    best = None
    ok169 = False
    for h in (169, 150, 130, 120, 110, 100, 90, 80, 70, 60, 50):
        p = subprocess.run([mp, "-X", "heapsize=%dk" % h, "-c", "import sys;sys.path.insert(0,'.');import CLIB;CLIB.H=9999;import CHEM"], input=data, capture_output=True, cwd=CALC, timeout=120)
        o = (p.stdout + p.stderr).decode(errors="replace")
        if h == 169:
            ok169 = want in o and "Traceback" not in o
        if "emory" in o or want not in o:
            break
        best = h
    print("%-24s works@169k=%s  min heap=%sk" % (name, ok169, best))
