# MicroPython/CPython runner (cwd = calculator/): micropython ../tests/mp_sim.py MODULE < inputs.txt
import sys
sys.path.insert(0, ".")
import CLIB
CLIB.H = 100000
def _inp(p=""):
    l = sys.stdin.readline()
    if not l:
        raise EOFError()
    l = l.rstrip("\n")
    print(p + l)
    return l
CLIB.inp = _inp
m = __import__(sys.argv[1])
try:
    m.run()
except EOFError:
    pass
import gc
gc.collect()
try:
    print("\n[free heap]", gc.mem_free(), "alloc", gc.mem_alloc())
except AttributeError:
    print("\n[free heap] n/a")
