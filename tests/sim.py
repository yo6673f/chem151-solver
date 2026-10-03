# Desktop simulator: feeds scripted keypresses to a calculator module.
# usage: python3 tests/sim.py MODULE "in1" "in2" ...
import sys, os, builtins
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "calculator"))
mod = sys.argv[1]
q = list(sys.argv[2:])
class Done(BaseException): pass
def fake(p=""):
    if not q:
        raise Done()
    v = q.pop(0)
    if p.strip() != "[enter]":
        print(p + v)
    return v
builtins.input = fake
import CLIB
CLIB.H = 10**6   # no paging in tests
m = __import__(mod)
try:
    m.run()
except Done:
    pass
