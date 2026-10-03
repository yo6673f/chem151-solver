# prints RAM used by importing each module (run in calculator/ under micropython)
import gc, sys
sys.path.insert(0, ".")
gc.collect()
base = gc.mem_alloc()
for m in ["CLIB","CELEM","IONIC","LEWIS","BAL","U0","U1A","U1B","U2A","U2B","U2C","U2D","U3A","U3B","U4A","U4B","U5A","U5B","U5C"]:
    gc.collect(); a = gc.mem_alloc()
    peak = gc.mem_free()
    __import__(m)
    gc.collect(); b = gc.mem_alloc()
    print(m, "kept", b - a)
gc.collect(); print("total", gc.mem_alloc() - base)
