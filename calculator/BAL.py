# BAL: equation parser and balancer (integer Gauss-Jordan)
import CELEM as E
import CLIB as L

def sides(eq):
  s=eq.replace(" ","").replace("->",">").replace("=",">")
  if ">" not in s:
    raise ValueError("use > between sides")
  a,b=s.split(">",1)
  r=[_strip(x) for x in a.split("+") if x]
  p=[_strip(x) for x in b.split("+") if x]
  return r,p

def _strip(x):
  # drop a typed coefficient and (s)/(aq)/(l)/(g)
  for t in ("(s)","(l)","(g)","(aq)"):
    x=x.replace(t,"")
  i=0
  while i<len(x) and x[i].isdigit():
    i+=1
  return x[i:]

def _els(sp):
  els=[]
  for s in sp:
    for e,n in E.parse(s)[0]:
      if e not in els:
        els.append(e)
  return els

def solve(eq):
  r,p=sides(eq)
  sp=r+p
  cs=[E.counts(s) for s in sp]
  els=_els(sp)
  n=len(sp)
  A=[[cs[j].get(e,0)*(1 if j<len(r) else -1) for j in range(n)] for e in els]
  piv=[]
  row=0
  for col in range(n):
    pr=None
    for i in range(row,len(A)):
      if A[i][col]:
        pr=i
        break
    if pr is None:
      continue
    A[row],A[pr]=A[pr],A[row]
    for i in range(len(A)):
      if i!=row and A[i][col]:
        a=A[row][col]
        b=A[i][col]
        A[i]=[A[i][k]*a-A[row][k]*b for k in range(n)]
        gg=0
        for v in A[i]:
          gg=L.gcd(gg,v)
        if gg>1:
          A[i]=[v//gg for v in A[i]]
    piv.append(col)
    row+=1
    if row==len(A):
      break
  free=[c for c in range(n) if c not in piv]
  if not free:
    return None
  x=[0]*n
  m=1
  for i in range(len(piv)):
    m=m*abs(A[i][piv[i]])//L.gcd(m,abs(A[i][piv[i]]))
  for fc in free:
    x[fc]=m
  for i in range(len(piv)):
    s=0
    for fc in free:
      s+=A[i][fc]*x[fc]
    x[piv[i]]=-s//A[i][piv[i]]
  gg=0
  for v in x:
    gg=L.gcd(gg,v)
  x=[v//gg for v in x]
  if x[0]<0:
    x=[-v for v in x]
  for v in x:
    if v<=0:
      return None
  return x

def show(eq):
  out=L.out
  r,p=sides(eq)
  sp=r+p
  co=solve(eq)
  out("1) Formulas (never change subscripts):")
  out(" "+" + ".join(r)+" -> "+" + ".join(p))
  if co is None:
    out("Could not balance (check formulas).")
    return None
  let="abcdefghij"
  cs=[E.counts(s) for s in sp]
  els=_els(sp)
  out("2) Coefficients "+", ".join([let[i]+"="+sp[i] for i in range(len(sp))]))
  out(" Atom balance (left = right):")
  for e in els:
    lt=[]
    rt=[]
    for j in range(len(sp)):
      k=cs[j].get(e,0)
      if k:
        t=(str(k) if k>1 else "")+let[j]
        (lt if j<len(r) else rt).append(t)
    out(" "+e+": "+"+".join(lt)+" = "+"+".join(rt))
  out("3) Set a = 1 and solve:")
  for j in range(len(sp)):
    gg=L.gcd(co[j],co[0])
    nu=co[j]//gg
    de=co[0]//gg
    out(" "+let[j]+" = "+str(nu)+("/"+str(de) if de>1 else ""))
  if co[0]>1:
    out(" fraction -> multiply all by "+str(co[0]))
  out("4) Smallest whole numbers:")
  st=[(str(co[j])+" " if co[j]>1 else "")+sp[j] for j in range(len(sp))]
  out(" "+" + ".join(st[:len(r)])+" -> "+" + ".join(st[len(r):]))
  out(" Check each atom:")
  for e in els:
    a=sum([co[j]*cs[j].get(e,0) for j in range(len(r))])
    b=sum([co[j]*cs[j].get(e,0) for j in range(len(r),len(sp))])
    out("  "+e+": "+str(a)+" = "+str(b)+(" OK" if a==b else " !!"))
  return sp,co,len(r)
