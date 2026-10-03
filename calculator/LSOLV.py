# LSOLV: places electrons by formal-charge search; build() = best Lewis structure
import CELEM as E
HAL=("F","Cl","Br","I")
SHORT=("B","Be","Al")

def per(e):
  return E.pg(E.z(e))[0]

def solve(el,ed,q):
  n=len(el)
  ch=[E.chi(e) for e in el]
  pr=[per(e) for e in el]
  nv=[E.nval(e) for e in el]
  tot=sum(nv)-q
  nb=[[] for i in range(n)]
  for k in range(len(ed)):
    i,j=ed[k]
    nb[i].append(k)
    nb[j].append(k)
  var=[k for k in range(len(ed)) if el[ed[k][0]] not in ("H",)+HAL and el[ed[k][1]] not in ("H",)+HAL]
  bo=[1]*len(ed)
  best=[None,[]]
  def bs(i):
    return sum([bo[k] for k in nb[i]])
  def mult(i):
    for k in nb[i]:
      if bo[k]>1:
        return True
    return False
  def over(i):
    s=bs(i)
    if el[i]=="H":
      return s>1
    if mult(i) or pr[i]<=2:
      return s>4
    return s>6
  def ev():
    B=2*sum(bo)
    R=tot-B
    if R<0:
      return
    lone=[]
    for i in range(n):
      if el[i]=="H":
        lone.append(0)
      else:
        lone.append(max(0,8-2*bs(i)))
    D=sum(lone)
    rad=-1
    if R>D:
      x=R-D
      cand=[i for i in range(n) if el[i]!="H" and pr[i]>=3 and not mult(i)]
      if not cand:
        return
      cand.sort(key=lambda i:-len(nb[i]))
      lone[cand[0]]+=x
    elif R<D:
      sh=D-R
      for i in range(n):
        if el[i] in SHORT and sh>0:
          t=min(sh,lone[i])
          lone[i]-=t
          sh-=t
      if sh==1:
        c=[i for i in range(n) if el[i]!="H" and lone[i]>=1]
        if not c:
          return
        c.sort(key=lambda i:ch[i])
        lone[c[0]]-=1
        rad=c[0]
        sh=0
      if sh:
        return
    fc=[nv[i]-lone[i]-bs(i) for i in range(n)]
    sc=(sum([abs(x) for x in fc])+sum([max(0,fc[i]) for i in range(n) if el[i] in ("N","O","F")]),round(sum([fc[i]*ch[i] for i in range(n)]),2))
    st=(bo[:],lone,fc,rad)
    if best[0] is None or sc<best[0]:
      best[0]=sc
      best[1]=[st]
    elif sc==best[0]:
      best[1].append(st)
  def rec(t):
    if t==len(var):
      ev()
      return
    k=var[t]
    i,j=ed[k]
    for o in (1,2,3):
      bo[k]=o
      if over(i) or over(j):
        break
      rec(t+1)
    bo[k]=1
  rec(0)
  if best[0] is None:
    return None
  return {"el":el,"ed":ed,"q":q,"tot":tot,"sc":best[0],"all":best[1],"nb":nb}

def build(s,q=0):
  # best Lewis structure over all candidate skeletons
  import LEWIS
  cands=LEWIS.skeletons(s)
  best=None
  for sk,tag,pen in cands:
    if pen>=9 and best is not None:
      break
    r=solve(sk.el,sk.ed,q)
    if r is None:
      continue
    r["tag"]=tag
    r["key"]=(r["sc"][0]+(pen if pen<9 else 0),r["sc"][1],pen)
    if best is None or r["key"]<best["key"]:
      best=r
  if best:
    b,l,fcs,rad=best["all"][0]
    best["bo"]=b
    best["lp"]=l
    best["fc"]=fcs
    best["rad"]=rad
    best["nres"]=len(best["all"])
  return best

# ---------- geometry & polarity ----------
