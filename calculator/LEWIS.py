# LEWIS: formula tokenizer and skeleton builder (see LSOLV, LGEO)
import CELEM as E

HAL=("F","Cl","Br","I")

SHORT=("B","Be","Al")

TV={"H":1,"C":4,"N":3,"O":2,"S":2,"P":3,"Si":4,"B":3,"Se":2,"F":1,"Cl":1,"Br":1,"I":1,"Be":2,"Al":3,"Xe":0,"Kr":0}

def tval(e):
  if e in TV:
    return TV[e]
  v=E.nval(e)
  return 8-v if v>=4 else v

def per(e):
  return E.pg(E.z(e))[0]

# ---------- tokenizer for formulas/condensed formulas ----------

def toks(s):
  s=s.replace("=","").replace("#","").replace("-","").replace(" ","")
  cs=False
  for ch in s:
    if ch.islower():
      cs=True
  out=[]
  st=[out]
  i=0
  while i<len(s):
    c=s[i]
    if c=="(":
      g=[]
      st[-1].append(["G",g,1])
      st.append(g)
      i+=1
    elif c==")":
      st.pop()
      i+=1
      j=i
      while j<len(s) and s[j].isdigit():
        j+=1
      if j>i:
        st[-1][-1][2]=int(s[i:j])
      i=j
    else:
      e=E._sym(s,i,cs)
      if e is None:
        raise ValueError("bad formula")
      i+=len(e)
      j=i
      while j<len(s) and s[j].isdigit():
        j+=1
      st[-1].append(["A",e,int(s[i:j]) if j>i else 1])
      i=j
  return out

def condensed(tk):
  seen=[]
  for t in tk:
    if t[0]=="G":
      return True
    if t[1] in seen:
      return True
    seen.append(t[1])
  return False

# ---------- skeleton builders ----------

class Sk:
  def __init__(s):
    s.el=[]
    s.ed=[]
  def add(s,e):
    s.el.append(e)
    return len(s.el)-1
  def bond(s,i,j):
    s.ed.append((i,j))
  def deg(s,i):
    return sum([1 for a,b in s.ed if a==i or b==i])
  def free(s,i):
    return tval(s.el[i])-s.deg(i)

def parse_cond(tk,sk=None):
  # returns skeleton; builds chain from condensed formula
  if sk is None:
    sk=Sk()
  prev=None
  last=None
  pendH=0
  pendB=[]
  n=len(tk)
  for idx in range(n):
    t=tk[idx]
    if t[0]=="G":
      sub=t[1]
      if len(sub)==2 and sub[0][1]=="C" and sub[1][1]=="H" and sub[1][2]==2 and sub[0][2]==1 and prev is not None:
        for r in range(t[2]):
          c=sk.add("C")
          sk.bond(prev,c)
          for h in range(2):
            sk.bond(c,sk.add("H"))
          prev=c
          last=c
        continue
      for r in range(t[2]):
        s2=Sk()
        parse_cond(sub,s2)
        off=len(sk.el)
        for e in s2.el:
          sk.add(e)
        for a,b in s2.ed:
          sk.bond(a+off,b+off)
        if prev is None:
          pendB.append(off)
        else:
          sk.bond(prev,off)
      continue
    e=t[1]
    k=t[2]
    if e=="H":
      if last is None:
        pendH+=k
      else:
        for h in range(k):
          sk.bond(last,sk.add("H"))
      continue
    if e in HAL and last is not None:
      for h in range(k):
        sk.bond(last,sk.add(e))
      continue
    nh=0
    j=idx+1
    if j<n and tk[j][0]=="A" and tk[j][1]=="H":
      nh=tk[j][2]
      j+=1
    later=j<n
    ids=[]
    for r in range(k):
      a=sk.add(e)
      ids.append(a)
      branch=False
      if prev is not None:
        sk.bond(prev,a)
        hh=nh if r==k-1 else 0
        more=r<k-1 or later
        if e in ("O","S") and hh==0 and more and sk.el[prev]=="C" and sk.free(prev)>=2:
          branch=True
      if pendH:
        for h in range(pendH):
          sk.bond(a,sk.add("H"))
        pendH=0
      for b in pendB:
        sk.bond(a,b)
      pendB=[]
      last=a
      if not branch:
        prev=a
    if k>1 and nh and e=="C":
      # spread H over the chain atoms (CnH2n+1 style)
      tk[idx+1]=["A","H",0]
      left=nh
      for a in ids:
        m=min(left,max(0,sk.free(a)-(1 if a==ids[-1] and later else 0)))
        for h in range(m):
          sk.bond(a,sk.add("H"))
        left-=m
      for h in range(left):
        sk.bond(ids[-1],sk.add("H"))
  return sk

def skeletons(s):
  # returns list of (Sk, tag)
  try:
    tk=toks(s)
  except:
    return []
  if condensed(tk):
    return [(parse_cond(tk),"condensed",0)]
  cnt={}
  order=[]
  for t in tk:
    if t[1] not in cnt:
      cnt[t[1]]=0
      order.append(t[1])
    cnt[t[1]]+=t[2]
  res=[]
  tot=sum(cnt.values())
  if tot==1:
    sk=Sk()
    sk.add(order[0])
    return [(sk,"atom",0)]
  if tot==2:
    sk=Sk()
    a=sk.add(order[0])
    b=sk.add(order[-1] if len(order)>1 else order[0])
    sk.bond(a,b)
    return [(sk,"diatomic",0)]
  elig=[e for e in order if e not in ("H","F")]
  if not elig:
    elig=[order[0]]
  uq=[e for e in elig if cnt[e]==1] or elig
  hc=uq[0]
  for e in uq:
    if (E.chi(e) or 9)<(E.chi(hc) or 9):
      hc=e
  for c in elig:
    pen=0 if c==hc else 1.5
    others=[]
    for e in order:
      m=cnt[e]-(1 if e==c else 0)
      others+= [e]*m
    nH=others.count("H")
    outer=[e for e in others if e!="H"]
    nO=outer.count("O")
    for k in range(0,min(nH,nO)+1):
      sk=Sk()
      ci=sk.add(c)
      oi=[]
      for e in outer:
        j=sk.add(e)
        sk.bond(ci,j)
        if e=="O":
          oi.append(j)
      for h in range(k):
        sk.bond(oi[h],sk.add("H"))
      for h in range(nH-k):
        sk.bond(ci,sk.add("H"))
      lim=4 if per(c)<=2 else 6
      if sk.deg(ci)<=lim:
        res.append((sk,"central "+c,pen))
  # chain candidates
  bb=[]
  for e in order:
    if e!="H" and e not in HAL:
      bb+= [e]*cnt[e]
  mono=nH+sum([cnt.get(x,0) for x in HAL])
  if len(bb)>=2:
    base=Sk()
    ids=[base.add(e) for e in bb]
    for i in range(len(ids)-1):
      base.bond(ids[i],ids[i+1])
    fr=[max(0,base.free(i)) for i in ids]
    hal=[]
    for x in HAL:
      hal+= [x]*cnt.get(x,0)
    dist=[]
    def rec(i,left,cur):
      if len(dist)>40:
        return
      if i==len(ids):
        if left==0:
          dist.append(cur[:])
        return
      rest=sum(fr[i+1:])
      for m in range(min(fr[i],left),-1,-1):
        if left-m<=rest:
          cur.append(m)
          rec(i+1,left-m,cur)
          cur.pop()
    rec(0,mono,[])
    for d in dist:
      sk=Sk()
      for e in bb:
        sk.add(e)
      for a,b in base.ed:
        sk.bond(a,b)
      hl=hal[:]
      for i in range(len(ids)):
        for m in range(d[i]):
          x=hl.pop(0) if hl else "H"
          sk.bond(i,sk.add(x))
      res.append((sk,"chain",9))
  return res

# ---------- electron placement by formal-charge search ----------
