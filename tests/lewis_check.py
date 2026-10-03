import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'calculator'))
import LSOLV, LGEO as LW
cases=[("CO2",0),("H2O",0),("NH3",0),("CH4",0),("NO3",-1),("CO3",-2),("SO4",-2),("O3",0),("SO2",0),("SO3",0),
("HCN",0),("CH2O",0),("BF3",0),("PCl5",0),("SF6",0),("SF4",0),("XeF4",0),("ClF3",0),("I3",-1),("NO2",0),("NO",0),
("N2",0),("O2",0),("CO",0),("NH4",1),("H3O",1),("H2SO4",0),("HNO3",0),("H3PO4",0),("N2O",0),("H2O2",0),("N2H4",0),
("C2H4",0),("C2H2",0),("C2H6",0),("CH3OH",0),("CH3CH2OH",0),("CH3COOH",0),("HCOOH",0),("CH3COCH3",0),("CH3OCH3",0),
("CH3CHO",0),("CH2CHCH3",0),("CH3CN",0),("CH3NH2",0),("CH2Cl2",0),("CHCl3",0),("CCl4",0),("C3H8",0),("C2H5OH",0),
("CH3CO2H",0),("CH3COO",-1),("(CH3)3CH",0),("CH3(CH2)4CH3",0),("CH2CHCHCH2",0),("C4H10",0),("HOCl",0),("OF2",0),("C3H4",0),("CH3COOCH3",0),("CH3CONH2",0),("H2S",0),("PH3",0),("HCl",0),("Cl2",0),("CS2",0),("NO2",-1),("ClO4",-1),("XeF2",0),("BeCl2",0),("OCl2",0),("SCN",-1),("NOCl",0),("POCl3",0),("SOCl2",0),("COCl2",0),("NH2OH",0),("BrF5",0),("ICl4",-1),("CN",-1),("C2O4",-2)]
for s,q in cases:
    m=LSOLV.build(s,q)
    if not m:
        print(s,q,"FAIL"); continue
    lab=LW.labels(m["el"])
    bonds=[]
    for k,(a,b) in enumerate(m["ed"]):
        if m["el"][a]=="H" or m["el"][b]=="H": continue
        bonds.append(lab[a]+"-="[0:1]*0+("-","=","#")[m["bo"][k]-1]+lab[b])
    lps={lab[i]:m["lp"][i] for i in range(len(lab)) if m["lp"][i]}
    pol,why=LW.polar(m)
    geo=[]
    for i in LW.centers(m):
        d,l=LW.dom(m,i)
        geo.append(lab[i]+":"+LW.GEO.get((d,l),("?","?"))[1])
    print(s,q,"|",m["tag"],"|"," ".join(bonds),"| lp",lps,"| res",m["nres"],"| sc",m["sc"],"|",",".join(geo),"| polar" if pol else "| nonpolar")
