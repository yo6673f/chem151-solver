# U1A: Section 1.1-1.5 menu
import CLIB as L

def run():
  L.run_menu("1 Phases/Separations",
    ["Intensive/extensive","Heating/cooling curve","Phase at T,P (diagram)","Phase changes along path","Phase diagram facts","b.p. from VP curves","Rank volatility by VP","Separation techniques","Design separation"],
    ["U1C.props","U1C.curve","U1P.phasept","U1P.path","U1P.pdinfo","U1V.vpboil","U1V.vprank","U1V.tech","U1V.sep2"])
