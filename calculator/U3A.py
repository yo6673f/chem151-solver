# U3A: Section 3.1-3.4 menu
import CLIB as L

def run():
  L.run_menu("3 Light/Atoms",
    ["lambda<->nu<->E photon","Photoelectric effect","Energy levels -> lines","lambda_max for E_a/bond","e- configuration+PES","Identify from config/PES","Rank radius/IE/chi","Rank atom/ion sizes","Valence/radical/notes"],
    ["U3L.photon","U3L.pe","U3L.levels","U3L.ebind","U3E.config","U3E.ident","U3T.trends","U3T.ions","U3T.vr"])
