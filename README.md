# Chem 151 Step-by-Step Solver (TI-84 Evo, Python)

A menu-driven Python app for the **TI-84 Evo** that solves every problem type on the
*Chem 151 Formula Sheet — Chemical Thinking Vol. I* and **prints every step**: the
formula used, the numbers plugged in, the unit conversions, and a one-line reason
for each conclusion. It is meant for studying and in-class activities.

It also runs on any computer with Python 3: `cd calculator && python3 CHEM.py`.

## Install on the calculator

1. Open **TI Connect Evo** (or TI Connect CE) with the calculator plugged in.
2. Send **all 36 `.py` files in the `calculator/` folder**. The software converts
   them to Python files on the calculator. Every file is needed: menus load the
   tool files only when you open them, which keeps memory free.
3. On the calculator, open the **Python** app and run **`CHEM`**.

File names are uppercase and at most 8 characters, which is what the calculator
expects. Don't rename them, because the modules import each other by name.

> Memory: the app needs at most ~110 KB of Python heap on its heaviest screens.
> The Evo has ~169 KB, so there's room to spare. The older TI-84 Plus CE Python
> (~20 KB heap) is **too small** to run it.

## How to type answers

| What | How |
|---|---|
| Menus | type the item number; `0` or just [enter] goes back |
| Chemical formulas | 2nd letter lowercase: `NaCl`, `Cl2`, `Fe2(SO4)3`, `CuSO4.5H2O` (on the TI, **alpha alpha** gives lowercase). If you type all caps the app tells you how it read it. |
| Ions | `NO3-`, `NH4+`, `SO4^2-` (use `^` for charges bigger than 1) |
| Organic molecules | condensed formulas: `CH3CH2OH`, `CH3COOH`, `CH2CHCH3`, `(CH3)3CH`, `CH3(CH2)4CH3` |
| Equations | `C3H8+O2>CO2+H2O` |
| Numbers | `6.022E23`, `6.022*10^23`, or math like `2*16.00` |
| Units | type a value with its unit (`750mmHg`, `25C`, `250mL`, `500nm`) or just the number and pick the unit from the list ([enter] = the default unit shown) |
| Solve-for-anything tools | leave the unknown **blank** |
| `[enter]` | shown when the screen is full; press enter for more steps |

**Settings** (main menu) changes sig figs (default 4), screen width and lines per page.

## What's on each menu (follows the formula sheet)

**Sec 0: Constants / units:** constants table; temperature C/K/F; pressure, energy,
volume, mass, length conversions (shown as factor-label); metric prefixes;
ppmv ↔ mole fraction; your own factor-label chain.

**Sec 1: Phases / particles**
- 1.1–1.5: intensive vs extensive; heating/cooling curve segments (with optional
  q = mcΔT + q = mΔH energy totals); phase at a given T,P for water, CO₂ or a custom
  diagram; every phase change along a constant-P or constant-T path; normal m.p./b.p.;
  b.p. at any pressure from vapor-pressure readings; volatility ranking; separation
  techniques; design a separation (cool a gas mixture / heat a liquid mixture);
  fractional-distillation trays and fractions; "which are liquid between T₁ and T₂".
- 1.6–1.10: ⟨Ek⟩ = ½m⟨v⟩² (any unknown, or from T and M); speed ratios at the same T;
  P = k_B NT/V (any unknown); P₁V₁/T₁ = P₂V₂/T₂ (any unknown, with sanity check);
  "what if N, T, V change" ratios; real-gas deviations; PEC written-justification
  template; particulate-model notes and common traps.

**Sec 2: Moles / composition**
- 2.1–2.3: classify matter (element/compound/mixture, molecular/ionic); name ↔ formula
  (ionic with Roman numerals and polyatomic ions, binary molecular with prefixes);
  molar mass with every atom counted; grams ↔ moles ↔ particles ↔ atoms (also liquid by
  volume + density); molarity forwards and backwards; mass %, mass % ↔ molarity;
  dilution.
- 2.4: PV = nRT (any unknown, moles or grams); M = mRT/PV, d = PM/RT, n/V = P/RT;
  identify an unknown gas (closest common gases + search by elements); STP molar
  volume and volume % of mixtures; pollutant conversions mg/m³ ↔ g/L ↔ mol/L ↔
  molecules/mL ↔ ppmv (plus temperature change); % increase; molecules of B per A;
  GWP per molecule.
- 2.5–2.7: p/n/e from symbol, A, charge (or the reverse); average atomic mass and
  2-isotope abundance; identify MS peaks (M⁺, lost fragments, fragment ion formula,
  M+1/M+2); isotope peak patterns for Cl, Br or custom (X₂ and more).
- 2.8: empirical formula from mass % or grams (with the ×2/×3/×4/×5 rule), combustion
  analysis from CO₂ and H₂O masses, molecular formula from M⁺, mass % from a formula.

**Sec 3: Light / Lewis / shape**
- 3.1–3.4: λ ↔ ν ↔ wavenumber ↔ E per photon ↔ kJ/mol (with region and color);
  photoelectric effect E_I = hν − E_k; all emission/absorption lines from energy levels;
  λ_max for an activation/bond energy; electron configurations for atoms and ions
  (noble-gas notation, unpaired e⁻, orbital boxes, Cr/Cu-type exceptions) with PES
  peaks; identify an element from a configuration or PES heights; rank radius / IE /
  electronegativity with reasons; rank atom/ion sizes (isoelectronic series);
  valence / bonding capacity; radical check; light & atom notes.
- 3.5–3.8: **Lewis structures from any formula** following the sheet's six steps
  (central atom, valence-electron count, single bonds, lone pairs, multiple bonds,
  check), with resonance forms, radicals, expanded and incomplete octets; VSEPR
  electron & molecular geometry and angles for every center; bond δ and dipole
  direction; polar/nonpolar with reason; IMFs present; geometry from domain counts;
  organic shortcuts; |δ| and μ (C·m and Debye) for any bond; IR band assignments;
  expected IR bands, IR activity/greenhouse and functional groups for a molecule.

**Sec 4: IMFs / ionic**
- 4.1–4.3: molecular vs ionic (Δχ, % ionic); IMFs in a substance; interactions and
  mixing/solubility between two substances; **rank boiling points** with the written
  justification; polymer and protein structure effects.
- 4.4–4.6: usual ion of an element; formula unit from two ions (criss-cross, reduce,
  parentheses) + name; unit-cell counting; Coulomb's law; rank ionic melting points;
  solubility rules with the dissolving equation; precipitation reactions (swap
  partners, net ionic and balanced molecular equation).

**Sec 5: Reactions / energy**
- 5.1–5.3: chemical vs physical change; six reaction assumptions & rate factors;
  energy diagrams (ΔH, Ea forward/reverse, exo/endo, ASCII sketch); multi-step
  diagrams; equilibrium extent from step probabilities; λ_max for photo-activated
  reactions; **balance any equation** (atom-balance algebra, fractions cleared,
  atom check).
- 5.4: grams A → grams B (A given as g, mol, liquid volume, gas PVT or solution);
  limiting reactant, theoretical yield, excess left, % yield; CO₂ emission factor,
  air-fuel ratio, rich/lean.
- 5.5: bond-energy table; **ΔH from a balanced equation automatically** (draws every
  Lewis structure, counts bonds, uses 799 for C=O in CO₂) or from bonds you enter;
  per gram / per amount of fuel; E_p of a molecule; q = mcΔT (any unknown); fuel
  heating water.

## Limits worth knowing

- Lewis structures: rings (benzene, cyclohexane) aren't supported. For organic
  molecules typed as molecular formulas (C₂H₆O), the app assumes a straight chain in
  formula order and says isomers are possible. Type a condensed formula to choose
  one. Formal charges are shown as an extra; the sheet's method decides the electron
  placement (for example, SO₄²⁻ is drawn with single bonds).
- Phase diagrams are modeled from the triple point, critical point and normal points
  (water and CO₂ built in), so crossings far from those points are approximate.
- Boiling-point ranking follows the sheet's method (size first if very different,
  then H-bond > dipole > dispersion, straight > branched). It is a reasoning aid,
  not measured data.
- Transition-metal ion charges must be typed in, or read from the anion in a formula.

## Tests (for developers)

```
python3 tests/run_all.py                 # every menu tool, scripted keypresses (CPython)
python3 tests/run_all.py path/to/micropython   # same under MicroPython with a 169 KB heap
python3 tests/heap_paths.py path/to/micropython # smallest heap each heavy screen needs
python3 tests/lewis_check.py             # Lewis/VSEPR/polarity for ~80 molecules
```
