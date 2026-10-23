# Power and Cell

**Programmable, self-dividing current cells in standard CMOS: how they divide, how they are programmed, how they exchange power, how far the idea reaches down to single atoms, and how it could grow into controllable energy conversion. A pre-registered model and simulation study, closed against published measurements.**

Research by **counfield** — <https://counfield.com/> · Türkçe: [README.tr.md](README.tr.md)

> **Status (October 2026).**
> - Everything in this repository is **model and simulation**. Nothing was built or measured by this project. The only measured data are
>   third-party: public cryogenic characterisation data of the process that was modelled, a published data set of silicon ESD current
>   filaments (Mamanee 2012) and results from the literature.
> - Every experiment keeps a `predictions_*.json` that was written **before** the run. Failed predictions are kept and reported. Some
>   reports start with a correction box; where a box conflicts with the original text, the box applies.
> - The per-experiment reports are in Turkish (`results/expN/summary_*.md`). This README is the English overview.
> - Experiments 27–30, 34 and 39–41 (the outlook line, Section 9) are **not yet part of this repository**; a few reports refer to them.

---

## Why "power and cell"

A **cell** here is a self-regulating channel that carries current. It grows with the current it carries, divides when it becomes too
wide, and its daughters divide again. **Power** is what the cells carry and exchange: when one cell dies, its neighbours take over its
current; when the drive changes, cells are born, split or disappear in a predictable way.

The work started from one question, whether a structure on an ordinary chip can grow and divide by itself under a steady electrical
drive, and grew into four:

1. What sets the size at which a cell divides, and can the answer be written as a law without fitted constants?
2. Can birth, number, survival and division be **programmed** with the electrical drive?
3. What is a cell worth as a unit that carries and recycles power, and how far down can programmability be pushed, to single atoms?
4. Which obstacles of controllable energy conversion could programmability remove? (Section 9, to be released.)

## Highlights

| Result | Number | Experiments | Evidence |
|---|---|---|---|
| A closed-form split law with no fitted constants, tested with predictions saved before the runs | predicted / observed drain-drop margin **1.001 ± 0.069** over 22 splits; 1.029 ± 0.071 (n = 20) after re-running with an independent load solver | 17, 20 | model, pre-registered |
| Birth site chosen by a +30 mV seed over 10 µm | **10 of 10** runs; the birth-site rule was right in **77 of 80** test runs | 18 | model |
| Cell number set by the height of a sudden drive step | 1, 2, 3–4 and 4 cells; a forced shrink splits all cells at once (**3 → 7**) | 19, 20 | model |
| Power handed to the neighbours when a cell is killed | **98.5–102.3 %** (recycling efficiency 0.985–1.023) | 23 | model |
| Floating-body kink in public 4 K characterisation data | present in **10 of 10** n-channel devices (ratio 2.8 against 1.15 at 300 K for one device) | 32, 33 | **measured (third party)** |
| Single-filament current limit | **0.889 A** predicted against 0.9 A measured | 31 | **measured (third party)** |
| Writing multi-atom patterns | global drive alone: 0 of 7 alternating and 0 of 7 random patterns; closed-loop write–read–correct: **21 of 21** | 35, 36 | model |
| Closure against independent measurements | 28 results compared: **16 agree, 6 mixed, 3 disagree**, 3 method-only | 29–32 | mixed |
| Parameter-free thrust law against a published ion-thruster test set | **40 of 40** operating points within 2σ, RMS **1.45 %** | 29 (outlook) | **published data** |
| Division at room temperature | **0 divisions in 21 runs** | 21 | model, **negative** |

---

## 1. Question and model

### 1.1 Two observations behind the question

- **(i) Stable, sized filaments.** Bistable (S-shaped) current–voltage characteristics are known to form current filaments in
  semiconductors (Schöll; Hüpper et al. 1993). A *global* constraint, such as a load resistor, selects a single filament whose width
  grows with current.
- **(ii) Replicating spots.** Reaction–diffusion systems with a *local* inhibitor show self-replicating spots (Pearson 1993).

Our question: can a standard CMOS layout provide the local inhibitor, so that a current filament in silicon behaves like a dividing
cell? If so, what sets the division size?

### 1.2 The device model

- **Device.** A wide n-channel transistor of a standard open 130 nm CMOS process, in a reduced model built from the foundry compact-model
  card (BSIM4, simulated with ngspice into tables of drain current, impact-ionisation current and threshold voltage) and extended with
  explicit, **swept** cryogenic assumptions: subthreshold slope, impact-ionisation multiplier, body-diode onset, frozen-well leakage and
  lateral body sheet resistance.
- **Local bistability.** Each segment of width $dy$ has a body potential $V_b$ with
  $C_b\,\dot V_b = \kappa\,\partial_y^2 V_b + F_{ii}\,M(V_D)\,j_d(V_b) - j_0 e^{(V_b - V_\text{on})/V_s} - g_\text{sub}V_b$.
  For a range of local drain voltage $V_D$ this has a low and a latched state: a floating-body latch that conducts in only part of the
  channel width, a **current filament**.
- **Inhibitor.** The drain is a resistive rail ($r_D$, Ω/µm). Each point of the rail is tied to a common node $V_c$ through a distributed
  ballast ($r_v$, Ω·µm). The drop $u = V_c - V_D$ obeys $u'' = u/\lambda^2 - r_D\,j_d(y)$ with $\lambda = \sqrt{r_v/r_D}$.
- **Reference set.** SS = 20 mV/dec, $F_{ii}$ = 10, gate overdrive −0.1 V, $V_\text{on}$ = 1.05 V, $V_s$ = 5 mV,
  $g_\text{sub}$ = $10^{-13}$ A/(V·µm); two further sets (U1, U2) test universality.

---

## 2. A filament that divides (Exp 12–13, 16–17, 20)

- A wide transistor biased just below threshold forms a current filament. If the drain is segmented and every segment is tied to the supply
  through its own ballast resistor, the filament **grows with current and divides in two, and its daughters divide again**.
- **The split law.** A filament splits when the drain voltage at its centre falls to $V_\text{down}$, the lowest voltage at which the
  conducting state can exist. For an isolated filament this gives a closed form,
  $w_s = 2\lambda\ln[\lambda/(\lambda-\lambda_c)]$. Its constants come from the device model; nothing is fitted.
- **Pre-registered test.** In the clean protocol (load $R_L$ = 5 kΩ) the centre-to-edge drain drop at 22 observed splits equals the
  predicted margin within a factor of **1.001 ± 0.069** (range 0.89–1.11). For two very different device-parameter sets it is
  0.98 ± 0.13 (9 splits). Across all protocols 63 splits were analysed. First-split widths of 14 clean runs follow the closed form with an
  RMS error of 9 %.
- **Independent re-check (Exp 20).** The explicit load line of the original solver was found to leave an unstable state after a sudden
  current excess. With an implicit load line the 40 Exp 17 runs were repeated:

  | quantity | original solver | implicit solver |
  |---|---|---|
  | mechanism ratio, clean | 1.001 ± 0.069 (n = 22) | 1.029 ± 0.071 (n = 20) |
  | mechanism ratio, stiffer load | 1.038 ± 0.062 (n = 22) | 1.003 ± 0.057 (n = 24) |
  | mechanism ratio, universality sets | 0.979 ± 0.133 (n = 9) | 0.973 ± 0.116 (n = 9) |
  | first-split widths, clean | n = 14, mean +1.6 %, RMS 9.1 % | n = 13, mean +3.2 %, RMS 8.6 % |

  13 of the 14 published first splits are reproduced. One run that split at 5.38 V (+20 %) no longer splits up to 5.6 V, and one
  "no division" statement for a very-low-leakage run was withdrawn (it splits at 175 µm, about twice the prediction).
- **One prediction failed.** At λ = 36 µm, 1 % above $\lambda_c$, the filament did not split.
- **Equation discovery (Exp 16).** A sparse-regression pipeline was first validated on known equations: it recovers a reaction–diffusion
  system exactly from clean data, fails at 1 % noise in the strong form, and recovers the impact-ionisation law in 90 of 90 curves.
  On the filament data it selects the width law $w = I/j_\text{on}$ ($j$ ≈ 2.55 ± 0.03 µA/µm, R² ≥ 0.9988). A fit with one free constant
  returns 1.01 × the theoretical value. Blind symbolic regression found no comparable compact law. Two checks failed and are reported:
  a finite-element sphere-capacitance convergence and a scaling test that the data could not decide.

## 3. Dimension decides division (Exp 11, 14, 15)

- **A window of dimensions.** In continuum models of chemically active droplets, division happens only inside a window of dimensions.
  Below the window (d ≈ 2–2.2) the droplet becomes a stripe; above it the droplet empties from the inside and becomes a liquid shell.
  The window moves up and widens as the internal consumption rate $\Delta\nu$ increases: d ≈ 2.25–2.5 at $\Delta\nu$ = −0.004,
  2.3–3.5 at −0.001, 2.5–4 at +0.004, and 6D droplets divide into two compact daughters for $\Delta\nu \ge +0.012$.
- **A genuine 6D run.** A Cahn–Hilliard active droplet was run on a 6D hypercube (16⁶ = 16.8 million points). An O(5)-symmetric
  reduction reproduces it to ≤ 2.5 % up to t = 40, and the high-dimensional division and shell results come from that validated reduction.
- **3D versus 2D.** A 3D active droplet divides into two equal daughters over the whole scanned drive range, and its division radius varies
  continuously with the drive. In 2D the same solver gives worms and labyrinths (binary division in 0 of 30 runs; sustained division
  2 of 30 against 0 of 28 for a supercapacitor-type model, Fisher p = 0.49, not significant).
- **Autocatalytic spots.** In a Gray–Scott scan, **5 of 42** (F, k) points keep replicating spots; 6 of 6 survive 5 % noise, and new
  spots are born by division (no births from the background).
- **A transistor channel is 2D and gated.** There the supercapacitor-type and active-droplet mechanisms do not give cell-like division;
  only autocatalytic systems remain candidates. In silicon the dividing geometry is one-dimensional: in 2D arrays the current spot
  elongates into a worm instead of dividing.
- **Phase-space holes (Exp 15).** Electron holes in a collisionless plasma merge in 1D and break up and disappear in 2D and 3D. No
  persistent division was found in any dimension; this is consistent with classical results (Morse and Nielson 1969).

## 4. Programming the cells (Exp 18–21)

| Lever | Result |
|---|---|
| **Where** a cell is born | The first cell is born at the single cell with the largest over-drive (19 of 20 calibration runs). The rule predicted the birth site in **77 of 80** test runs. A +30 mV seed over 10 µm made the first cell appear only there in **10 of 10** runs; the needed margin is about 2.5 times the random per-cell threshold spread ($\sigma$ = 12 mV per 0.5 µm). |
| **How many** are born | A sudden drive step of height $V_s$ sets the count: $V_s \le 1.0$ V gives 1 cell, 1.05–1.1 V exactly 2 (10 of 10 fields), 1.2 V gives 3–4, 1.4 V gives all 4 seeds (10 of 10). No cell was ever born outside a seed. Under a slow sweep one seed always wins (strongest input 19 of 20 on a ring device). |
| **Which survive** | Two knobs: the pulse height sets the births, the hold voltage sets the survivors (a carrying capacity). On-chip common-node capacitances ($R_L C \le 5$ ms) do not change the result. |
| **Damage and regeneration** | A slowly dying cell hands its current to the neighbours, which widen (3 → 2 cells at 3.3 V); no new cell is born (0 of 16 runs). A sudden death can spawn 2–3 new cells at the damage edges (3 → 4–5). Weakening the neighbour's centre by 150 mV makes it divide at once (3 → 2 → 3). At 5.0 V the widened neighbour regenerates the lost cell by itself, with a division-order ratio of 0.93. |
| **Moving a cell** | A moving gate pattern (+40 mV ahead, −40 mV behind) carries a cell to the device end in 5 steps (3.3 V) or 9 steps (5.0 V). |
| **Splitting at once** | A sudden forced shrink (releasing a gate pit, or a $V_{DD}$ drop of at least 0.2 V) splits the cells simultaneously: **3 → 7** at 5.0 V, −60 mV. |

**Limits found along the way.**
- The drive is a *single-winner* selector at birth, not a *k-winner* one: two equal seeds were born together in only 2 of 10 runs, and the
  winners of a step are not the k strongest seeds. A later-found predictor (the total sub-threshold weight) is better but did not meet its
  registered threshold in every new field.
- **At 300 K there was no division in 21 runs (Exp 21).** Bistability is weak (window 43–68 mV against about 770 mV at 4 K; on/off
  current ratio 5–12 against about $10^4$); the well does not freeze, so the whole 600 µm body acts as one node; and device mismatch eats
  the thin window. Four-kelvin division needs a wide bistable window and independent micro-cells.

## 5. Energy between cells (Exp 22–23)

- **Sub-bands.** Inversion-layer sub-bands exist in every transistor of the process (47–72 % of the electrons sit in the lowest
  sub-band at 300 K). Coulomb blockade would need ≤ 0.4–0.8 K; the filaments themselves are classical. A dissipative cell in this device
  class cannot reach the size a quantum dot needs (a quantum target at 4 K asks for a cell radius of 6–20 nm; the most optimistic
  dissipative cell is about 123 nm in a 2D layer).
- **Power exchange.** When a cell is killed, **98.5–102.3 %** of its power reappears in the surviving cells (recycling efficiency
  η = 0.985–1.023 over slow and sudden death at 3.3 V and 5.0 V; the worst case, slow death at 3.3 V, loses 1.5 %). The share each neighbour
  receives depends on proximity (about 73 % / 27 % for two neighbours), and a gate pre-bias changes it by only ±4 %.
- **Clamping and a frank caveat.** The cell voltage is clamped at $V_M$, so a cell does not concentrate energy. Recycling is nearly
  lossless, **but the architecture is not**: 74.6 % (3.3 V) and 83 % (5.0 V) of the source power is dissipated in the load resistor, and
  only 24 % and 16 % reaches the cells. A cell also stores almost no energy (body energy of order $10^{-14}$ J): it is a *flow* structure,
  not an energy store. Reducing that load loss is an open design problem (Section 9).

## 6. From cells to atoms (Exp 24–26, 35–36)

**6.1 The cell cannot shrink to atom size (negative).** Three floors were derived and tested: the frozen-body coupling length
(about 2.2 µm in this model), the depletion width (about 83 nm) and a 10-year stability width of at least 36 nm. Programmability
disappeared below the coupling length. "The programmable cell shrinks to an atom" is false for the electronic filament; it holds for an
ionic filament in an oxide (single-atom width in kinetic Monte Carlo), which is a different device class.

**6.2 Atoms trapped in a cell can be set to 0/1 (model).** A frozen cell body holds acceptor atoms that are either neutral or ionised;
each atom has its own switching voltage, its "address". At 1 K a 100 nm body has about 35 atoms with a margin of at least 8 kT/e on both
sides, and about 31 are set with 88 % success. At 4 K the method that works is a body of about 60 nm with addresses calibrated at 4 K:
about 3.3 atoms are programmed with ≈ 99 % accuracy. At 40 nm the addresses drift out of the window. The registered thermal predictions
(mean accuracy ≥ 0.9) failed: 0.82 at 4 K and 0.88 at 1 K.

**6.3 Fields programme logic (model, 4 K).** The cell field (0.4–1.8 V/µm) shifts a donor atom's hyperfine and Zeeman energies by kHz–MHz,
enough to select the atom by frequency. A NAND gate built from three natural acceptor atoms works in 27 of 28 triples and 9 of 10
islands when the thresholds are calibrated for every input state. The first method had a window error (3 of 28) that stays on record.

**6.4 A universal atom (model).** A central donor electron interacts with the nearby nuclear spins and can drive them with a complete gate
set. A single RF coil builds electron + 4 nuclei GHZ entanglement with fidelity **0.999**, 0.987 when about 2000 further nuclei are included,
and about 0.91 with an electron coherence of at least 10 ms. Published silicon GHZ fidelities are 92.5–96.2 %, so the model is optimistic
(one of the three disagreements in Section 7).

**6.5 Natural versus purified silicon.** In natural silicon only 1–3 % of the strongly coupled neighbours can be selected one by one
(spectral crowding), and storage collapses under measured conditions. A universal programmable atom that reaches **every** neighbour is
possible in purified ²⁸Si (at 50 ppm all ²⁹Si within 6 nm in 98 % of configurations) and not in natural silicon; the obstacle is the
nearest-neighbour ²⁹Si–²⁹Si pairs, which purification removes. Published experiments in purified silicon and diamond (cited in the
Exp 26 report) agree with the material condition.

**6.6 The control stack (Exp 35–36).**
- **A single global drive cannot write an arbitrary pattern:** monotone patterns 7 of 7, alternating 0 of 7, random 0 of 7.
- **Closed-loop write–read–correct** with a global correction reaches the right pattern in **21 of 21** tasks (selected 3–4-atom
  patterns, modal), marginal 0.992, full-pattern 0.973. Open-loop address calibration alone failed (2 of 14 hard patterns).
- **The price is background disturbance:** 6.25 % of non-target atoms flip against a registered limit of 5 %. The result is full control
  of a *selected target subset*, not isolation of the whole body.
- **Spatial resolution decides:** with a physical surface kernel and single-shot readout a 6×6 grid passes the registered full-control
  rule (13 of 14 modal, background 4.70 %), and an 8×8 grid with a finer step gives 14 of 14 (background 3.95 %).

## 7. Closure against measured data (Exp 31–33)

- **Current filaments in silicon (Mamanee 2012, 300 K).** At fixed filament number the current density per filament is constant
  (0.1 % and 5.2 %), the number law $N \approx I/(0.75\,j\,w_s)$ falls inside the measured ranges (0.84–0.90 of the measured mean),
  and the single-filament limit is **0.889 A against 0.9 A** measured. **The law fails in one respect:** the current density falls by
  22 % from N = 2 to N = 4 (the model had a single $j_\text{on}$).
- **Public 4 K characterisation of the modelled process.** The floating-body kink is present at 4 K in **10 of 10** devices (ratio 2.8
  against 1.15 at 300 K for a W7/L0.15 device). With the measured impact-ionisation current and transfer curve the cell constants become
  $V_M$ = 0.70 V and $j_\text{on}$ = 64 µA/µm (model constants, not a measurement of the filament). The body effect is 2–3 times smaller than
  the Exp 12 model assumed (measured in the reverse direction only).
- **Not measured:** the filament itself and its hysteresis. The widest public device is 7 µm, has no ballast, and was swept in one direction
  only.
- **Scoreboard (Exp 32).** All 28 results of Exp 0–27 were compared with an independent published measurement: **16 agree** (including
  those corrected by the measurement), **6 are mixed or partial**, **3 disagree** (the sub-band solver of Exp 22, the simulated GHZ of Exp 25
  and the simulated coherence of Exp 26) and 3 are method experiments.

## 8. What did not work (kept on purpose)

- The early attempt at a room-temperature **quantum** cell (Exp 0–10) is currently unsuccessful: scale, coupling and coherence do not
  reach the required regime.
- No division at 300 K (Exp 21); no persistent division of phase-space holes (Exp 15); no division of a transistor-channel droplet in 2D
  (Exp 11).
- The split law fails in its filament-number dependence (Exp 31) and in the λ = 36 µm prediction (Exp 17); one first split and one
  "no division" claim changed after the independent re-check (Exp 20).
- Most of the source power is dissipated outside the cells in the present architecture (Exp 23).
- The programmable cell cannot shrink to atom size in silicon (Exp 24); the registered thermal accuracy of atom programming was missed.
- Three model results disagree with measurements (Exp 22, 25, 26); open-loop addressing failed and the global-correction loop exceeds the
  registered background limit (Exp 35).
- Several analysis steps in the equation-discovery study failed or could not be tested (Exp 16).

## 9. Outlook: controllable "fuel cells" (to be released)

> **Status.** Experiments 27–30, 34 and 39–41 form a separate line that treats programmable cells as **controllable energy cells**
> ("fuel cells" in the sense of units that convert and release a controlled amount of power, not electrochemical cells). It is
> theory-level; only parts of it are closed against published data. The reports will be released with their predictions and data in a later
> update. The numbers below are quoted from those reports.

**9.1 The idea.** Cells supply positive energy; atoms (ions) carry the momentum. A cell is a self-regulating current channel with a
characteristic current of order 0.2 mA; the number of cells sets the delivered power in steps; a dead cell hands its power to its
neighbours.

**9.2 What is already established.**
- **Routes without reaction mass are closed.** A closed vehicle cannot exceed a thrust per power of 1/c = 3.34 nN/W; vacuum-fluctuation
  routes give at most 2 × 10⁻²¹ N from the best measured source; negative-energy metric routes would need 39–71 orders of magnitude more
  energy than any accessible vacuum energy (Exp 27, negative result).
- **The atom route carries far more momentum per unit of energy.** An exhausted atom carries 2/v_e of momentum per unit of spent energy
  and a photon 1/c, a ratio 2c/v_e: $10^4$ at 60 km/s and $6 \times 10^4$ at 10 km/s (Exp 28).
- **A parameter-free thrust law, $F/P = 2\eta/v_e$ (with $F = \gamma J_b\sqrt{2mV/q}$ for the beam), was tested against published data**
  (Exp 29–30). It predicts **40 of 40** bench operating points of a published ion thruster inside 2σ with no free parameter (mean
  deviation +1.0 %, RMS **1.45 %**; 1.18 % for the 33 points above 1000 V), reproduces a flight-measured point (92 mN, 2.3 kW) with
  η = 0.61, and matches the directly measured thrust of a liquid-metal field-emission thruster to 1.6–4.3 %.
- **Atom-level links are measured in published experiments:** on-demand single ions (88.5 %, detection-limited), the $\sqrt{m}$ scaling
  of flight time (1.97σ) and $I = e f$ for an artificial atom (0.27 ppm).
- **Model envelope.** In the model each cell yields 7–26 µN (1–5 kV), the cell number sets thrust in steps, and thrust changes by at most
  0.8 % when one cell dies (Exp 28). An energy-optimal operating envelope was derived from all earlier data (Exp 34), including
  sub-predictions that failed (a registered range that was off by $10^3$ because of an arithmetic slip is corrected in the report).
- **Exploratory studies** of mission-scale scaling and of extreme scale descent (Exp 39–41) are speculative and will be released with
  that label.

**9.3 Obstacles that programmability could remove.**

| Obstacle for controllable energy conversion | Programmable lever | Evidence here | Status |
|---|---|---|---|
| Power cannot be throttled in fine, reversible steps | cell number set by the step height; hold voltage sets survivors; forced shrink splits all cells (3 → 7) | Exp 19–20 | model; no division at 300 K (Exp 21) |
| One local failure takes down the whole unit | a killed cell hands 98.5–102.3 % of its power to neighbours; slow death spawns no new cell (0 of 16); deliberate damage is repaired by division (3 → 2 → 3) | Exp 18, 23 | model |
| Current crowding and uncontrolled filamentation | a global constraint selects one filament; distributed ballast makes it divide at a predictable width, $w = I/j_\text{on}$ | Exp 12–17 | model; ESD filaments measured (Exp 31) |
| The place where activity starts is random | a +30 mV seed fixes the birth site (10 of 10; rule 77 of 80) | Exp 18 | model |
| Energy is lost when a unit is switched off | recycling between cells is nearly lossless (η 0.985–1.023) and the cell voltage is clamped at $V_M$; **but 75–83 % of the source power is lost in the load resistor** | Exp 23 | model; **open** |
| Programmability needs feedback | write–read–correct loop reaches the target pattern in 21 of 21 tasks; a global drive alone cannot (0 of 7) | Exp 35–36 | model |
| A unit cannot be made small enough | the electronic cell has floors (≈ 36 nm stability, ≈ 2.2 µm body coupling); programmable acceptor atoms inside the body are the smaller unit | Exp 24 | model; floors are limits, not overcome |
| Material purity limits control | universal control of all neighbours needs purified material | Exp 26 | model + cited measurements |
| Room-temperature operation | the electronic filament does not divide at 300 K; an ionic filament in an oxide is the candidate class | Exp 21, 24 | negative; simulation only |
| Validation gap | 28 results compared with independent measurements; the thrust link is closed (40 of 40) | Exp 29–32 | partly closed |

**9.4 What is still open (and next).**
1. The only unmeasured link of the chain is the **cell itself at cryogenic temperature**: current steps of a seeded, ballasted device.
2. A seeded birth-site test on a real device (Exp 18 is a model; seeded devices are reported to fix the birth site in measurements of
   silicon filaments).
3. A way to cut the load loss of Section 5, a device class that divides at room temperature (ionic filaments are the candidate), and a
   closed-loop readout fast enough for Section 6.6.
4. Release of Experiments 27–30, 34 and 39–41 with their data, so that every number in this section can be reproduced.

## 10. Method: pre-registration and honesty

- **Predictions before runs.** Each experiment writes a `predictions_*.json` before the simulation; the file is write-once, and later
  scripts must not overwrite it. Failed predictions are kept and reported next to the passing ones.
- **Disclosure.** When a prediction was not blind (for example when the sources were read before the criteria were written), the report
  says so.
- **Corrections stay visible.** Audits and re-checks are added as boxes at the top of the affected reports; the text below a box is kept for
  the record.
- **Closure.** Each model result is compared with an independent measurement when one exists, and marked "no measurement" or "method" when
  none does. Nothing is fitted to the comparison.
- **Negative results** are part of the result set (Section 8), not a footnote.

## 11. Repository layout and reproduction

| Path | Contents |
|---|---|
| `src/qcell/` | models: finite elements, Poisson–Schrödinger, phase field, active matter, device tables |
| `experiments/` | one folder per experiment (`expN_*`), runnable scripts |
| `results/` | `expN/`: Turkish summary reports, `predictions_*.json` (written before the run), data, figures |
| `tests/` | 74 pytest tests |
| `data/` | records and licences of the third-party device models and measurement data |

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
PYTHONPATH=src pytest -q                              # 74 tests
```

Equation-discovery experiments need `requirements-optional.txt` in a separate environment. Some experiments need `ngspice` (device
tables; the tables are included) or FEniCSx/dolfinx (optional cross-checks); each experiment's report lists what it used. Large raw
arrays are not stored; they are regenerated by re-running the scripts.

## 12. Limitations

- **Model, not measurement.** The 4 K device behaviour is an extrapolation from 300 K compact-model tables plus swept assumptions. Not
  included: self-heating, trap kinetics, the parasitic bipolar transistor and source-side resistance.
- **Drive and load.** The split law assumes a quasi-static sweep. A stiffer load ($R_L$ = 2 kΩ) delays splits by about 4 %, or 14 % with low
  mismatch; very low body leakage delays or suppresses them.
- **Scatter and threshold.** At individual splits the criterion scatters by ±7 % (1σ, clean protocol). Close to threshold the predicted
  $w_s$ is extremely sensitive to the constants.
- **Untested pieces.** The periodic-array formula and the filament-count guide line are not validated.
- **2D results.** Most voltage steps in the 2D runs did not reach the steady-state tolerance; the 2D conclusion (no division) is weaker
  evidence than the 1D results. Active-droplet results in $d \ge 4$ come from the O(d−1)-symmetric reduction, validated in 3D and against a
  6D hypercube during growth only.
- **Atom-level results are models.** They assume idealised control fields, and the closure table shows where they overestimate (Section 7).
- **Not quantum.** The cell results are classical, dissipative pattern formation.

## 13. Selected references

- A. Hüpper, E. Schöll et al., Phys. Rev. B 47, 15515 (1993); E. Schöll and D. Drasdo (1990): S-shaped current–voltage characteristics
  and current filaments.
- J. E. Pearson, *Complex patterns in a simple system*, Science 261, 189 (1993).
- D. Zwicker et al., *Growth and division of active droplets provides a model for protocells*, Nat. Phys. 13, 408 (2017).
- W. Mamanee, *Investigation of current filament behavior in smart power ESD devices*, TU Wien thesis (2012).
- B. Morse and C. Nielson, Phys. Rev. Lett. 23, 1087 (1969): phase-space holes.
- T. D. Ladd et al., Phys. Rev. B 71, 014401 (2005): nuclear-spin coherence in silicon.
- Published atom-level experiments cited in the Exp 26 report: Pla et al. (2014), Asaad et al. (2020), Mądzik et al. (2022),
  Thorvaldson et al. (2025), Edlbauer et al. (2025), Bradley et al. (2019).

## 14. Licence, third-party data and citation

MIT licence (see `LICENSE`). The device-model files and the public measurement records are third-party material under their own licences;
see `NOTICE.md`. Published papers are cited, not redistributed. To cite this work, see `CITATION.cff`.
