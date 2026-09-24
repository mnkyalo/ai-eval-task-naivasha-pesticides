# Design notes

These notes cover how the task evolved, why each change was made, and what the evidence says about the current version. I publish the reasoning with the task because the reasoning is where the design skill lies.

## 1. Starting point: a real campaign and a real analytical problem

Abbasi & Mannaerts (2018) deployed silicone rubber sheets and Speedisks at sites in the Lake Naivasha basin in June–July 2016 to measure organochlorine pesticides. Converting a passive sampler's accumulated mass into a water concentration is routine work for environmental analytical chemists. It becomes a real modelling problem when two things are true, and both are true here:

- the samplers span a very wide range of equilibration times, and
- the water concentration did not stay constant during the deployment.

The task keeps the published campaign's geometry, polymer, PRC suite, sampler types and period. It then replaces every measurement with seeded synthetic data. That means the generating truth is known exactly, so the answer key is the truth itself rather than a published estimate.

## 2. Version history

**v2: nested windows only.** Four sheet sets, all retrieved together but deployed on different dates, so each set weights the declining concentration differently. The difficulty was supposed to come from recognising that the sets disagree in order of the compound's memory length (τ = Kpw·m/Rs). Probe solvers found it too easy. Once the brief says "the sets differ only in deployment length", fitting one shared decline across them is a natural move.

**Biofouling addendum, rejected before building.** I tested whether biofouling, as a second mechanism pulling the other way, could mask the ordering of the sets. The feasibility check found no version that was both realistic and hidden, so no data were generated for it. Stopping at this stage costs an afternoon. Stopping after a build costs a review round.

**v3: co-deployed Speedisks.** Adsorptive disks integrate linearly for every compound but carry no PRCs, so they have to be calibrated in situ against the sheets. They also supply the one piece of information the sheets cannot: the 30-day average of compounds that equilibrate on a sheet within a day, whose sheets remember nothing earlier than that. The disks exist for a real reason in the campaign design, which matters (see section 4).

**v5: the dissolved-organic-matter term.** An adsorptive disk also collects compound carried on dissolved organic matter, which a partitioning sheet does not see. I modelled the excess as `1 + phi·Kpw` per site. It is about 1.06–1.1× at HCB and 3.5–4.9× at mirex, which is the size reported for adsorption-based samplers. I added five hydrophobic legacy organochlorines (HCB, cis- and trans-chlordane, trans-nonachlor and mirex) so that every site has enough compounds the sheets are still integrating to separate G from phi. Without them the two parameters cannot be identified.

## 3. How the grading was sized

- **Tolerance:** within 20% relative on at least 54 of 57 concentrations, and within 10% on all three F values.
- **Honest routes:** the reference has a worst-case error of 5.7%. Five defensible alternatives have worst cases of 8.5–11.8% and pass 23–24 of 24 regenerated datasets.
- **Wrong routes:** the best wrong answer from the design record reaches 45/57, and every wrong answer passes 0 of 24 regenerated datasets.
- **The 54/57 bar** tolerates up to three outliers from a legitimately different but noisier method. The 20% band sits above every honest route and below every wrong one.
- **The F band** of 10% absorbs blank-handling slips, measured at up to 9%. The unit error of reporting sheet mass in grams instead of kilograms is 1000× and fails outright.

The validation harness in this repository adds independent checks. My own honest route, which reports the disk for the equilibrated compounds, lands at 57/57 with an 11.8% worst case. My own wrong answers land at 23/57, 32/57, 38/57, 9/57 and 0/57.

## 4. What the evidence says now, and what I would change

**The task's weak point is that the binding term shows up as a clean trend.** In my probe, all three fresh agents named the Kpw trend in the disk-to-sheet ratio as uptake bound to DOC or colloids within their first fitting pass. Two scored 57/57. The third had every mechanism right and failed only on a unit conversion. A trend that appears on the first diagnostic plot gets found rather than discovered, and tuning its size down would only shrink the error it produces. The lesson I took from this and other tasks: difficulty has to come from a mechanism that leaves no pattern in a standard first look at the data, not from one that is merely quiet.

**Blank correction does not discriminate.** The blanks are small compared with the accumulated masses. Skipping them gives a 14.4% worst-case concentration error and a 4.5% error on F, both inside the bands. The harness reports this as NOT CAUGHT rather than claiming the blanks are a trap. If blank handling is meant to count, the blanks need to be large enough to matter for the low-concentration compounds, and the grading should look closely at exactly those compounds.

**Where I would take a next version.** The disk calibration should depend on something the solver has to establish for a reason of its own. One candidate is a per-site DOC measurement that exists in the dataset because it is standard water-chemistry output, and that only reveals the binding term once the disk-to-sheet comparison has been set up correctly. Before building anything, I would run the plan-from-documentation test on it: a fresh model reads only the brief and the methods, with no data, and writes its analysis plan. If that plan already contains the correction, the mechanism is not hidden, and no amount of data design will hide it.

## Reference

Abbasi Y, Mannaerts CM (2018). Evaluating organochlorine pesticide residues in the aquatic environment of the Lake Naivasha River basin using passive sampling techniques. *Environmental Monitoring and Assessment* 190:349. https://doi.org/10.1007/s10661-018-6713-4
