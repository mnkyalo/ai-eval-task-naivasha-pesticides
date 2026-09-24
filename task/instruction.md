# Freely dissolved organochlorine pesticides, Lake Naivasha basin

Passive sampler data from a June-July 2016 campaign at three sites in the Lake Naivasha catchment, Kenya, is in `/app/data/`.

## Deployment

Three sites were instrumented with silicone rubber passive samplers: Upper Malewa (river), Middle Malewa (river) and Lake Naivasha (open water).

Each site received four sampler sets, designated A, B, C and D, of three replicate sheets each. The sets were deployed on different dates and all were retrieved together on a single visit.

| Set | Deployed     | Retrieved    |
|-----|--------------|--------------|
| D   | 20 June 2016 | 20 July 2016 |
| C   | 30 June 2016 | 20 July 2016 |
| B   | 10 July 2016 | 20 July 2016 |
| A   | 17 July 2016 | 20 July 2016 |

The staggered deployment was a logistical arrangement: boat and vehicle access to the three sites was available on four dates, and all sheets were collected together on the final visit.

Sheet dimensions are 55 x 90 x 0.5 mm. Silicone density is 1.15 g/cm3.

Every deployed sheet and every reference sheet was spiked before deployment with a suite of performance reference compounds (PRCs) at a nominal 120 ng per sheet per compound, for in-situ calibration. Reference sheets were carried to the field, kept sealed, and returned unexposed.

Sampling rates for the silicone sheets follow

    Rs = F * A / M^0.47

where Rs is in L/day, A is the sheet surface area in cm2 (both faces, edges excluded), M is molar mass in g/mol, and F is a site-specific coefficient absorbing local hydrodynamic conditions. F is not known in advance and differs between sites, the PRCs are there to determine it. It does not differ between sets at the same site.

Each site also received three Speedisk samplers (H2O-philic divinylbenzene extraction disks, 0.6 g sorbent), mounted beside the sheets. They were deployed with set D on 20 June 2016 and retrieved with everything else on 20 July 2016. Speedisks carry no PRCs. Sorbent loadings in this campaign were far below capacity.

Speedisk uptake is controlled by the water boundary layer, so Speedisk sampling rates follow the same molar-mass dependence,

    Rs_SD = G / M^0.47

where Rs_SD is in L/day and G is a site-specific coefficient. G is not known in advance and differs between sites.

## Files in /app/data/

| File                             | Contents                                                                                              |
|----------------------------------|-------------------------------------------------------------------------------------------------------|
| `prc_exposed_samplers.csv`       | PRC mass remaining on each exposed sheet (ng), by site and set                                        |
| `prc_reference_samplers.csv`     | PRC mass on the 3 unexposed reference sheets (ng)                                                     |
| `target_compound_masses.csv`     | Pesticide mass accumulated on each exposed sheet (ng), by site and set                                |
| `compound_properties.csv`        | Silicone-water partition coefficient (log Kpw, L/kg) and molar mass                                   |
| `procedural_blanks.csv`          | 3 sheets taken through the full handling, extraction and analysis chain, never spiked, never deployed |
| `speedisk_compound_masses.csv`   | Pesticide mass accumulated on each exposed Speedisk (ng), by site                                     |
| `speedisk_procedural_blanks.csv` | 3 Speedisks taken through the full handling, extraction and analysis chain, never deployed            |

All masses are ng per sampler (per sheet or per Speedisk). Partition coefficients are on a L/kg basis, so a sheet's uptake capacity depends on its mass, not its area.

## Task

Report the time-weighted average freely dissolved water concentration of each of the 19 target organochlorine pesticides at each of the 3 sites, over the full 20 June - 20 July 2016 campaign period.

Use all the sampler data available to you.

## Outputs

Write `/app/answer.csv` with exactly these columns:

    site,compound,cw_ng_per_L

- One row per site x compound: 57 rows plus the header.
- `site` and `compound` spelled exactly as they appear in the input files.
- `cw_ng_per_L` is the time-weighted average freely dissolved concentration in ng/L over the full 30-day campaign period, 20 June to 20 July 2016.
- Give at least 4 significant figures. Do not round to integers.

Write `/app/sampling_rates.csv` with exactly these columns:

    site,F

giving the fitted site coefficient F of the silicone sheets, one row per site, 3 rows plus the header. Use A in cm2 (both faces, edges excluded) and M in g/mol, so that F carries the units implied by Rs = F * A / M^0.47 with Rs in L/day.

You have 7200 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
