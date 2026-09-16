# Phase-4 label schema

Phase-4 errors are generated at UAV nodes and propagated to interferometric
edges. Arrays use MATLAB order; phase/coherence stacks are `[K,H,W]`, node
profiles are `[N,H]`, and graph connectivity is `edge_index[2,K]`.

## Observation fields

| Field | Shape | Unit | Meaning |
|---|---:|---:|---|
| `wrappedphase_node_noise_only` | `[K,H,W]` | rad | Complex node-SLC observation before Phase-4 errors |
| `wrappedphase_withnoise` | `[K,H,W]` | rad | Final observation after Phase-4 errors |
| `wrappedphase_multilook` | `[K,H,W]` | rad | Complex multilooked final observation |
| `valid_edge_mask` | `[K]` | boolean | Edges available to the model |
| `coregistration_valid_mask` | `[K,H,W]` | boolean | Pixels valid after node warping; inactive edges are zero |

Training losses must use both `valid_edge_mask` and
`coregistration_valid_mask`. A zero-filled unavailable edge is not a valid
zero-phase measurement.

## Node-level Phase-4 truth

Stored under `phase4_errors`:

| Field | Shape | Unit |
|---|---:|---:|
| `sync_phase_error_rad` | `[N,H]` | rad |
| `trajectory_phase_error_rad` | `[N,H]` | rad |
| `los_range_error_m` | `[N,H]` | m |
| `range_displacement_px` | `[N,H]` | pixel |
| `azimuth_displacement_px` | `[N,H]` | pixel |
| `node_parameters` | `[N]` struct | mixed, named units |

For edge `(u,v)`, the phase-error target is node `u` minus node `v`.
Relative displacement can be derived with the same incidence convention.

## Failure labels

Stored under `failure_labels` and duplicated in patch metadata:

- `initial_edge_mask`: random variable-K mask before explicit failure.
- `final_edge_mask`: final observation mask after UAV dropout.
- `edge_dropout_mask`: edges removed by an unavailable UAV.
- `uav_available_mask`: physical node availability.
- `sync_jump_node_mask`: nodes whose jump exceeds the configured threshold.
- `low_coherence_fraction`: fraction of active pixels below the threshold.
- `dropout_present`, `sync_anomaly_present`, `low_coherence_present`.
- `shortest_baseline_active`, `longest_baseline_active`.
- `compound_failure_present`: all configured failure conditions coexist.
- `failure_factor_count`: number of active factors among dropout, sync anomaly,
  low coherence, and longest-baseline observation.

The profile names end in `unvalidated`: thresholds and distributions are
simulation hypotheses until calibrated against flight data.
