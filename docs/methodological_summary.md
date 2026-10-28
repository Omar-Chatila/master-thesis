# Movement segmentation and behavioral annotation for state-dependent random walks

## 1. Purpose and scope

This document describes the methods implemented in **hmmcma**, with emphasis on **BCPA, hidden Markov models (HMM), the BCPA–HMM hybrid, and criteria-based segmentation**. It explains what each method does, how it works, and what information passes through its modules. The intended use is as a basis for a thesis or paper, rather than as a programming tutorial. Pseudocode is language agnostic.

The motivating application is to group observed displacements with comparable movement properties, use these groups to estimate **state-dependent random-walk transition kernels**, and apply those kernels to interpolation and utilization-distribution (UD) calculations. The package implements the segmentation and annotation part of that workflow. Kernel estimation, path conditioning, and UD estimation remain downstream tasks.

Descriptions of the current implementation are based on the source files in this checkout. General theory is identified separately and linked to references. A fixed-budget Oystercatcher implementation ablation is reported in Section 12; it is distinct from general method descriptions and is not held-out validation.

**Naming:** the implemented hybrid class is `BCPAHMM`, and the configuration method is `bcpahmm`. This is the method referred to as HMMBCPA in the motivating request. The package's `BCPA` class uses PELT with an RBF segment cost; it does not reproduce the original likelihood-based BCPA algorithm of Gurarie and colleagues.

## 2. Motivation: grouping displacements rather than merely cutting tracks

For locations $x_i$ observed at times $t_i$, define the incoming displacement and elapsed time as

$$
d_i=x_i-x_{i-1},\qquad \Delta t_i=t_i-t_{i-1}.
$$

A single distribution fitted to every displacement can mix several movement regimes: small local movements, persistent travel, repeated turns, or transitions between these patterns. Segmentation and annotation provide a way to identify observations that may reasonably inform different movement distributions.

Three objects must be distinguished:

- A **segment** is a consecutive part of a trajectory.
- A **state** is a recurring movement class. Multiple separated segments can belong to the same state and contribute to the same kernel.
- A **behavioral name**, such as resting or foraging, is an interpretation of a state. Statistical separation alone does not establish that interpretation.

The relevant modeling goal is therefore both temporal and distributional: identify where movement changes and determine which intervals share sufficiently similar dynamics to be pooled. HMM, BCPA, the hybrid, and criteria rules address these two questions differently.

| Method | How boundaries arise | How recurring states arise | Main organizing principle |
| --- | --- | --- | --- |
| HMM | Runs of the decoded hidden state | Learned emission distributions and state-transition probabilities | Explain feature sequences through latent movement states. |
| Package BCPA | PELT detects changes in feature distributions | Cluster summaries of detected segments | Detect episodes first, then group comparable episodes. |
| BCPA–HMM hybrid | PELT supplies boundaries | Aggregate HMM state evidence within each detected segment | Retain BCPA episodes and assign each one an HMM state. |
| Criteria-based segmentation | A global partition satisfies labeled temporal/geometric rules | Segments satisfying the same labeled rule share a state | Define acceptable movement properties explicitly. |

## 3. Shared input and feature preparation

### 3.1 Input → preparation module → output

```text
INPUT
    trajectories with individual/trajectory IDs and observation times
    projected x and y coordinates in metres
    selected movement features (terrain IDs excluded)

PREPARATION MODULE
    retain row identities and original geometry
    preserve invalid-coordinate/time fixes as barriers
    identify blocks with approximately constant sampling intervals
    calculate features inside blocks and retain contiguous finite rows
    retain sequence identity and a mapping back to output rows

OUTPUT
    ordered feature sequences for fitting
    original observations with derived feature columns
    row mapping for attaching states and segment boundaries
```

The modern HMM, BCPA, and hybrid share this preparation. The criteria engine has its own preparation path, described in Section 7.

The package expects projected columns named `utm_x` and `utm_y`; it does not select or construct their projection. They must describe a coherent metric coordinate system within each trajectory. The original geometry is preserved. Features use planar differences of these metric coordinates; the output geometry is not reprojected.

### 3.2 Available features and their interpretation

| Feature | Meaning | Relevance for displacement grouping | Implementation detail |
| --- | --- | --- | --- |
| `distance` | Incoming step length | Movement magnitude | Metres; alternative to speed. |
| `speed` | Step length / elapsed seconds | Movement rate | Metres per second; not path speed between fixes. |
| `direction` | Clockwise heading from north | Absolute orientation | Displayed in degrees, modelled as sine/cosine on the unit circle. |
| `angular_difference` | Unsigned change of heading | Turning magnitude | Degrees; an alias of turn_angle, not an independent dimension. |
| `turn_angle` | Unsigned change of heading | Turning magnitude | Radians, bounded by pi. |
| `persistence_velocity` | Current velocity projected onto the preceding heading | Continuation and reversal | Metres per second; positive for continuation. |
| `turning_velocity` | Signed lateral component relative to the preceding heading | Sideways movement | Metres per second; positive for clockwise turning. |
| `acceleration` | Backward change of step speed / current elapsed seconds | Changes in movement rate | Optional additional dimension, m/s². |

With signed clockwise heading change $\theta_i\in[-\pi,\pi)$, the persistence and lateral
components are

$$
p_i=\frac{\|d_i\|}{\Delta t_i}\cos\theta_i,
\qquad
q_i=\frac{\|d_i\|}{\Delta t_i}\sin\theta_i.
$$

Persistence alone is the compact default. The optional pair represents speed
and signed turning without separately adding distance, speed or angle again. `speed` plus unsigned `turn_angle` is another compact
configuration. Known duplicate representations are rejected. Terrain class IDs
are rejected in the modern and legacy feature interfaces; no environmental data
are fetched or inferred.

Headings of zero-length steps are undefined. `min_direction_distance` can impose
an accuracy-based minimum length; its default zero excludes only zero steps.
Undefined first-step speeds, turns without two valid headings, and acceleration
without two valid speeds remain NaN. A zero current displacement can have zero
projected velocity if the preceding heading is defined. Missing values are never
interpreted as observed stillness or straight movement.

### 3.3 Sampling, gaps, missing values, and scaling

The modern pipeline partitions intervals around the first cadence of each block,
allowing `sampling_rtol=0.1` relative jitter plus `sampling_atol=1` second.
Anchoring prevents gradual drift from joining substantially different rates.
Rate changes below a factor of five also create boundaries. Sustained slower
rates are retained. An isolated interval above `max_gap_factor=5` times the
median, unsupported by either neighbouring interval, is treated as a gap.
No interpolation is performed.

Invalid coordinates/times and duplicate timestamps are barriers before feature
calculation. Each sampling block uses its immediately preceding fix only as an
anchor for its first displacement; turns and acceleration do not borrow history
from the preceding block. Only contiguous finite rows of the selected feature
representation enter modelling sequences. Persistence therefore needs two within-block displacements. Rows without model assignments remain
`state=-1`, with their original locations and timestamps retained where the
MovingPandas input/output representation supports them.

HMM and hybrid globally scale linear features. The library supports rate-group
fits (`sampling_model_scope="by_rate"`) and one shared fit across all prepared
sequences (`"shared"`). The current experiment uses only the main cadence and
a shared model, as detailed in Section 12. `sampling_interval_seconds` anchors
the acceptance band globally so clock jitter does not resplit a selected cadence. Sine/cosine heading coordinates retain equal unit-circle
scales. PELT scales each sequence separately, then globally scales segment
descriptors for clustering. These scopes differ. Tiny rate groups without enough
observations or transitions remain unmodelled and have an explicit status.

```text
PREPARE_DATA(trajectories, selected_features, cadence_tolerances):
    preserve original row identities and invalid-fix barriers
    for each valid contiguous portion of a trajectory:
        identify blocks matching the fixed target cadence, or constant-rate blocks if no target is set
        calculate incoming steps, valid headings and temporal derivatives per block
        retain contiguous finite feature rows, with block and sequence identities
    encode absolute directions by sine/cosine
    return original points, finite sequences, cadence records and row mapping
```

**Source:** [preprocessing.py](../preprocessing.py): `prepare_trajectory_collection`, `feature_matrix`, `normalise_features`, feature calculation and gap helpers.

## 4. Method 1 — Hidden Markov model

### 4.1 What an HMM is

An HMM describes a sequence of observed feature vectors through an unobserved sequence of discrete states. Each state has a distribution of observable features, and transitions between states follow a first-order Markov chain: the next state depends on the current state. Conditional on the state sequence, the standard HMM treats emissions at different times as independent. See [Rabiner's HMM tutorial](https://web.mit.edu/6.435/www/Rabiner89.pdf).

Its main components are the initial-state probabilities $\pi$, transition probabilities $A$, and state-specific emission distributions $f_s$:

$$
P(S_1=s)=\pi_s,\qquad
P(S_i=b\mid S_{i-1}=a)=A_{ab},\qquad
Y_i\mid S_i=s\sim f_s.
$$

The joint probability for a state and feature sequence is

$$
p(s_{1:n},y_{1:n})=
\pi_{s_1}f_{s_1}(y_1)
\prod_{i=2}^{n}A_{s_{i-1},s_i}f_{s_i}(y_i).
$$

State persistence arises through self-transition probabilities. For a time-homogeneous HMM, the implied residence time in a state is geometrically distributed in observation steps; with self-transition probability $A_{ss}<1$, its mean is $1/(1-A_{ss})$. This is an observation-step assumption, not a continuous-time duration model. [Rabiner, 1989](https://web.mit.edu/6.435/www/Rabiner89.pdf).

### 4.2 How fitting and state inference work

Three related computations serve different purposes:

1. **Forward computation:** sum over possible hidden-state paths to evaluate the observation likelihood.
2. **Baum–Welch / expectation–maximization:** alternate between expected state occupancies/transitions and parameter updates to improve that likelihood.
3. **Viterbi decoding:** find the single most probable complete state sequence under the fitted model. Forward–backward smoothing can instead provide marginal probabilities for each state at each observation.

These are standard HMM operations; marginally most probable states and the most probable complete path need not agree. [hmmlearn tutorial](https://hmmlearn.readthedocs.io/en/stable/tutorial.html).

### 4.3 What this package implements

The implementation uses a **Gaussian HMM with diagonal state covariance matrices**. Each state has a mean and a variance for every selected feature. Within a state, the model does not estimate covariance between feature dimensions.

In shared mode one model is fitted across all selected prepared sequences. The library also retains separate rate-group fits for unfiltered inputs; the current benchmark selects only the main cadence and uses shared mode. Sequence lengths prevent the fit from treating the end of one sequence and the start of the next as an observed transition. Multiple initializations are tried, using seeds 0–9 by default, and the model with the highest training log-likelihood is selected. Defaults include three states, at most 500 fitting iterations, likelihood tolerance $10^{-3}$, and a covariance-floor setting of $10^{-2}$.

The method decodes a state for each prepared observation. Labels are then ordered by the mean speed of the observations assigned to each decoded state; if speed was not selected, the first feature determines the order. This relabeling assists interpretation but does not change the fitted model.

All modern data-driven states have neutral names `state_0`, `state_1`, etc. Ordering occurs within the shared model, or separately within each rate group when using that library option. The main-cadence benchmark has one state basis across all retained blocks; equal ranks from separate rate-group models do not establish biological equivalence.

### 4.4 Input → HMM module → output

```text
INPUT
    prepared feature sequences
    chosen number of states and fitting settings

HMM MODULE
    optionally standardize features across all sequences
    fit several Gaussian HMM initializations
    retain highest-likelihood model
    decode a state sequence
    relabel states by speed or the first selected feature

OUTPUT
    state and cluster for each modeled observation
    sequence_id and contiguous state-run identity in segment_id
    feature columns, state mapping, descriptive names, fit metadata
```

```text
HMM_ANNOTATE(trajectories, features, K):
    prepared ← PREPARE_DATA(trajectories, features, gap_factor)
    Z ← optionally standardize all prepared feature sequences together
    candidates ← fit diagonal-Gaussian HMMs with K states in the selected shared/rate-group scope
                 using multiple initializations and sequence lengths
    model ← candidate with highest observation log-likelihood
    states ← Viterbi_decode(model, Z, respecting sequence boundaries)
    mapping ← order decoded states by mean speed, else first feature
    attach mapped states to original rows
    sequence_id ← prepared-sequence identity
    segment_id ← contiguous final state-run identity
    change_point ← within-sequence decoded state switch
    return common annotation result
```

**Output distinction:** `sequence_id` identifies a finite prepared sequence;
`state_run_id` and HMM `segment_id` identify contiguous decoded-state runs.
Runs restart at all sequence boundaries. `change_point` and the change-point
dictionary record within-sequence state switches, without treating rate/gap
boundaries as observed behavioural transitions.

### 4.5 Assumptions relevant to random walks

The Gaussian distributions describe selected **features**, not a spatial transition kernel automatically usable for interpolation. Nonnegative lengths and unsigned turns are Gaussian approximations with bounded support not enforced by the emissions. Absolute headings use a sine/cosine embedding, so 359° and 1° are adjacent. This is not a von-Mises model. Terrain is excluded.

The separate `MovementHMM` implements gamma step lengths and von Mises signed turning angles, with optional zero inflation and marginalization of undefined turns. These emissions are not options of the Gaussian fitting helper. See the [model comparison and API](movement-hmm-comparison.md), the [oystercatcher benchmark](benchmarks/oystercatcher_movement_hmm/report.md), and [Michelot, Langrock, and Patterson, 2016](https://theomichelot.github.io/files/2016_MEE_moveHMM.pdf).

An observation-step HMM also does not automatically adjust state-transition probabilities for unequal elapsed times. Separate cadence models avoid pooling their observation-step transition matrices. This still does not supply a continuous-time transition law or justify applying a matrix at a different time step.

**Sources:** [state_annotation.py](../state_annotation.py): `HMM.annotate`; [models.py](../models.py): `fit_gaussian_hmm`, `converged_cleanly`. Each selected fit records its seed, iteration count, finite likelihood, stop reason, likelihood history and unvisited states. The convergence flag describes that selected fit, requiring a near-monotone history and a final nonnegative increment below tolerance; an iteration-limit stop alone is not convergence.

## 5. Method 2 — BCPA: change detection followed by segment clustering

### 5.1 Original BCPA and the meaning of the name

The original behavioral change-point analysis treats movement-derived time series as observations of a continuous stochastic process and seeks structural changes in its parameters. This provides a statistical treatment of temporal correlation and irregular observation intervals. [Gurarie, Andrews, and Laidre, 2009](https://doi.org/10.1111/j.1461-0248.2009.01293.x).

The associated BCPA software searches for breaks in analysis windows, compares changes in mean, variance, and correlation structure using a modified BIC, and aggregates the window-wise evidence. This is useful conceptual background, but it differs from the implementation here. [Original BCPA package documentation](https://search.r-project.org/CRAN/refmans/bcpa/html/bcpa-package.html).

For this repository, the precise method description is **PELT-based feature change-point detection followed by clustering of segment summaries**.

### 5.2 Detecting temporal changes with PELT

PELT, short for Pruned Exact Linear Time, optimizes a penalized partitioning objective. If a candidate partition has $m$ segments, the objective has the form

$$
\sum_{j=1}^{m}C(Y_{\tau_{j-1}:\tau_j})+\beta m.
$$

Here, $C$ measures within-segment heterogeneity and $\beta$ penalizes additional segments. Charging per internal change rather than per segment differs only by a constant for a nonempty sequence. Dynamic programming compares possible preceding boundaries, while pruning removes candidates that cannot improve the optimum. Its computational advantage depends on the cost and data; the name is not an unconditional linear-time guarantee. [Killick, Fearnhead, and Eckley, 2012](https://arxiv.org/abs/1101.1438).

The package uses an **RBF kernel cost**, with similarity

$$
k(y,z)=\exp(-\gamma\|y-z\|^2).
$$

An idealized kernel-dispersion cost for a segment $I$ is

$$
C(I)=\sum_{i\in I}k(y_i,y_i)
-\frac{1}{|I|}\sum_{i,j\in I}k(y_i,y_j).
$$

In plain language, a segment is inexpensive when its feature observations are mutually similar. Splitting can lower this heterogeneity enough to compensate for a new-segment penalty. RBF embedding can detect distributional changes beyond a change in the ordinary arithmetic mean; it does not explicitly fit movement autocorrelation. Bandwidth selection is delegated to the library's median heuristic. [ruptures RBF cost documentation](https://centre-borelli.github.io/ruptures-docs/user-guide/costs/costrbf/).

**The RBF similarity kernel is a change-detection device. It is a different object from the downstream random-walk transition kernel.**

Detection runs separately on each prepared sequence after standardizing its features. The default penalty is 5 and the default minimum segment size is three observations. A sequence shorter than twice that minimum is retained as one segment without searching for internal changes. The wrapper explicitly sets and stores candidate-index spacing `jump=1` by default; callers can choose a larger grid to reduce computation. See [PELT implementation reference](https://centre-borelli.github.io/ruptures-docs/code-reference/detection/pelt-reference/).

### 5.3 Grouping detected segments into states

Detected segments initially describe temporal episodes, without recurring state identities. For each segment, the package summarizes every selected **unscaled** feature by its mean, standard deviation, and median. It concatenates these into a segment descriptor:

$$
q_j=[\operatorname{mean},\operatorname{std},
\operatorname{median}]_{\text{selected features}}.
$$

All segment descriptors from all trajectories are pooled, standardized, and clustered using K-means. K-means represents each group by a centroid and seeks to reduce the sum of squared distances from descriptors to their assigned centroids. It alternates assigning descriptors to nearby centroids and updating those centroids; multiple initializations reduce dependence on a single starting arrangement. [scikit-learn K-means documentation](https://scikit-learn.org/stable/modules/clustering.html#k-means).

The assigned cluster becomes the state of every point in that segment. The number of clusters is capped by the number of available segments. Defaults are three clusters and ten K-means initializations with a fixed random seed.

No sums or duration features are included by default, so repeating the same movement distribution for longer does not change the typical-movement descriptor. Headings are summarised through their sine/cosine coordinates rather than raw degree values. Each segment contributes one descriptor, rather than being weighted explicitly by its observation count or duration. Cluster labels are arbitrary identifiers and are not speed-ordered as in the HMM.

### 5.4 Input → BCPA module → output

```text
INPUT
    prepared feature sequences
    change penalty, minimum segment size, desired cluster count

BCPA MODULE
    standardize each sequence separately
    detect boundaries with PELT and RBF segment cost
    summarize each detected segment using original feature values
    standardize pooled segment descriptors
    cluster descriptors and propagate segment clusters to points

OUTPUT
    detected temporal segments and change-point times
    one recurring cluster/state label per segment
    common point annotation and selected movement features
```

```text
BCPA_ANNOTATE(trajectories, features, penalty, min_size, K):
    prepared ← PREPARE_DATA(trajectories, features, gap_factor)
    segments ← empty list
    for each prepared sequence Y:
        if length(Y) < 2 × min_size:
            boundaries ← none
        else:
            boundaries ← PELT(standardize(Y), RBF_cost, penalty, min_size)
        for each interval between boundaries:
            descriptor ← concatenate mean, std, median of Y in interval
            retain interval, descriptor, and row identities
    Q ← standardize all retained descriptors together
    labels ← K_MEANS(Q, clusters = min(K, number_of_segments))
    assign each interval's cluster to all its rows
    mark the first observation after each internal detected boundary
    return common annotation result
```

Two adjacent BCPA segments can have the same final cluster. A detected change point therefore need not be a final state transition. For kernel construction, recurring clusters allow pooling across episodes, while segment IDs preserve the temporal structure.

**Source:** [BCPA.py](../BCPA.py): `BCPA.change_points`, `segment_features`, `BCPA.annotate`.

## 6. Method 3 — BCPA–HMM hybrid

### 6.1 What the hybrid combines

The hybrid separates episode detection from recurring state assignment. BCPA supplies temporal boundaries; a Gaussian HMM supplies state evidence. One final state is assigned to each BCPA segment.

The HMM is fitted to the prepared **point-level feature sequences**, using the same fitting process as the standalone HMM. It is not trained on segment summaries, and it is not refitted after segment smoothing. The hybrid uses only BCPA's boundary detection; it does **not** use the standalone BCPA's K-means segment clustering.

The result is a staged smoothing procedure, rather than a jointly optimized HMM with fixed boundaries or a hidden semi-Markov model.

### 6.2 Segment assignment methods

For segment $I_j$, the default `posterior_sum` selects

$$
\hat s_j=\arg\max_s\sum_{i\in I_j}
P(S_i=s\mid Y_{1:n},\hat\theta).
$$

The probabilities are HMM-smoothed marginal state probabilities over the containing prepared sequence. Their sum measures expected occupancy in that segment. Choosing the largest sum selects the constant segment label with the greatest expected number of matching point states under those marginals. It is **not** the posterior probability that the entire segment remained in that state, and it is not constrained Viterbi decoding.

The alternative `majority_vote` selects the most frequent hard Viterbi state in the segment. Both options use an argmax; exact ties select the first raw model-state index before the speed-based relabeling.

### 6.3 Input → hybrid module → output

```text
INPUT
    prepared feature sequences
    BCPA boundary settings and HMM state/fitting settings
    segment assignment choice: posterior sum or majority vote

HYBRID MODULE
    detect BCPA boundaries separately in each sequence
    fit a Gaussian HMM on the full prepared feature sequences
    obtain point-level HMM paths and posterior state probabilities
    aggregate HMM evidence within each BCPA segment
    assign a constant final state throughout each segment

OUTPUT
    BCPA segment IDs and detected change points
    HMM-derived state for every segment and its points
    segment-to-state mapping, state ordering, and fit metadata
```

```text
HYBRID_ANNOTATE(trajectories, features, BCPA_settings, HMM_settings):
    prepared ← PREPARE_DATA(trajectories, features, gap_factor)
    segments ← detect PELT/RBF boundaries in each prepared sequence
    model ← fit multi-start Gaussian HMM on globally standardized features
    paths ← Viterbi_decode(model, each sequence)
    posteriors ← forward_backward_state_probabilities(model, each sequence)
    for each BCPA segment I:
        if assignment is POSTERIOR_SUM:
            raw_state ← argmax_s sum(posteriors[i, s] for i in I)
        else:
            raw_state ← most frequent paths[i] for i in I
        label ← speed-ordered mapping of raw_state
        assign label to every observation in I
    retain BCPA boundaries even when neighboring labels are equal
    return common annotation result
```

### 6.4 Consequences for state-dependent kernels

The hybrid groups point-level movement evidence into coherent episodes before supplying observations to state-specific kernels. Short HMM state excursions can disappear when the surrounding segment favors another state. Conversely, a BCPA episode containing a genuine internal state transition still receives only one final state.

The HMM's fitted transition matrix describes the unsmoothed latent model, whereas the returned segment labels describe the smoothed annotation. They are not interchangeable state-transition descriptions. If random-walk simulation uses transitions between final annotated states, those transitions need their own treatment.

**Source:** [hybrid_annotation.py](../hybrid_annotation.py): `BCPAHMM.annotate`, `_bcpa_segments`, `_segment_state`. The segment aggregation is a property of this implementation; no claim of a separately established hybrid algorithm is implied by the class name.

## 7. Method 4 — Criteria-based segmentation

### 7.1 What criteria-based segmentation is

Criteria-based segmentation defines a predicate on a candidate interval: does that interval satisfy a specified movement property? Rules associate predicates with labels, and a partitioning algorithm seeks a small number of valid labeled segments.

The geometric framework distinguishes conditions preserved when an interval is shortened, such as a bounded attribute range, from conditions preserved when it is extended, such as sufficient duration. Their combination motivates global partitioning rather than a simple greedy longest-prefix rule. See [Alewijnse et al., 2014](https://ad-research.cs.uni-freiburg.de/proceedings/SIGSPATIAL%202014/paper_full_alewijnse.pdf). This reference supplies theoretical context; the statements below about available rules and boundary behavior come from the local implementation.

Unlike the learned methods, criteria encode the analyst's definition of acceptable movement. A rule can be interpreted directly, but its label is only as defensible as its thresholds and underlying movement assumptions.

### 7.2 Available criteria

| Criterion | Meaning of the candidate-interval condition | Modeling use |
| --- | --- | --- |
| `minimum_hops(h)` | Requires at least $h$ consecutive transitions between fixes, subject to endpoint conventions. | Ensure an episode has enough movement observations. |
| `maximum_hops(h)` | Bounds the difference between endpoint observation indices. | Limit episode size in observation steps. |
| `minimum_duration(T)` | Requires endpoint elapsed time **strictly greater than** $T$ in the current implementation. | Exclude episodes too brief for the intended interpretation. |
| `speed(r)` | Bounds the range of interval speeds: $\max(v)-\min(v)\leq r$. | Identify intervals with relatively homogeneous speed. |
| `speed_with_outliers(r,q)` | An outlier-tolerant speed-range implementation with a supplied outlier count. | Reduce sensitivity to exceptional speed observations. |
| `heading(a)` | Bounds the smallest circular arc containing interval headings by $a$ degrees. | Identify movement with stable absolute orientation. |
| `heading_with_outliers(a,q)` | An outlier-tolerant circular heading-range condition. | Allow exceptional headings within an otherwise stable episode. |
| `turning_angle(a)` | Bounds the circular range of relative headings/turns by $a$ degrees. | Identify consistency in turning behavior. |
| `turning_angle_with_outliers(a,q)` | An outlier-tolerant circular turning-range condition. | Allow exceptional turns. |
| `disk(r)` | Tests spatial confinement using an approximate enclosing-disk radius at most $r$ metres. | Identify location-confined episodes. |
| `all_of(...)` | All constituent predicates must hold for the same interval. | Combine duration and confinement, for example. |
| `any_of(...)` | At least one constituent predicate must hold for the interval. | Admit alternative movement patterns under one label. |
| `negate(...)` | Complement the predicate over candidate intervals. | Express the failure of another interval condition. |

**Speed range is not a speed ceiling.** Constant speeds of 1 m/s and 20 m/s both have range zero. The speed criterion alone can therefore group either as homogeneous. Similarly, a narrow turning-angle range means consistent turning, not necessarily small turns or straight movement. Negating a homogeneous-speed condition does not turn it into a high-speed threshold.

The criteria engine calculates forward speed and heading from consecutive projected fixes. Relative headings use circular angle differences, whereas the shared data-driven `turn_angle` feature is unsigned. Threshold meanings must therefore be compared at the predicate level, rather than inferred from similar feature names.

Disk testing maintains extreme points in five directions and calculates an enclosing circle from those representatives. It is an approximation to confinement of the full point set, rather than a guaranteed exact minimum enclosing disk of all observations.

The outlier variants express the intention to tolerate a fixed number of exceptional values without removing the corresponding output rows. Their detailed interval and exclusion conventions differ from the ordinary variants; the speed-outlier loop does not enumerate every possible allocation of exclusions between the low and high ends. For publication-level claims, describe these as the implemented outlier-tolerant predicates rather than asserting exact equivalence to an ideal trimmed-range test.

### 7.3 Combining criteria into named rules

A rule consists of a label, a predicate or Boolean combination, and a segment-length preference. A conceptual rule for a confined episode is

```text
LABEL: location_confined
CONDITION:
    duration exceeds a chosen threshold
    AND approximate enclosing-disk radius does not exceed a chosen radius
PREFERENCE: prefer a longer interval among equivalent solutions
```

This defines a movement pattern. Calling it resting requires an additional biological justification: an animal may remain spatially confined while actively moving.

`LONGEST` and `SHORTEST` guide choices among equal-score candidates by preferring earlier or later starts within the admissible ranges. They do not replace the primary objective with globally longest or globally shortest segments. When multiple rules offer the same best score, earlier rules can win because the solver only replaces a stored solution on a strict improvement. Rule order is therefore relevant to interpretation.

### 7.4 How the partition is found

For every possible segment endpoint, the core represents the sets of starts that satisfy each rule. Boolean combinations operate on these admissible-start sets. A dynamic program chooses successive valid intervals, charging the same cost per segment; equivalently, it seeks the fewest segments over the solver's admissible domain. Efficient range-maximum queries select the best preceding partition and apply start-index tie preferences.

An idealized recurrence is

$$
D(j)=1+\min_{i,\ell:\ C_\ell(i,j)=\mathrm{true}}D(i),
$$

where $C_\ell$ is the predicate associated with label $\ell$. The source expresses the equivalent objective as maximizing a score that decreases by 0.5 per segment.

The current core solves endpoint positions through the **penultimate retained fix**, and the output adapter extends the final segment's label to the last fix. Internal boundary fixes are assigned to the segment beginning there. Consequently, the output should not be described as a proof that every final labeled point interval, including the terminal edge, satisfies an ideal inclusive-endpoint predicate.

If the solver cannot construct a valid partition, it returns a fallback unclassified segmentation. It does not generally label isolated local failures and then continue fitting every remaining interval independently.

### 7.5 Input → criteria module → output

```text
INPUT
    trajectory IDs, times, and metric coordinates
    ordered labeled rules with thresholds and tie preferences

CRITERIA MODULE
    sort valid fixes by time within each trajectory
    use the first fix at duplicate timestamps for computation
    construct admissible-start sets for each endpoint and rule
    combine predicates and optimize a valid labeled partition
    map segment labels back to original rows, including duplicates

OUTPUT
    segment_id, semantic segment_name, state, cluster, change_point
    label-to-state metadata
    unmodeled or fallback statuses where applicable
```

```text
CRITERIA_ANNOTATE(trajectories, ordered_rules):
    for each trajectory:
        valid ← finite coordinates and parsable times, sorted by time
        unique ← first valid fix at each timestamp
        if fewer than three unique valid fixes:
            keep original rows without an assignment
            continue
        for each rule:
            admissible_starts ← evaluate rule over candidate intervals
        partition ← dynamic_program_minimum_segments(
                        admissible_starts,
                        ordered_rules,
                        longest_or_shortest_tie_preferences,
                        terminal_endpoint = penultimate_fix)
        if no valid partition:
            assign fallback unclassified label
        else:
            assign selected rule and segment identity to each interval
            extend final segment label to final fix
        copy labels to duplicate-time rows
    map rule labels to numeric states and derive boundary flags
    restore original row order and return common annotation result
```

### 7.6 Differences from shared preprocessing

The criteria path does **not** apply the HMM/BCPA long-gap splitting rule, feature standardization, or shared feature extraction. It computes geometric predicates directly from coordinates and times. It requires at least three unique valid fixes for a trajectory to be segmented. Duplicate-time fixes receive the label computed from the first valid fix at that time.

State 0 means `unclassified`; supplied rule labels map to states 1, 2, and so on in rule order. Rows without assignments have state `-1`, which has a different meaning from unclassified. Segment numbering restarts for each trajectory, so use trajectory identity together with `segment_id` as a key.

The `observation_variance` setting is stored on the core's trajectory points. The exposed temporal/geometric predicates and minimum-segment solver do not use it to fit a measurement-error model. The compiled source also contains Brownian-bridge likelihood and variance-segmentation machinery, but the public criteria workflow does not invoke those fitting routines and does not produce a UD.

**Sources:** [criteria.py](../criteria.py); [criterion construction and row handling](../cpp/segmentation/csv_segmentation.cpp); [range and disk criteria](../cpp/segmentation/inw_monotone_crit.cpp); [criterion classes and endpoint settings](../cpp/segmentation/inw_monotone_crit.h); [minimum duration and hop count](../cpp/segmentation/outw_monotone_crit.cpp); [geometric feature definitions](../cpp/segmentation/trajectoryfunction.cpp); [partition solver](../cpp/segmentation/brownian_criterion.cpp); [tie handling](../cpp/segmentation/maxtree.cpp); [disk approximation](../cpp/segmentation/dynamic_smallest_encl_disk.cpp).

## 8. Common outputs and methodological diagnostics

All modern annotators return `AnnotationResult`, containing an annotated trajectory collection, numeric labels, a change-point dictionary, and metadata. The core fields have method-dependent meanings:

| Output | HMM | BCPA | Hybrid | Criteria |
| --- | --- | --- | --- | --- |
| `state` | Decoded HMM state | Segment cluster | HMM-derived segment state | Rule label encoded numerically |
| `cluster` | Copy of state | Copy of state | Copy of state | Copy of state |
| `segment_id` | Contiguous state run | Detected PELT segment | Detected PELT segment | Rule-solver segment, local to trajectory |
| `change_point` | Within-sequence state switch | Internal detected boundary | Internal detected boundary | Change in assigned segment ID |
| `segment_name` | Absent | Absent | Absent | Semantic rule label |

`evaluate()` describes coverage, state frequencies, transitions, boundary counts, segment lengths, and durations. It uses state runs as the HMM's primary segment basis and segment IDs for the other methods. It also reports both families of summaries. Length means observation count; duration uses first-to-last observation time. State proportions are observation proportions, not automatically time budgets.

State runs restart at prepared-sequence boundaries and unmodelled rows. PELT
detection segments remain separate even when their final states coincide;
`state_run_id` provides the additional final-label grouping. The hybrid retains
raw ordered labels, raw model indices, PELT segment IDs and aggregation method.
Both aggregation rules break ties by the smallest raw HMM state index. Segment
posterior shares describe expected occupancy, not a calibrated probability of
remaining in one state throughout the segment.

`state_feature_summary()` and `compare_state_feature_summaries()` describe feature count, mean, median, and standard deviation within states. These help explain what the states represent, but they are neither kernel estimators nor biological-validation scores. Plotting functions display states and segment boundaries without changing the annotation.

`SegmentationBuilder` and `SegmentationConfig` configure and instantiate the four methods. They do not introduce another statistical procedure. The legacy `HMMStateAnnotator` accepts precomputed feature columns, including caller-defined features, and returns annotated points, trajectory tuples, a typical interval in minutes, and state mappings. Its point states are zero-based and tuple states one-based; the modern annotation result uses a different contract.

The fitted HMM objects, scaling and raw/ordered transition matrices are retained under `metadata["sampling_models"]`. Ordered point posteriors are exported in `hmm_posterior_*` columns. Empirical adjacent final-label transition counts/matrices are stored separately, excluding sequence interruptions. Raw HMM transitions do not describe transitions between final PELT episodes or define a continuous-time process.

### Public interfaces and supporting functions

| Interface | Input → responsibility → output |
| --- | --- |
| `HMM`, `BCPA`, `BCPAHMM`, `CriterionAnnotator` | Collection and method settings → selected annotation algorithm → `AnnotationResult`. |
| `CriteriaSegmentation` | Alias for `CriterionAnnotator`; the same criteria method. |
| `Feature` | Selected feature names → validated movement-descriptor choices → feature configuration. |
| `CriterionConfig`, `CriterionType` | Thresholds and child predicates → declarative rule condition → criterion configuration. |
| `SegmentationRule`, `SegmentPreference` | Label, criterion, and tie preference → named interval class → rule supplied to the criteria solver. |
| `annotate_trajectory_csv`, `CsvAnnotationOptions` | CSV and rules → low-level criteria segmentation → output CSV path, with a `segmentation` column containing segment number and label. |
| `HMMStateAnnotator`, `ColumnConfig` | Precomputed-feature table and column mapping → legacy HMM workflow → annotated table and historical trajectory/interval/mapping outputs. |
| `SegmentationBuilder`, `SegmentationConfig`, `SegmentationMethod` | Method name and parameters → reusable configuration → fresh annotator when created. |
| `BaseAnnotator` | Common method interface → preparation, result handling, summaries, and diagnostics → consistent annotation workflow. |
| `state_feature_summary`, `compare_state_feature_summaries` | Annotated results → descriptive aggregation by state → feature-summary tables. |
| `plot_states`, `plot_segments` | Annotated result/collection → state or boundary visualization → saved figure. |

### Main parameter choices

| Setting | Methods | Meaning and default |
| --- | --- | --- |
| Feature set | HMM, BCPA, hybrid | Time-normalized persistence velocity by default; lateral velocity and acceleration are optional. |
| Number of states/clusters | HMM, BCPA, hybrid | Three by default; BCPA uses `num_clusters`, the others `num_states`. Selected by the caller rather than estimated automatically. |
| Change penalty | BCPA, hybrid | `penalty=5`; discourages extra segments. It is not a statistical significance level or the classical BCPA modified-BIC setting. |
| Minimum segment size | BCPA, hybrid | `min_segment_size=3` observations; a minimum count rather than elapsed duration. |
| Gap threshold multiplier | HMM, BCPA, hybrid | `max_gap_factor=5`; guards isolated long intervals, alongside explicit constant-cadence blocks. |
| Global feature scaling | HMM, hybrid | `scale=True`; BCPA additionally has its separate sequence and descriptor scaling. |
| HMM fitting budget | HMM, hybrid | Up to 500 iterations per initialization, seeds 0–9 by default. |
| Segment state assignment | Hybrid | `posterior_sum` by default; `majority_vote` is available. |
| Ordered predicates and preferences | Criteria | Supplied by the analyst; `LONGEST` is the default tie preference. |
| Observation variance | Criteria | 25 by default; retained in the core data but not used by the exposed geometric/temporal partition objective. |

**Sources:** [annotation.py](../annotation.py), [configuration.py](../configuration.py), [state_annotation.py](../state_annotation.py), [plot.py](../plot.py).

## 9. Connecting the annotation to state-dependent kernels

### 9.1 What segmentation contributes

Let $\hat s_i$ be the state assigned to observed displacement $d_i$. The annotated displacement sample for state $s$ is

$$
\mathcal D_s=\{(d_i,\Delta t_i,\text{context}_i):\hat s_i=s\}.
$$

BCPA creates these groups by clustering episode summaries; HMM creates them through decoded emission states; the hybrid pools BCPA episodes with the same HMM-derived label; criteria pool intervals matching the same rule. Retaining segment identity allows recurrence to be pooled without discarding temporal boundaries.

Before forming these samples, define **which displacement a point label annotates**. The shared data-driven features are largely attached to the arriving observation, while criteria use forward interval quantities and their own boundary conventions. A kernel-fitting adapter must specify incoming-edge or outgoing-edge assignment, transitions crossing a boundary, artificial initial values, invalid rows, and gaps. Point states cannot be copied blindly onto every pair of consecutive output rows.

### 9.2 What still needs to be estimated

A spatial transition kernel describes the distribution of the next location under a state and chosen duration:

$$
K_s(x'\mid x,\Delta t,h,E),
$$

where $h$ represents relevant history, such as previous heading, and $E$ any environmental context. A simple translation-invariant model could estimate a state-specific distribution of displacement vectors. A model with directional persistence would additionally condition on orientation or rotate displacements relative to the previous heading. These are downstream modeling choices.

An HMM emission distribution is a distribution of the selected feature vector. It does not automatically specify a valid distribution of geographic displacement vectors: lengths and angles need consistent units, support, orientation, and time-scale treatment. Likewise, PELT's RBF kernel measures feature similarity rather than movement probability.

State-specific fitting must determine how much to pool across individuals and recurring segments, how to accommodate small state samples, and how to treat uncertain assignments. The package's hard labels alone do not quantify that uncertainty. Using the same movement observations to define states and then fit their kernels also makes the partition an estimated component of the movement model, rather than independent biological ground truth.

### 9.3 Input → segmentation package → kernel module → output

```text
observed timestamped trajectories
    → selected segmentation/annotation method in hmmcma
    → states, segment identities, boundaries, movement descriptors
    → explicit point-label-to-displacement alignment
    → state-grouped displacement samples at known elapsed times
    → downstream state-specific movement-distribution estimation
    → spatial transition kernels plus a state-evolution model
```

```text
BUILD_STATE_KERNELS(annotation, original_observations, model_choices):
    edges ← construct consecutive observed displacements within valid sequences
    discard or separately model artificial initial values and gap-crossing edges
    states ← align point/segment labels to edges using a declared convention
    for each recurring state:
        samples ← corresponding displacement vectors, durations, and context
        kernel[state] ← estimate chosen movement law at the intended time scale
        normalize destination probabilities or continuous density
    specify persistence/state-transition assumptions for modeled paths
    return state-dependent kernels and state-evolution specification
```

This pseudocode expresses the intended connection; `BUILD_STATE_KERNELS` is not an implemented function in this repository.

## 10. Random-walk interpolation and utilization distributions

Interpolation asks for plausible paths during observation gaps **conditioned on the observations**. A freely simulated walk from the starting fix generally does not satisfy the next observed fix. A downstream bridge or other conditioning procedure must combine the transition model with endpoint information and, where appropriate, observation error.

State assignments inferred at observed times do not identify the exact state history within a long gap. The interpolation model must decide whether a state remains fixed, changes at specified boundaries, or follows a probabilistic state process. In particular, neither the HMM's observation-step transition matrix nor the hybrid's segment labels automatically defines continuous-time dynamics.

A finite-period observation-conditioned UD can be expressed conceptually as

$$
u(x)=\frac{1}{T}\int_0^T p(X_t=x\mid\mathcal O)\,dt,
$$

interpreting $p$ as a density or spatial probability mass as appropriate, and $\mathcal O$ as the observations. State-dependent kernels affect this occupancy by changing where the process can travel, how dispersed it becomes, and how long different movement dynamics apply.

```text
state-dependent kernels + state-evolution assumptions + observation model
    → endpoint-conditioned movement paths or location distributions
    → occupancy weighted by represented time
    → normalized spatial probability mass or density
    → utilization distribution for the stated period
```

With irregular observation times, counting annotated fixes does not estimate time spent in a state or place without further weighting assumptions. Distinguish grid probability mass from density per unit area, and distinguish the UD over an observation period from a stationary distribution of a random walk. These formulations describe the intended downstream framework, not implemented estimators in hmmcma.

## 11. Suggested thesis or paper description

> Trajectories were annotated to identify recurring movement regimes whose observed displacements could inform state-dependent transition kernels. Four segmentation strategies were considered: a Gaussian hidden Markov model; PELT-based distributional change detection followed by clustering of segment summaries; a hybrid assigning HMM-derived states to PELT segments; and a criteria-based partition defined by labeled temporal and geometric predicates. Segments describe consecutive episodes, whereas recurring state labels allow displacement information to be pooled across episodes. The annotation methods differ in whether states are learned from point-level feature distributions, from episode-level summaries, or from explicit movement constraints.
>
> The segmentation layer supplies movement-state assignments and temporal structure. Subsequent kernel estimation requires an explicit alignment of point labels with displacement intervals and a treatment of elapsed time, observation gaps, and directional history. The resulting transition kernels can support observation-conditioned random-walk interpolation and time-weighted occupancy estimation. Biological activity names require justification separately from the statistical or geometric definitions of movement states.

For a completed study, replace general alternatives with the actual feature choices, thresholds, state counts, interval conventions, fitting settings, and downstream movement model used. Call the implemented BCPA method **PELT/RBF change-point segmentation with segment clustering** when distinguishing it from classical BCPA.

## 12. Oystercatcher comparison at the main observation interval

The current experiment first selects **one dominant observation cadence** and
benchmarks only its contiguous sections. Main cadence is determined from times
alone: find the most populated positive-interval band under a 10% relative plus
one-second tolerance, then use that band's median as the target. This avoids
mistaking the low-jitter exact mode of 16 seconds for the largest cadence group.
The current target is **615 seconds (10.25 minutes)**, with accepted incoming
intervals from 552.5 to 677.5 seconds. `--main-dt-seconds` can instead specify a
known nominal interval, for example 600 seconds; the tolerances are configurable.

From 11,625 original fixes in 15 source episodes, this retains **633 contiguous
blocks with 3,354 fixes**. Every block has at least three fixes. Its first two
fixes provide step/heading history and are excluded from every variant's scoring
cohort, leaving **2,088 identical scored observations**. Excluded intervals,
invalid fixes and changes away from the target rate are barriers. No new edge is
created by dropping an interval, and no interpolation is used. Source-row and
trajectory identities, retained/discarded runs and the exact scoring mask are
exported for audit.

HMM and hybrid each fit **one shared model across all retained main-cadence
sequences**, using sequence lengths to suppress transitions across section
boundaries. Blocks are neither fitted individually nor assigned independent
state semantics. There is one common state basis for each feature configuration.
PELT still detects boundaries within each contiguous sequence, then clusters its
segment descriptors globally. Both hybrid label rules reuse identical HMM fits
and PELT boundaries.

The frozen pre-change implementation is **refitted on this same main-cadence
selection**, including its required context fixes; it is not scored by reusing a
model trained on all original rates. Comparisons fix three states, seeds 0–2,
150 iterations and PELT penalty 8. PELT minimum size is 3; its old jump is 5 and
its new jump is 1. No parameter search is used. Reference label permutation
maximizes macro state IoU on the same scored observations: these are descriptive
in-sample comparisons, not independent validation. Other sampling rates are
outside the benchmark denominator; unmodelled scored fixes would still count as
errors. All 48 completed runs cover 100% of the common scoring cohort.

The table uses the previous **application** configuration as its baseline. The
new `speed + turn_angle` feature set and persistence-only library default are
reported separately. The API-default baseline and all isolated changes also
appear in the complete report.

| Method | Previous application IoU | New speed + turn IoU | Change | New persistence IoU | Change |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gaussian HMM | 0.2359 | 0.4165 | +0.1806 | 0.3947 | +0.1588 |
| PELT/RBF + clustering | 0.2193 | 0.2738 | +0.0545 | 0.1982 | -0.0210 |
| Hybrid posterior sum | 0.2316 | 0.3735 | +0.1419 | 0.3718 | +0.1403 |
| Hybrid Viterbi majority | 0.2377 | 0.3937 | +0.1560 | 0.3733 | +0.1356 |

Removing distance alone under the old preparation does not improve the HMM here
(0.2359 to 0.2305). The revised finite-feature preparation with the compact
speed/turn set raises it to 0.4165. For compact PELT, dropping length-dependent
sums changes IoU from 0.2546 to 0.2620; the revised complete compact pipeline
reaches 0.2738. The isolated jump change has no observed score effect on these
mostly short blocks. Linear and sine/cosine absolute-heading HMM variants both
score 0.4102; circular geometry remains the reason to avoid linear degrees.
Thus the previous mixed-rate experiment's conclusions about individual changes
are not carried over to this new cohort.

For the new compact hybrid, Viterbi majority reaches 0.3937 versus 0.3735 for
posterior sums. With persistence alone, their values are 0.3733 and 0.3718.
The exact pre-permutation label differences and per-class confusion tables are
reported alongside the method results. The fitted raw HMM transition matrix
remains distinct from empirical transitions between final smoothed labels.

Extra dimensions remain optional. Lateral velocity plus persistence gives HMM
IoU 0.2370; adding acceleration to that pair raises it to 0.3703, while
persistence alone reaches 0.3947. PELT responds differently: its persistence-only
score decreases against the application baseline, whereas speed/turn improves.
These findings do not support a universal feature expansion or a universal
method ranking.

The scored reference contains **1,239 inactive/maintenance, 828 foraging and
only 21 flight/travel fixes**. Selection uses no behavioural labels, but it
changes the represented behaviour distribution; the flight-class result has
particularly limited support. State count, cadence tolerance and biological
interpretation still require independent validation. Terrain is absent, so no
terrain or land/water effect is fabricated.

The [current report](benchmarks/oystercatcher_main_dt/report.md),
[per-run records](benchmarks/oystercatcher_main_dt/runs.json),
[score table](benchmarks/oystercatcher_main_dt/scores.csv),
[manifest](benchmarks/oystercatcher_main_dt/manifest.json),
[selected blocks](benchmarks/oystercatcher_main_dt/selected_blocks.csv) and
[scoring cohort](benchmarks/oystercatcher_main_dt/scoring_cohort.csv) preserve the
selection, parameters, per-class results and input/source hashes. The
[runner](../benchmarks/oystercatcher_ablation.py), application entry point and
parameter search now all select the main cadence and use shared HMM fits.
Identical deterministic fits are cached, so recorded stage times are not
independent runtime benchmarks. The [earlier mixed-rate report](benchmarks/oystercatcher_ablation/report.md)
is retained as historical evidence; its numbers use a different cohort and are
not the baseline for the current table.

## 13. References and implementation map

### Scientific and algorithmic references

1. **Gurarie, E., Andrews, R. D., and Laidre, K. L. (2009).** A novel method for identifying behavioural changes in animal movement data. *Ecology Letters*, 12, 395–408. [DOI](https://doi.org/10.1111/j.1461-0248.2009.01293.x). Original BCPA motivation and continuous stochastic-process treatment; not the algorithm implemented by this package's BCPA class.
2. **Killick, R., Fearnhead, P., and Eckley, I. A. (2012).** Optimal detection of changepoints with a linear computational cost. *Journal of the American Statistical Association*, 107, 1590–1598. [Author manuscript](https://arxiv.org/abs/1101.1438). PELT optimization and pruning.
3. **Rabiner, L. R. (1989).** A tutorial on hidden Markov models and selected applications in speech recognition. *Proceedings of the IEEE*, 77, 257–286. [Paper](https://web.mit.edu/6.435/www/Rabiner89.pdf), [DOI](https://doi.org/10.1109/5.18626). HMM structure, inference, and parameter estimation.
4. **Michelot, T., Langrock, R., and Patterson, T. A. (2016).** moveHMM: an R package for the statistical modelling of animal movement data using hidden Markov models. *Methods in Ecology and Evolution*, 7, 1308–1315. [Author copy](https://theomichelot.github.io/files/2016_MEE_moveHMM.pdf), [DOI](https://doi.org/10.1111/2041-210X.12578). Movement-specific HMM emissions and interpretation beyond the Gaussian implementation here.
5. **Alewijnse, S. P. A., Buchin, K., Buchin, M., Kölzsch, A., Kruckenberg, H., and Westenberg, M. A. (2014).** A framework for trajectory segmentation by stable criteria. *ACM SIGSPATIAL*. [Paper](https://ad-research.cs.uni-freiburg.de/proceedings/SIGSPATIAL%202014/paper_full_alewijnse.pdf), [DOI](https://doi.org/10.1145/2666310.2666415). Theoretical context for admissible intervals, stable criteria, and partition optimization; does not establish that every facility in the paper is exposed by this package.

### Primary software documentation

- [hmmlearn tutorial](https://hmmlearn.readthedocs.io/en/stable/tutorial.html) and [GaussianHMM API](https://hmmlearn.readthedocs.io/en/stable/api.html#hmmlearn.hmm.GaussianHMM): Gaussian emissions, sequence fitting, decoding, and posterior calculations.
- [ruptures RBF cost](https://centre-borelli.github.io/ruptures-docs/user-guide/costs/costrbf/), [PELT guide](https://centre-borelli.github.io/ruptures-docs/user-guide/detection/pelt/), and [PELT implementation](https://centre-borelli.github.io/ruptures-docs/code-reference/detection/pelt-reference/): segment cost and admissible-boundary settings.
- [Original BCPA package documentation](https://search.r-project.org/CRAN/refmans/bcpa/html/bcpa-package.html): window-sweep likelihood and modified-BIC analysis, to distinguish classical BCPA from this implementation.

### Local implementation map

| Module | Methodologically relevant responsibilities |
| --- | --- |
| [preprocessing.py](../preprocessing.py) | Feature definitions, normalization, trajectory preparation, gap splitting, matrices, and row mapping. |
| [models.py](../models.py) | Multi-initialization Gaussian HMM fitting and fit metadata. |
| [state_annotation.py](../state_annotation.py) | Modern HMM annotation, state ordering, and legacy precomputed-feature interface. |
| [BCPA.py](../BCPA.py) | PELT/RBF boundaries, segment descriptors, K-means grouping, and annotation. |
| [hybrid_annotation.py](../hybrid_annotation.py) | BCPA segment construction and aggregation of HMM states/posteriors. |
| [criteria.py](../criteria.py) | Rule-annotation adapter, semantic/numeric state mapping, and common outputs. |
| [criteria core](../cpp/segmentation/csv_segmentation.cpp) | Predicate construction, time sorting, duplicate handling, and CSV/row adaptation. |
| [criterion classes](../cpp/segmentation/inw_monotone_crit.h) and [implementations](../cpp/segmentation/inw_monotone_crit.cpp) | Speed and angular ranges, outlier variants, disk confinement, and interval conventions. |
| [duration/hop criteria](../cpp/segmentation/outw_monotone_crit.cpp) | Minimum-duration and minimum-observation-step conditions. |
| [start-stop operations](../cpp/segmentation/startstop.cpp) | Boolean operations on admissible intervals. |
| [partition solver](../cpp/segmentation/brownian_criterion.cpp) | Global criteria segmentation; also contains separate internal Brownian-bridge routines. |
| [annotation.py](../annotation.py) | AnnotationResult, shared result handling, state summaries, and descriptive diagnostics. |
| [configuration.py](../configuration.py) | Common configuration and method builder. |
| [plot.py](../plot.py) | State and segment visualization. |

External references describe theory and dependencies. Source-linked implementation statements describe this checkout; they should be rechecked if the algorithms or dependency versions change.
