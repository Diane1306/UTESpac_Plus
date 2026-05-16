# SSITC Quality-Flagging Algorithm Used in UTESpac

## 1. Purpose

The code generates quality flags for eddy-covariance fluxes using a Foken/EddyPro-like framework. Two versions of the quality flag are stored:

| Output variable | Meaning | Intended use |
|---|---|---|
| `*_SSITC_TEST` | Combined steady-state + integral turbulence characteristics flag | AmeriFlux-style reporting |
| `*_SS_ONLY_TEST` | Steady-state-only flag | Scientific analysis of slope-flow/canopy regimes where ITC assumptions break down |

The following fluxes are flagged:

| Flux | Variable tested |
|---|---|
| `TAU` | `u'w'` and `v'w'` |
| `H` | `w'θ_v'` |
| `LE` | `w'H2O_WPL'` |
| `FC` | `w'CO2_WPL'` |

---

# 2. Averaging Period and Subperiods

Each averaging period, usually 30 min, is divided into 5-min subperiods.

If the averaging period is:

```math
T = 30 \ \mathrm{min}
```

and the subperiod length is:

```math
T_s = 5 \ \mathrm{min}
```

then the number of subperiods is:

```math
N_s = \frac{T}{T_s} = 6
```

The code requires:

```math
\mathrm{mod}(T,T_s)=0
```

---

# 3. Steady-State Test

## 3.1 Full-period covariance

For a scalar or velocity variable $begin:math:text$x$end:math:text$, the full-period covariance with vertical velocity $begin:math:text$w$end:math:text$ is:

```math
\overline{w'x'} =
\frac{1}{N}
\sum_{i=1}^{N}
\left(w_i-\overline{w}\right)
\left(x_i-\overline{x}\right)
```

where:

```math
\overline{w}=
\frac{1}{N}\sum_{i=1}^{N}w_i
```

```math
\overline{x}=
\frac{1}{N}\sum_{i=1}^{N}x_i
```

The code uses population-style Reynolds averaging, i.e. normalization by $begin:math:text$N$end:math:text$, not $begin:math:text$N\-1$end:math:text$.

---

## 3.2 Subperiod covariance

For each 5-min subperiod $begin:math:text$k$end:math:text$, the covariance is:

```math
\overline{w'x'}_k =
\frac{1}{N_k}
\sum_{i=1}^{N_k}
\left(w_{i,k}-\overline{w}_k\right)
\left(x_{i,k}-\overline{x}_k\right)
```

The mean subperiod covariance is:

```math
\left<\overline{w'x'}_k\right>
=
\frac{1}{N_s}
\sum_{k=1}^{N_s}
\overline{w'x'}_k
```

---

## 3.3 Steady-state deviation

The relative steady-state deviation is:

```math
SS_{\mathrm{dev}}
=
\left|
\frac{
\left<\overline{w'x'}_k\right>
-
\overline{w'x'}
}{
\overline{w'x'}
}
\right|
\times 100
```

where:

| Symbol | Meaning |
|---|---|
| $begin:math:text$SS\_\{\\mathrm\{dev\}\}$end:math:text$ | steady-state deviation in percent |
| $begin:math:text$\\overline\{w\'x\'\}$end:math:text$ | full-period covariance |
| $begin:math:text$\\left\<\\overline\{w\'x\'\}\_k\\right\>$end:math:text$ | mean of subperiod covariances |

---

# 4. Steady-State Variables Tested

## 4.1 Momentum flux

For momentum, two covariance components are tested:

```math
SS_{\mathrm{dev},uw}
=
\left|
\frac{
\left<\overline{u'w'}_k\right>
-
\overline{u'w'}
}{
\overline{u'w'}
}
\right|
\times 100
```

```math
SS_{\mathrm{dev},vw}
=
\left|
\frac{
\left<\overline{v'w'}_k\right>
-
\overline{v'w'}
}{
\overline{v'w'}
}
\right|
\times 100
```

The code uses the worse of the two:

```math
SS_{\mathrm{dev},\tau}
=
\max
\left(
SS_{\mathrm{dev},uw},
SS_{\mathrm{dev},vw}
\right)
```

---

## 4.2 Sensible heat flux

```math
SS_{\mathrm{dev},H}
=
\left|
\frac{
\left<\overline{w'\theta_v'}_k\right>
-
\overline{w'\theta_v'}
}{
\overline{w'\theta_v'}
}
\right|
\times 100
```

---

## 4.3 Latent heat flux

The WPL-corrected water-vapor perturbation is:

```math
\rho_{v,\mathrm{WPL}}'
=
\rho_v'
+
\rho_{v,\mathrm{external}}'
```

In the code:

```math
H2O_{\mathrm{WPL}}'
=
H2O'
+
\rho_{v,\mathrm{external}}' \times 10^3
```

The steady-state deviation is:

```math
SS_{\mathrm{dev},LE}
=
\left|
\frac{
\left<\overline{w'H2O_{\mathrm{WPL}}'}_k\right>
-
\overline{w'H2O_{\mathrm{WPL}}'}
}{
\overline{w'H2O_{\mathrm{WPL}}'}
}
\right|
\times 100
```

---

## 4.4 CO2 flux

The WPL-corrected CO2 perturbation is:

```math
\rho_{c,\mathrm{WPL}}'
=
\rho_c'
+
\rho_{c,\mathrm{external}}'
```

The steady-state deviation is:

```math
SS_{\mathrm{dev},FC}
=
\left|
\frac{
\left<\overline{w'\rho_{c,\mathrm{WPL}}'}_k\right>
-
\overline{w'\rho_{c,\mathrm{WPL}}'}
}{
\overline{w'\rho_{c,\mathrm{WPL}}'}
}
\right|
\times 100
```

---

# 5. Integral Turbulence Characteristics Test

## 5.1 Current implementation

The code uses only the vertical-velocity ITC test:

```math
\frac{\sigma_w}{u_*}
```

Scalar ITC tests are not used.

This choice avoids imposing scalar similarity assumptions on heat, water vapor, and CO2 fluxes in forested slope-flow regimes.

---

## 5.2 Friction velocity

The code computes:

```math
u_*
=
\sqrt{
\sqrt{
\overline{u'w'}^2
+
\overline{v'w'}^2
}
}
```

because the stored momentum-flux variable is:

```math
\tau_{\mathrm{kin}}
=
\sqrt{
\overline{u'w'}^2
+
\overline{v'w'}^2
}
```

so:

```math
u_* = \sqrt{\tau_{\mathrm{kin}}}
```

---

## 5.3 Vertical velocity standard deviation

The vertical velocity standard deviation is:

```math
\sigma_w =
\sqrt{
\frac{1}{N}
\sum_{i=1}^{N}
\left(w_i-\overline{w}\right)^2
}
```

For strict Reynolds-averaging consistency, this should use population normalization $begin:math:text$N$end:math:text$, not $begin:math:text$N\-1$end:math:text$.

In MATLAB, this corresponds to:

```matlab
std(wPF_P, 1, 'omitmissing')
```

---

## 5.4 Measured ITC quantity

```math
ITC_{\mathrm{measured}}
=
\frac{\sigma_w}{u_*}
```

---

# 6. Stability Parameter

For above-canopy ITC, the stability parameter is:

```math
\zeta =
\frac{z-d}{L}
```

where:

| Symbol | Meaning |
|---|---|
| $begin:math:text$z$end:math:text$ | sonic measurement height |
| $begin:math:text$d$end:math:text$ | displacement height |
| $begin:math:text$L$end:math:text$ | Obukhov length |

The code uses:

```math
\zeta =
\frac{
\max(z-d,0.1)
}{
L
}
```

The lower bound of 0.1 m prevents negative or zero effective height.

For current FM Dolly use, recommended configuration is:

```math
d = 0
```

because displacement height is uncertain in steep, heterogeneous forested terrain.

---

# 7. Above-Canopy ITC Model

For above-canopy measurements, the expected value of $begin:math:text$\\sigma\_w\/u\_\*$end:math:text$ is modeled as:

```math
\left(\frac{\sigma_w}{u_*}\right)_{\mathrm{model}}
=
c_1 |\zeta|^{c_2}
```

The code currently uses:

| Stability regime | Condition | $begin:math:text$c\_1$end:math:text$ | $begin:math:text$c\_2$end:math:text$ | Model |
|---|---:|---:|---:|---|
| Unstable | $begin:math:text$\\zeta \< \-0\.032$end:math:text$ | 2.0 | 1/8 | $begin:math:text$2\.0\|\\zeta\|\^\{1\/8\}$end:math:text$ |
| Near-neutral unstable | $begin:math:text$\-0\.032 \\le \\zeta \< 0$end:math:text$ | 1.3 | 0 | $begin:math:text$1\.3$end:math:text$ |
| Stable fallback | $begin:math:text$\\zeta \\ge 0$end:math:text$ | 1.3 | 0 | $begin:math:text$1\.3$end:math:text$ |

Thus:

```math
\left(\frac{\sigma_w}{u_*}\right)_{\mathrm{model}}
=
2.0|\zeta|^{1/8},
\quad \zeta < -0.032
```

```math
\left(\frac{\sigma_w}{u_*}\right)_{\mathrm{model}}
=
1.3,
\quad \zeta \ge -0.032
```

The stable case is treated conservatively using the near-neutral value because stable ITC assumptions are especially uncertain in forested complex terrain.

---

# 8. Canopy-Aware ITC Model

For within-canopy measurements, where:

```math
z \le h_c
```

the code uses a canopy-aware parameterization:

```math
\left(\frac{\sigma_w}{u_*}\right)_{\mathrm{model}}
=
a_i
\left[
\exp
\left(
-\alpha_i
\left(
1-\frac{z}{h_c}
\right)^{\beta_i}
\right)
(1-\gamma_i)
+
\gamma_i
\right]
```

where:

| Parameter | Value |
|---|---:|
| $begin:math:text$a\_i$end:math:text$ | 1.25 |
| $begin:math:text$\\alpha\_i$end:math:text$ | 0.9 |
| $begin:math:text$\\beta\_i$end:math:text$ | 1.2 |
| $begin:math:text$\\gamma\_i$end:math:text$ | -0.63 |

and:

| Symbol | Meaning |
|---|---|
| $begin:math:text$z$end:math:text$ | sonic height |
| $begin:math:text$h\_c$end:math:text$ | canopy height |
| $begin:math:text$z\/h\_c$end:math:text$ | normalized canopy height |

The normalized height is constrained as:

```math
0 \le \frac{z}{h_c} \le 1
```

---

# 9. ITC Deviation

The ITC deviation is:

```math
ITC_{\mathrm{dev}}
=
\left|
\frac{
ITC_{\mathrm{model}}
-
ITC_{\mathrm{measured}}
}{
ITC_{\mathrm{model}}
}
\right|
\times 100
```

or:

```math
ITC_{\mathrm{dev}}
=
\left|
\frac{
\left(\sigma_w/u_*\right)_{\mathrm{model}}
-
\left(\sigma_w/u_*\right)_{\mathrm{measured}}
}{
\left(\sigma_w/u_*\right)_{\mathrm{model}}
}
\right|
\times 100
```

---

# 10. SSITC Final Flag

The combined SSITC flag is:

| Condition | `*_SSITC_TEST` |
|---|---:|
| $begin:math:text$SS\_\{\\mathrm\{dev\}\} \< 30\\\%$end:math:text$ and $begin:math:text$ITC\_\{\\mathrm\{dev\}\} \< 30\\\%$end:math:text$ | 0 |
| $begin:math:text$SS\_\{\\mathrm\{dev\}\} \< 100\\\%$end:math:text$ and $begin:math:text$ITC\_\{\\mathrm\{dev\}\} \< 100\\\%$end:math:text$ | 1 |
| otherwise | 2 |

Thus:

```math
SSITC =
\begin{cases}
0, & SS_{\mathrm{dev}} < 30\% \ \mathrm{and} \ ITC_{\mathrm{dev}} < 30\% \\
1, & SS_{\mathrm{dev}} < 100\% \ \mathrm{and} \ ITC_{\mathrm{dev}} < 100\% \\
2, & \mathrm{otherwise}
\end{cases}
```

---

# 11. Steady-State-Only Final Flag

The steady-state-only flag is:

| Condition | `*_SS_ONLY_TEST` |
|---|---:|
| $begin:math:text$SS\_\{\\mathrm\{dev\}\} \< 30\\\%$end:math:text$ | 0 |
| $begin:math:text$SS\_\{\\mathrm\{dev\}\} \< 100\\\%$end:math:text$ | 1 |
| otherwise | 2 |

Thus:

```math
SS_{\mathrm{only}}
=
\begin{cases}
0, & SS_{\mathrm{dev}} < 30\% \\
1, & SS_{\mathrm{dev}} < 100\% \\
2, & \mathrm{otherwise}
\end{cases}
```

---

# 12. Instrument and Diagnostic Overrides

If relevant instrument or diagnostic flags indicate invalid data, the final quality flag is forced to 2.

## 12.1 Momentum

```math
TAU_{\mathrm{flag}} = 2
```

if:

```math
\mathrm{rotatedSonFlag} = 1
```

---

## 12.2 Sensible heat

```math
H_{\mathrm{flag}} = 2
```

if:

```math
\mathrm{rotatedSonFlag} = 1
```

or:

```math
\mathrm{TsonFlag} = 1
```

---

## 12.3 Latent heat

```math
LE_{\mathrm{flag}} = 2
```

if:

```math
\mathrm{rotatedSonFlag} = 1
```

or:

```math
\mathrm{H2OFlag} = 1
```

---

## 12.4 CO2 flux

```math
FC_{\mathrm{flag}} = 2
```

if:

```math
\mathrm{rotatedSonFlag} = 1
```

or:

```math
\mathrm{CO2Flag} = 1
```

---

# 13. Output Variables

For each sonic height, the code stores eight columns:

| Output variable | Description |
|---|---|
| `TAU_SSITC_TEST` | momentum combined SS + ITC flag |
| `TAU_SS_ONLY_TEST` | momentum steady-state-only flag |
| `H_SSITC_TEST` | sensible heat combined SS + ITC flag |
| `H_SS_ONLY_TEST` | sensible heat steady-state-only flag |
| `LE_SSITC_TEST` | latent heat combined SS + ITC flag |
| `LE_SS_ONLY_TEST` | latent heat steady-state-only flag |
| `FC_SSITC_TEST` | CO2 flux combined SS + ITC flag |
| `FC_SS_ONLY_TEST` | CO2 flux steady-state-only flag |

Example output headers:

```matlab
'time'
'4.42m:TAU_SSITC_TEST'
'4.42m:TAU_SS_ONLY_TEST'
'4.42m:H_SSITC_TEST'
'4.42m:H_SS_ONLY_TEST'
'4.42m:LE_SSITC_TEST'
'4.42m:LE_SS_ONLY_TEST'
'4.42m:FC_SSITC_TEST'
'4.42m:FC_SS_ONLY_TEST'
```

---

# 14. Recommended Interpretation for FM Dolly

For AmeriFlux reporting:

```math
*_SSITC\_TEST
```

should be reported for transparency and compatibility.

For scientific analysis of slope-flow and within-canopy turbulence:

```math
*_SSITC\_TEST
```

should not be used as a strict exclusion filter because the analysis targets regimes where stationarity and integral turbulence similarity assumptions are expected to break down.

Instead, SSITC flags should be interpreted as diagnostic indicators of:

| Flag behavior | Physical interpretation |
|---|---|
| many flag-2 periods inside canopy | expected nonstationarity / canopy intermittency |
| ITC failure during slope flow | MOST or flux-variance similarity breakdown |
| SS failure during transitions | intermittent coherent transport or regime change |
| high-quality instrument diagnostics but poor SSITC | physically meaningful non-ideal turbulence |

Therefore, for FM Dolly:

```math
SSITC \neq \mathrm{instrument\ failure}
```

Rather:

```math
SSITC = \mathrm{diagnostic\ of\ stationarity\ and\ similarity\ assumptions}
```

---

# 15. Recommended Statement for Reports

Foken-style steady-state and integral turbulence characteristics flags were computed for AmeriFlux-compatible reporting. The combined SSITC flag uses the steady-state covariance test together with a canopy-aware, vertical-velocity-based ITC test. For within-canopy measurements, a Rannik-style canopy parameterization was used for $begin:math:text$\\sigma\_w\/u\_\*$end:math:text$, while above-canopy measurements used a Foken-style stability-dependent reference. Because the scientific analysis focuses on forested slope-flow regimes where stationarity and similarity assumptions are expected to break down, SSITC flags were not used as hard exclusion criteria. Instead, they were interpreted as diagnostic indicators of nonstationarity and departures from similarity theory.