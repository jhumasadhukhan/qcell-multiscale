"""Manuscript figures from results.json."""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.labelsize": 8.5, "axes.titlesize": 9.2, "axes.titleweight": "bold",
    "xtick.labelsize": 7.8, "ytick.labelsize": 7.8, "legend.fontsize": 7.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "lines.linewidth": 1.5,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
})

INK = "#1b2838"
C = {"clock": "#9aa5b1", "closed_raw": "#d1495b",
     "closed_od": "#edae49", "closed_multiscale": "#00798c"}
LBL = {"clock": "Clock (open loop)", "closed_raw": "Closed, raw ratio",
       "closed_od": "Closed, OD-corrected", "closed_multiscale": "Closed, multiscale"}

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results", "results.json")
if not os.path.exists(RESULTS):            # fall back to the working directory
    RESULTS = "results.json"
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)
R = json.load(open(RESULTS))


# ---------------------------------------------------------------- Figure 1
def fig1_framework():
    """The coupled framework. Labels are kept short deliberately: this panel is
    a map of the two layers and their interface, not a statement of the model,
    and long expressions here collide with the coupling arrows."""
    fig, ax = plt.subplots(figsize=(7.1, 4.9))
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.9); ax.axis("off")

    def box(x, y, w, h, fc, ec, title, lines, fs=7.0, tfs=8.0, dy=0.265):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.07",
                                    fc=fc, ec=ec, lw=1.2))
        if title:
            ax.text(x + w / 2, y + h - 0.26, title, ha="center", va="top",
                    fontsize=tfs, fontweight="bold", color=INK)
        top = y + h - (0.62 if title else 0.34)
        for i, s_ in enumerate(lines):
            ax.text(x + 0.16, top - i * dy, s_, ha="left", va="top",
                    fontsize=fs, color=INK)

    box(0.10, 3.70, 4.05, 3.05, "#e8f4f6", "#00798c",
        "QUANTUM LAYER   (fs - ns)",
        [r"Coupling  $J=\kappa\mu_D\mu_A/4\pi\varepsilon_0 n^2R^3$",
         r"Golden rule  $k_{EET}=4\pi^2c\,J^2\mathcal{J}$",
         r"Efficiency  $E=k_{EET}/(k_{EET}+\tau_D^{-1})$",
         r"Orientation  $\kappa^2_{eff}(\theta_{rot}/\tau_D)$",
         r"Index  $n(\phi_P)$  (Lorentz-Lorenz)",
         r"Transport  $e^{-\mu_s(\lambda)L}$ ,  $\lambda_D<\lambda_A$",
         r"Observable  $\mathcal{R}=I_A/I_D$"])

    box(5.90, 3.70, 4.05, 3.05, "#fdf3e3", "#c8872c",
        "CELLULAR LAYER   (min - h)",
        [r"BES formate supply  $R_{HCOOH}(V_{app})$",
         r"Substrate / nitrogen  $\dot C_S,\ \dot C_N$",
         r"Biomass / polymer  $\dot C_X,\ \dot C_P$",
         r"Haldane gates on $C_S$ and $C_N$",
         r"2-oxoglutarate node  $\dot G$",
         r"Active fraction  $\dot\phi_a=-b_{dec}(C_N)\phi_a+\ldots$",
         r"Granule loading  $\phi_P(C_X,C_P)$"])

    # coupling arrows, with their variable lists stacked so they clear the boxes
    ax.add_patch(FancyArrowPatch((5.85, 5.95), (4.20, 5.95),
                                 arrowstyle="-|>", mutation_scale=15,
                                 lw=1.8, color="#c8872c"))
    ax.text(5.0, 6.12, "downward", ha="center", fontsize=7.4,
            fontweight="bold", color="#c8872c")
    ax.text(5.0, 5.72, r"$G,\ \phi_P$", ha="center", fontsize=7.0, color=INK)
    ax.text(5.0, 5.42, r"$C_X,\ C_P,\ T$", ha="center", fontsize=7.0, color=INK)

    ax.add_patch(FancyArrowPatch((4.20, 4.62), (5.85, 4.62),
                                 arrowstyle="-|>", mutation_scale=15,
                                 lw=1.8, color="#00798c"))
    ax.text(5.0, 4.79, "upward", ha="center", fontsize=7.4,
            fontweight="bold", color="#00798c")
    ax.text(5.0, 4.39, r"$\mathcal{R}\rightarrow\hat G$", ha="center",
            fontsize=7.0, color=INK)
    ax.text(5.0, 4.09, r"burden $\Xi$", ha="center", fontsize=7.0, color=INK)

    box(1.05, 1.72, 7.9, 1.62, "#eef0f3", "#55606e",
        "ESTIMATOR  +  CONTROLLER",
        [r"Inverse map  $\hat G=h^{-1}(\mathcal{R};\ \phi_P,C_X,C_P)$  - undoes transport, $n$ and $\kappa^2$ drift",
         r"Nitrogen law  $F_N=\mathrm{sat}[-K_c((G^*-\hat G)+\frac{1}{\tau_I}\!\int\!(G^*-\hat G)dt)]$",
         r"Carbon substrate-stat on $C_S$   -   BES potential on reservoir level $M_{res}$"],
        dy=0.33)

    ax.add_patch(FancyArrowPatch((2.5, 3.64), (2.5, 3.40), arrowstyle="-|>",
                                 mutation_scale=13, lw=1.5, color="#00798c"))
    ax.add_patch(FancyArrowPatch((7.4, 1.66), (7.4, 1.22), arrowstyle="-|>",
                                 mutation_scale=13, lw=1.5, color="#55606e"))
    ax.text(7.55, 1.42, r"$F_N,\ F_S,\ V_{app}$", fontsize=7.0, color=INK)

    box(1.05, 0.15, 7.9, 1.00, "#f7f8f9", "#55606e", "",
        [r"Singular perturbation:  $\epsilon=\tau_Q/\tau_C\sim10^{-12}$  $\Rightarrow$  the quantum layer is solved",
         r"as an algebraic constraint $h(x,u)=0$ at each cellular step (adiabatic elimination)"],
        fs=7.0, dy=0.33)

    fig.savefig(os.path.join(FIGDIR, "fig1_framework.png"))
    plt.close(fig)


def fig2_quantum():
    d = np.array(R["drift"])
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.35))

    G = np.array(R["sensor_curve"]["G"])
    ax = axs[0]
    for k, col in zip(["0.00", "0.15", "0.30"], ["#00798c", "#edae49", "#d1495b"]):
        ax.semilogx(G, R["sensor_curve"]["curves"][k], color=col,
                    label=rf"$\phi_P={float(k):.2f}$")
    ax.axvline(0.16, color=INK, ls=":", lw=1.0)
    ax.text(0.17, ax.get_ylim()[0] + 0.05, "$G^*$", fontsize=7.2, color=INK)
    ax.set_xlabel("2-oxoglutarate, $G$ (mM)")
    ax.set_ylabel(r"ratio  $\mathcal{R}=I_A/I_D$")
    ax.set_title("(a) Sensor response")
    ax.legend(frameon=False, loc="upper left")

    ax = axs[1]
    ax.plot(d[:, 0], d[:, 2], color="#00798c", label=r"$n_{cyt}$")
    ax.set_xlabel("polymer fraction $f_P$ (g/g)")
    ax.set_ylabel(r"refractive index $n_{cyt}$", color="#00798c")
    ax.tick_params(axis="y", labelcolor="#00798c")
    ax2 = ax.twinx()
    ax2.plot(d[:, 0], d[:, 4], color="#d1495b")
    ax2.set_ylabel(r"$R_0$ (Å)", color="#d1495b")
    ax2.tick_params(axis="y", labelcolor="#d1495b")
    ax2.spines["top"].set_visible(False)
    ax.set_title("(b) Granules shift the optics")

    ax = axs[2]
    ax.axhline(0, color=INK, lw=0.8)
    ax.plot(d[:, 0], d[:, 7], color=C["closed_raw"], label="raw ratio")
    ax.plot(d[:, 0], d[:, 8], color=C["closed_od"], label="OD-corrected")
    ax.plot(d[:, 0], d[:, 9], color=C["closed_multiscale"], label="multiscale")
    ax.set_xlabel("polymer fraction $f_P$ (g/g)")
    ax.set_ylabel(r"error in $\hat G$ (%)")
    ax.set_title("(c) Calibration error")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_quantum.png"))
    plt.close(fig)


# ---------------------------------------------------------------- Figure 3
def fig3_trajectories():
    tr = R["trajectories"]
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.2))
    order = ["clock", "closed_raw", "closed_od", "closed_multiscale"]

    def get(m):
        t = np.array(tr[m]["t"]); y = np.array(tr[m]["y"]); a = np.array(tr[m]["aux"])
        return t, y, a

    ax = axs[0, 0]
    for m in order:
        t, y, a = get(m)
        ax.plot(t, y[:, 2], color=C[m], label=LBL[m])
    ax.set_ylabel("$C_X$ (g/L)"); ax.set_title("(a) Biomass")

    ax = axs[0, 1]
    for m in order:
        t, y, a = get(m)
        ax.plot(t, y[:, 3], color=C[m])
    ax.set_ylabel("$C_P$ (g/L)"); ax.set_title("(b) Polymer")

    ax = axs[0, 2]
    for m in order:
        t, y, a = get(m)
        ax.plot(t, y[:, 8], color=C[m])
    ax.set_ylabel(r"$\phi_a$"); ax.set_title("(c) Active fraction")

    ax = axs[1, 0]
    for m in order:
        t, y, a = get(m)
        ax.semilogy(t, np.maximum(y[:, 4], 1e-3), color=C[m])
    ax.axhline(0.16, color=INK, ls=":", lw=1.0)
    ax.set_ylabel("$G$ (mM)"); ax.set_title("(d) 2-oxoglutarate")

    ax = axs[1, 1]
    for m in order[1:]:
        t, y, a = get(m)
        msk = np.isfinite(a[:, 5]) & (y[:, 4] > 1e-4)
        ax.plot(t[msk], 100 * (a[msk, 5] - y[msk, 4]) / y[msk, 4], color=C[m])
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_ylabel(r"error in $\hat G$ (%)"); ax.set_title("(e) Estimator error")

    ax = axs[1, 2]
    for m in order:
        t, y, a = get(m)
        ax.plot(t, a[:, 0], color=C[m])
    ax.set_ylabel(r"$\phi_P$"); ax.set_title("(f) Granule loading")

    for ax in axs.ravel():
        ax.set_xlabel("time (h)")
    # One legend for the whole grid, in a reserved strip under the panels.
    # Inside panel (a) it sat on top of the biomass curve, and there is no
    # corner of that panel the sigmoid does not pass through.
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.legend(handles, labels, frameon=False, fontsize=6.8, ncol=4,
               loc="lower center", bbox_to_anchor=(0.5, -0.004),
               handlelength=1.8, columnspacing=1.6)
    fig.savefig(os.path.join(FIGDIR, "fig3_trajectories.png"))
    plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def fig4_performance():
    fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.4))

    sc = np.array(R["switch_scan"])
    ax = axs[0]
    ax.plot(sc[:, 0], sc[:, 4], color=INK, marker="o", ms=2.6)
    k = int(np.argmax(sc[:, 4]))
    ax.plot(sc[k, 0], sc[k, 4], "o", ms=6, mfc="none", mec="#d1495b", mew=1.5)
    ax.set_xlabel("nitrogen switch time (h)"); ax.set_ylabel("objective $J$ (g/L)")
    ax.set_title("(a) Open-loop optimum is sharp")

    sp = np.array(R["setpoint_scan"])
    ax = axs[1]
    ax.plot(sp[:, 0], sp[:, 4], color="#00798c", marker="o", ms=2.6)
    k = int(np.argmax(sp[:, 4]))
    ax.plot(sp[k, 0], sp[k, 4], "o", ms=6, mfc="none", mec="#d1495b", mew=1.5)
    ax.set_xlabel("2-OG setpoint $G^*$ (mM)"); ax.set_ylabel("objective $J$ (g/L)")
    ax.set_title("(b) Closed-loop setpoint")

    mc = R["monte_carlo"]
    ax = axs[2]
    keys = ["clock"] + ["closed_raw", "closed_od", "closed_multiscale"]
    data = [np.array(mc[k]) for k in keys]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False)
    for patch, k in zip(bp["boxes"], keys):
        patch.set_facecolor(C[k]); patch.set_alpha(0.75); patch.set_edgecolor(INK)
    for med in bp["medians"]:
        med.set_color(INK); med.set_linewidth(1.3)
    ax.set_xticks(range(1, len(keys) + 1))
    ax.set_xticklabels(["clock", "raw", "OD", "multi"], rotation=0)
    ax.set_ylabel("objective $J$ (g/L)")
    ax.set_title("(c) 60 perturbed strains")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig4_performance.png"))
    plt.close(fig)


# ---------------------------------------------------------------- Figure 5
def fig5_design():
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.7))

    ss = np.array(R["scattering"])
    ax = axs[0]
    ax.plot(ss[:, 0], ss[:, 1], color=C["closed_raw"], marker="o", ms=3, label="raw ratio")
    ax.plot(ss[:, 0], ss[:, 4], color=C["closed_od"], marker="s", ms=3, label="OD-corrected")
    ax.plot(ss[:, 0], ss[:, 7], color=C["closed_multiscale"], marker="^", ms=3,
            label="multiscale")
    ax.axvline(0.12, color=INK, ls=":", lw=1.0)
    ax.text(0.13, ax.get_ylim()[0] + 0.1, "nominal", fontsize=7, color=INK)
    ax.set_xlabel(r"granule scattering contrast  $a_P$ (L g$^{-1}$cm$^{-1}$)")
    ax.set_ylabel("objective $J$ (g/L)")
    ax.set_title("(a) When the quantum layer is needed")
    ax.legend(frameon=False, loc="lower left")

    dm = R["design_map"]
    KD = np.array(dm["KD"]); Rh = np.array(dm["Rh"]); S = np.array(dm["S"])
    ax = axs[1]
    pc = ax.pcolormesh(KD, Rh, S, cmap="viridis", shading="auto")
    ax.set_xscale("log")
    cs = ax.contour(KD, Rh, S, levels=6, colors="w", linewidths=0.5, alpha=0.6)
    ax.clabel(cs, inline=True, fontsize=5.5, fmt="%.2f")
    i, j = np.unravel_index(np.argmax(S), S.shape)
    ax.plot(KD[j], Rh[i], "*", ms=12, mfc="#ffd166", mec=INK, mew=0.8)
    ax.set_xlabel(r"binding constant $K_D$ (mM)")
    ax.set_ylabel(r"bound-state separation $R_{holo}$ (Å)")
    ax.set_title("(b) Sensor design map")
    fig.colorbar(pc, ax=ax, label="fractional ratio span")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig5_design.png"))
    plt.close(fig)


if __name__ == "__main__":
    fig1_framework()
    fig2_quantum()
    fig3_trajectories()
    fig4_performance()
    fig5_design()
    print("figures written")
