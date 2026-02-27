import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.widgets import Slider
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection
from math import gcd
from functools import reduce


# ═══════════════════════════════════════════════════════════════════
#  MATRIX CONSTRUCTORS
# ═══════════════════════════════════════════════════════════════════

def dft_matrix(n):
    j, k = np.meshgrid(np.arange(n), np.arange(n), indexing='ij')
    return np.exp(2j * np.pi * j * k / n)


def sylvester_hadamard(m):
    H = np.array([[1.0 + 0j]])
    for _ in range(m):
        H = np.block([[H, H], [H, -H]])
    return H


def is_complex_hadamard(matrix, tol=1e-8):
    n = matrix.shape[0]
    if matrix.shape != (n, n):
        return False
    if not np.allclose(np.abs(matrix), 1, atol=tol):
        return False
    product = matrix @ matrix.conj().T
    return np.allclose(product, n * np.eye(n), atol=tol)


# ═══════════════════════════════════════════════════════════════════
#  LCM UTILITIES
# ═══════════════════════════════════════════════════════════════════

def lcm(a, b):
    return a * b // gcd(a, b)


def multi_lcm(values):
    return reduce(lcm, values)


# ═══════════════════════════════════════════════════════════════════
#  OVERLAY PRODUCT
# ═══════════════════════════════════════════════════════════════════

def overlay_product_fast(A, B):
    """Overlay two matrices on [0,1]², multiply where they overlap."""
    nA, nB = A.shape[0], B.shape[0]
    N = lcm(nA, nB)
    idx = np.arange(N)
    a_idx = idx * nA // N
    b_idx = idx * nB // N
    return A[np.ix_(a_idx, a_idx)] * B[np.ix_(b_idx, b_idx)]


def overlay_multi(*matrices):
    """Overlay any number of matrices on [0,1]²."""
    result = matrices[0]
    for M in matrices[1:]:
        result = overlay_product_fast(result, M)
    return result


# ═══════════════════════════════════════════════════════════════════
#  EFFICIENT MULTI-DFT OVERLAY
#
#  Instead of building huge intermediate matrices, we work directly
#  with the argument sum:
#
#    arg(result[i,j]) = Σ_{m=2}^{n}  2π · floor(i·m/N) · floor(j·m/N) / m
#
#  where N = lcm(2, 3, ..., n).
#  We accumulate angles as floats and only convert to RGB at the end.
# ═══════════════════════════════════════════════════════════════════

def fourier_overlay_angles(n_max, N=None):
    """
    Compute the argument matrix for DFT(2) ◇ DFT(3) ◇ ... ◇ DFT(n_max).

    Returns (angles, N) where angles is N×N float array of total arguments
    and N = lcm(2, ..., n_max).

    For large N we use a resolution cap and sample uniformly.
    """
    dims = list(range(2, n_max + 1))
    true_N = multi_lcm(dims)

    # Cap resolution for display — beyond ~4000 pixels you can't see more
    MAX_RES = 4000
    if N is None:
        N = min(true_N, MAX_RES)

    angles = np.zeros((N, N), dtype=np.float64)

    for m in dims:
        idx = np.arange(N)
        # Map pixel i to row index in the m×m DFT
        m_idx = (idx * m // N) if N == true_N else (idx * m / N).astype(int) % m
        # DFT(m) entry at (r, c) has argument 2π·r·c/m
        rows = m_idx[:, None]  # (N, 1)
        cols = m_idx[None, :]  # (1, N)
        angles += 2 * np.pi * (rows * cols) / m

    return angles, N, true_N


def fourier_overlay_rgb(n_max, N=None):
    """Compute the RGB image for the first n_max Fourier overlay."""
    angles, N_used, true_N = fourier_overlay_angles(n_max, N)
    hues = (angles / (2 * np.pi)) % 1.0
    ones = np.ones_like(hues)
    hsv  = np.stack([hues, ones, ones], axis=-1)
    rgb  = mcolors.hsv_to_rgb(hsv)
    return rgb, N_used, true_N


# ═══════════════════════════════════════════════════════════════════
#  arg(z) → HUE COLORING
# ═══════════════════════════════════════════════════════════════════

def complex_to_rgb(matrix):
    angles = np.angle(matrix)
    hues   = (angles / (2 * np.pi)) % 1.0
    ones   = np.ones_like(hues)
    hsv    = np.stack([hues, ones, ones], axis=-1)
    return mcolors.hsv_to_rgb(hsv)


# ═══════════════════════════════════════════════════════════════════
#  COLOUR-WHEEL LEGEND
# ═══════════════════════════════════════════════════════════════════

def draw_color_wheel(ax, resolution=360):
    dtheta = 2 * np.pi / resolution
    patches, colors = [], []
    for i in range(resolution):
        t = i * dtheta
        patches.append(
            Wedge((0, 0), 1.0, np.degrees(t), np.degrees(t + dtheta), width=0.35)
        )
        colors.append(mcolors.hsv_to_rgb([(t / (2 * np.pi)) % 1, 1, 1]))
    ax.add_collection(
        PatchCollection(patches, facecolors=colors, edgecolors='none')
    )
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect('equal')
    labels = {0: '0', np.pi/2: 'π/2', np.pi: 'π', 3*np.pi/2: '3π/2'}
    for theta, txt in labels.items():
        ax.text(1.22 * np.cos(theta), 1.22 * np.sin(theta), txt,
                ha='center', va='center', fontsize=10)
    ax.set_title('arg(z) → hue', fontsize=11, pad=10)
    ax.axis('off')


# ═══════════════════════════════════════════════════════════════════
#  SINGLE MATRIX VISUALIZER WITH SLIDER
# ═══════════════════════════════════════════════════════════════════

class HadamardVisualizer:
    def __init__(self, matrix_func=dft_matrix,
                 n_min=2, n_max=512, initial_n=16,
                 title="DFT Matrix"):
        self.matrix_func = matrix_func
        self.title = title

        self.fig = plt.figure(figsize=(10, 9))
        self.ax_img   = self.fig.add_axes([0.08, 0.13, 0.62, 0.78])
        self.ax_wheel = self.fig.add_axes([0.76, 0.56, 0.22, 0.34])
        self.ax_info  = self.fig.add_axes([0.76, 0.15, 0.22, 0.30])
        self.ax_slide = self.fig.add_axes([0.08, 0.03, 0.62, 0.03])
        self.ax_info.axis('off')

        self.slider = Slider(self.ax_slide, 'n', n_min, n_max,
                             valinit=initial_n, valstep=1)
        self.slider.on_changed(self._update)
        draw_color_wheel(self.ax_wheel)
        self._update(initial_n)

    def _update(self, n):
        n = int(n)
        M   = self.matrix_func(n)
        rgb = complex_to_rgb(M)
        ax = self.ax_img
        ax.clear()
        ax.imshow(rgb, origin='lower', extent=[0, 1, 0, 1],
                  interpolation='nearest', aspect='equal')
        ax.set_title(f'{self.title}   n = {n}', fontsize=14)
        ax.set_xlabel('column  k / n', fontsize=11)
        ax.set_ylabel('row  j / n', fontsize=11)

        self.ax_info.clear()
        self.ax_info.axis('off')
        hadamard_ok = is_complex_hadamard(M)
        info = (
            f"dim = {n} x {n}\n"
            f"Hadamard? {'Y' if hadamard_ok else 'N'}\n\n"
            f"Corners (arg / pi):\n"
            f" (0,0):  {np.angle(M[0,0])/np.pi:+.4f}\n"
            f" (n-1,n-1): {np.angle(M[n-1,n-1])/np.pi:+.4f}"
        )
        self.ax_info.text(0.05, 0.95, info, transform=self.ax_info.transAxes,
                          fontsize=10, verticalalignment='top',
                          family='monospace',
                          bbox=dict(boxstyle='round', facecolor='#f0f0f0'))
        self.fig.canvas.draw_idle()

    def overlay(self, plot_func):
        plot_func(self.ax_img)
        self.fig.canvas.draw_idle()

    def show(self):
        plt.show()


# ═══════════════════════════════════════════════════════════════════
#  OVERLAY VISUALIZER — TWO MATRICES
# ═══════════════════════════════════════════════════════════════════

class OverlayVisualizer:
    def __init__(self,
                 matrix_func_A=dft_matrix,
                 matrix_func_B=dft_matrix,
                 n_min=2, n_max=30,
                 initial_nA=3, initial_nB=4,
                 title_A="DFT(A)", title_B="DFT(B)"):

        self.mfA, self.mfB = matrix_func_A, matrix_func_B
        self.title_A, self.title_B = title_A, title_B

        self.fig = plt.figure(figsize=(14, 10))
        self.ax_A     = self.fig.add_axes([0.04, 0.55, 0.28, 0.38])
        self.ax_B     = self.fig.add_axes([0.36, 0.55, 0.28, 0.38])
        self.ax_wheel = self.fig.add_axes([0.72, 0.62, 0.18, 0.28])
        self.ax_combo = self.fig.add_axes([0.04, 0.08, 0.54, 0.42])
        self.ax_info  = self.fig.add_axes([0.64, 0.08, 0.34, 0.42])
        self.ax_sA    = self.fig.add_axes([0.04, 0.51, 0.28, 0.02])
        self.ax_sB    = self.fig.add_axes([0.36, 0.51, 0.28, 0.02])
        self.ax_info.axis('off')

        self.slider_A = Slider(self.ax_sA, 'n_A', n_min, n_max,
                               valinit=initial_nA, valstep=1)
        self.slider_B = Slider(self.ax_sB, 'n_B', n_min, n_max,
                               valinit=initial_nB, valstep=1)
        self.slider_A.on_changed(self._update)
        self.slider_B.on_changed(self._update)

        draw_color_wheel(self.ax_wheel)
        self._update(None)

    def _update(self, _):
        nA = int(self.slider_A.val)
        nB = int(self.slider_B.val)
        A  = self.mfA(nA)
        B  = self.mfB(nB)
        N  = lcm(nA, nB)
        C  = overlay_product_fast(A, B)

        self.ax_A.clear()
        self.ax_A.imshow(complex_to_rgb(A), origin='lower',
                         extent=[0, 1, 0, 1], interpolation='nearest',
                         aspect='equal')
        self.ax_A.set_title(f'{self.title_A}  n={nA}', fontsize=11)

        self.ax_B.clear()
        self.ax_B.imshow(complex_to_rgb(B), origin='lower',
                         extent=[0, 1, 0, 1], interpolation='nearest',
                         aspect='equal')
        self.ax_B.set_title(f'{self.title_B}  n={nB}', fontsize=11)

        self.ax_combo.clear()
        self.ax_combo.imshow(complex_to_rgb(C), origin='lower',
                             extent=[0, 1, 0, 1],
                             interpolation='nearest', aspect='equal')
        self.ax_combo.set_title(
            f'Overlay   {nA} ◇ {nB} → {N}×{N}', fontsize=13)

        # Grid lines
        for t in np.linspace(0, 1, nA + 1):
            self.ax_combo.axhline(t, color='white', alpha=0.4, lw=0.7)
            self.ax_combo.axvline(t, color='white', alpha=0.4, lw=0.7)
        for t in np.linspace(0, 1, nB + 1):
            self.ax_combo.axhline(t, color='black', alpha=0.4, lw=0.7,
                                  ls='--')
            self.ax_combo.axvline(t, color='black', alpha=0.4, lw=0.7,
                                  ls='--')

        self.ax_info.clear()
        self.ax_info.axis('off')
        info = (
            f"A: {nA}x{nA}  Had? {'Y' if is_complex_hadamard(A) else 'N'}\n"
            f"B: {nB}x{nB}  Had? {'Y' if is_complex_hadamard(B) else 'N'}\n"
            f"\nlcm({nA},{nB}) = {N}\n"
            f"Result: {N}x{N}\n"
            f"  Had? {'Y' if is_complex_hadamard(C) else 'N'}\n"
            f"\nWhite  = A grid\n"
            f"Dashed = B grid"
        )
        self.ax_info.text(0.05, 0.95, info,
                          transform=self.ax_info.transAxes,
                          fontsize=10, verticalalignment='top',
                          family='monospace',
                          bbox=dict(boxstyle='round', facecolor='#f0f0f0'))
        self.fig.canvas.draw_idle()

    def show(self):
        plt.show()


# ═══════════════════════════════════════════════════════════════════
#  ★  CUMULATIVE FOURIER OVERLAY  — DFT(2) ◇ DFT(3) ◇ ... ◇ DFT(n)
# ═══════════════════════════════════════════════════════════════════

class FourierOverlayVisualizer:
    """
    Single slider for n.  Shows the overlay of DFT(2) through DFT(n)
    on the unit square.

    The true grid size is lcm(2,...,n), which explodes quickly:
        n=5  → 60
        n=7  → 420
        n=10 → 2520
        n=13 → 360360

    So we cap pixel resolution at MAX_RES and sample the continuous
    limit function:
        θ(x,y) = Σ_{m=2}^{n}  2π · floor(m·x) · floor(m·y) / m

    This gives the exact answer at rational points and a faithful
    picture everywhere.
    """

    MAX_RES = 3000   # pixels per side — bump up if you have RAM/patience

    def __init__(self, n_min=2, n_max=30, initial_n=5):
        self.n_min = n_min

        self.fig = plt.figure(figsize=(11, 9))
        self.ax_img   = self.fig.add_axes([0.06, 0.12, 0.62, 0.80])
        self.ax_wheel = self.fig.add_axes([0.74, 0.58, 0.22, 0.32])
        self.ax_info  = self.fig.add_axes([0.74, 0.12, 0.24, 0.40])
        self.ax_slide = self.fig.add_axes([0.06, 0.03, 0.62, 0.03])
        self.ax_info.axis('off')

        self.slider = Slider(self.ax_slide, 'n (overlay DFT 2..n)',
                             n_min, n_max,
                             valinit=initial_n, valstep=1)
        self.slider.on_changed(self._update)

        draw_color_wheel(self.ax_wheel)
        self._update(initial_n)

    def _compute(self, n_max):
        """
        θ(i,j) = Σ_{m=2}^{n_max}  2π · r_m(i) · r_m(j) / m

        where r_m(i) = floor(i · m / N) and N = display resolution.
        """
        dims = list(range(2, n_max + 1))
        true_N = multi_lcm(dims)
        N = min(true_N, self.MAX_RES)

        # Use float coordinates in [0, 1)
        t = np.linspace(0, 1, N, endpoint=False)

        angles = np.zeros((N, N), dtype=np.float64)
        for m in dims:
            # floor(m * t) gives the row/col index in the m×m DFT
            m_idx = np.floor(m * t).astype(np.int64)
            # Clamp just in case of floating point at boundary
            m_idx = np.clip(m_idx, 0, m - 1)
            rows = m_idx[:, None]  # (N,1)
            cols = m_idx[None, :]  # (1,N)
            angles += (2 * np.pi * rows * cols) / m

        return angles, N, true_N

    def _update(self, n):
        n = int(n)
        angles, N, true_N = self._compute(n)

        hues = (angles / (2 * np.pi)) % 1.0
        ones = np.ones_like(hues)
        hsv  = np.stack([hues, ones, ones], axis=-1)
        rgb  = mcolors.hsv_to_rgb(hsv)

        ax = self.ax_img
        ax.clear()
        ax.imshow(rgb, origin='lower', extent=[0, 1, 0, 1],
                  interpolation='nearest', aspect='equal')
        ax.set_title(f'DFT(2) ◇ DFT(3) ◇ ··· ◇ DFT({n})', fontsize=14)
        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('y', fontsize=12)

        self.ax_info.clear()
        self.ax_info.axis('off')

        dims_str = ' × '.join(str(m) for m in range(2, n + 1))
        info = (
            f"Overlaid: DFT(2) ... DFT({n})\n"
            f"  {n - 1} matrices total\n\n"
            f"True grid: {true_N} × {true_N}\n"
            f"  lcm(2,...,{n}) = {true_N}\n"
            f"Display:   {N} × {N} px\n\n"
            f"θ(x,y) = Σ 2π⌊mx⌋⌊my⌋/m\n"
            f"  summed over m = 2..{n}\n\n"
            f"Upper-right arg/π:\n"
            f"  {(angles[-1,-1]/np.pi) % 2:+.4f}"
        )
        self.ax_info.text(0.02, 0.98, info,
                          transform=self.ax_info.transAxes,
                          fontsize=9, verticalalignment='top',
                          family='monospace',
                          bbox=dict(boxstyle='round', facecolor='#f0f0f0'))
        self.fig.canvas.draw_idle()

    def show(self):
        plt.show()


# ═══════════════════════════════════════════════════════════════════
#  STATIC HELPERS
# ═══════════════════════════════════════════════════════════════════

def visualize_matrix(matrix, title="Complex Hadamard Matrix", ax=None):
    own_fig = (ax is None)
    if own_fig:
        fig, axes = plt.subplots(1, 2, figsize=(11, 7),
                                 gridspec_kw={'width_ratios': [3, 1]})
        ax, ax_w = axes
        draw_color_wheel(ax_w)

    rgb = complex_to_rgb(matrix)
    n   = matrix.shape[0]
    ax.imshow(rgb, origin='lower', extent=[0, 1, 0, 1],
              interpolation='nearest', aspect='equal')
    ax.set_title(f'{title}  (n = {n})', fontsize=14)
    ax.set_xlabel('column  k / n', fontsize=11)
    ax.set_ylabel('row  j / n', fontsize=11)

    if own_fig:
        plt.tight_layout()
        plt.show()
    return ax


# ═══════════════════════════════════════════════════════════════════
#  DEMO
# ═══════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    print("=" * 55)
    print("  Complex Hadamard Matrix Visualizer")
    print("=" * 55)
    print()
    print("  1 — Single DFT (dimension slider)")
    print("  2 — Overlay two DFTs (two sliders)")
    print("  3 — ★ Cumulative Fourier overlay:")
    print("        DFT(2) ◇ DFT(3) ◇ ··· ◇ DFT(n)")
    print("  4 — Static: DFT(3) ◇ DFT(4) ◇ DFT(5) → 60×60")
    print("  5 — Compare: single DFTs side by side")
    print()

    choice = input("Enter 1-5: ").strip()

    if choice == '1':
        HadamardVisualizer(
            matrix_func=dft_matrix,
            n_min=2, n_max=2187,
            initial_n=16,
            title="DFT Matrix"
        ).show()

    elif choice == '2':
        OverlayVisualizer(
            matrix_func_A=dft_matrix,
            matrix_func_B=dft_matrix,
            n_min=2, n_max=60,
            initial_nA=3, initial_nB=4,
            title_A="DFT(A)", title_B="DFT(B)"
        ).show()

    elif choice == '3':
        FourierOverlayVisualizer(
            n_min=2, n_max=60, initial_n=5
        ).show()

    elif choice == '4':
        C = overlay_multi(dft_matrix(3), dft_matrix(4), dft_matrix(5))
        visualize_matrix(C, title="DFT(3) ◇ DFT(4) ◇ DFT(5)")

    elif choice == '5':
        fig, axes = plt.subplots(1, 5, figsize=(20, 4.5),
                                 gridspec_kw={'width_ratios': [3,3,3,3,1]})
        for i, n in enumerate([4, 8, 16, 64]):
            M = dft_matrix(n)
            axes[i].imshow(complex_to_rgb(M), origin='lower',
                           extent=[0, 1, 0, 1], interpolation='nearest',
                           aspect='equal')
            axes[i].set_title(f'DFT({n})', fontsize=11)
        draw_color_wheel(axes[4])
        plt.tight_layout()
        plt.show()

    else:
        print("Running default: Cumulative Fourier overlay")
        FourierOverlayVisualizer(n_min=2, n_max=30, initial_n=5).show()