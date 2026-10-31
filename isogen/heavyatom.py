"""Estimate heavy-isotope enrichment from a measured isotope distribution.

This is the inverse of ``isodist(..., heavy_isotopes={...})``. The forward model
replaces every atom of one element with the mixture
``(1 - f) * natural + f * heavy``. The distribution of the molecule is then the
convolution of the n enriched atoms with the fixed distribution of all other
elements, so each peak is a degree-n polynomial in the single unknown ``f``
(the mass-isotopomer model of Kelleher & Masterson 1992 and Hellerstein &
Neese 1992). Distributions are reported relative to the base peak, so the
intensity scale is fitted together with ``f`` by bounded least squares
(``scipy.optimize.least_squares``) with an analytic Jacobian.
"""

import functools
import re

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import convolve

if __package__:
    from .isogen import isodist
    from .isogenwrapper import atom_formula_to_vector, isotope_mass
    from .mass import atom_masses_monoisotopic
    from .protein_mods import _ELEMENTS
else:
    from isogen import isodist
    from isogenwrapper import atom_formula_to_vector, isotope_mass
    from mass import atom_masses_monoisotopic
    from protein_mods import _ELEMENTS

# Heavy isotope assumed when only an element symbol is given.
_DEFAULT_HEAVY = {"H": 2, "C": 13, "N": 15, "O": 18, "S": 34}
_GRID_POINTS = 21


class HeavyPercentResult(tuple):
    """Result of :func:`calc_heavy_percent`.

    Unpacks as ``(percent, fit_dist)``. Extra information is available as
    attributes: ``percent``, ``fit_dist``, ``fraction``, ``stderr`` (standard
    error of the percentage, ``nan`` if it cannot be estimated), ``scale``
    (observed intensity per unit of ``fit_dist``), ``rmse`` (root-mean-square
    residual on the base-peak-normalized scale), ``heavy_atom`` (for example
    ``"2H"``) and ``n_atoms`` (number of atoms of the element in the formula).
    """

    def __new__(cls, percent, fit_dist, **info):
        self = super().__new__(cls, (percent, fit_dist))
        self.percent = percent
        self.fit_dist = fit_dist
        self.fraction = percent / 100.0
        for key, value in info.items():
            setattr(self, key, value)
        return self

    def __getnewargs_ex__(self):
        return (self.percent, self.fit_dist), dict(
            stderr=self.stderr, scale=self.scale, rmse=self.rmse,
            heavy_atom=self.heavy_atom, n_atoms=self.n_atoms)

    def __repr__(self):
        return ("HeavyPercentResult(percent={:.6g}, stderr={:.3g}, rmse={:.3g}, "
                "heavy_atom={!r})").format(self.percent, self.stderr, self.rmse, self.heavy_atom)


def _parse_heavy_atom(heavy_atom):
    if not isinstance(heavy_atom, str):
        raise TypeError("heavy_atom must be a string such as '2H' or '13C'")
    text = heavy_atom.strip()
    if text == "D":
        text = "2H"
    if text in _DEFAULT_HEAVY:
        text = "{}{}".format(_DEFAULT_HEAVY[text], text)
    match = re.fullmatch(r"(\d+)([A-Z][a-z]?)", text)
    if match is None:
        raise ValueError(
            "Invalid heavy atom {!r}; use a mass number and element symbol such as '13C'".format(heavy_atom))
    element = match.group(2)
    if element not in _ELEMENTS:
        raise ValueError("Unknown element in heavy atom {!r}".format(heavy_atom))
    return element, int(match.group(1))


def _convolve(a, b, length):
    """Convolution truncated or zero-padded to ``length``."""
    out = np.zeros(length)
    full = convolve(a, b)[:length]
    out[:full.size] = full
    return out


def _power(pmf, n, length):
    """``pmf`` convolved with itself ``n`` times, truncated to ``length``."""
    result = np.zeros(length)
    result[0] = 1.0
    base = pmf
    while n:
        if n & 1:
            result = _convolve(result, base, length)
        n >>= 1
        if n:
            base = _convolve(base, base, length)
    return result


@functools.lru_cache(maxsize=None)
def _atom_pmf(element, length):
    """Natural single-atom pmf indexed by mass shift from the lightest isotope."""
    dist = np.asarray(isodist(element, type="ATOM", isolen=max(length, 64), dist_only=True), dtype=float)
    pmf = dist[:length] / dist.sum()
    pmf.setflags(write=False)
    return pmf


def _prepare_observed(dist):
    arr = np.asarray(dist, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == 2:
        arr = arr[:, 1]
    if arr.ndim != 1 or arr.size < 2:
        raise ValueError("dist must be a 1-D intensity vector (or an (m, 2) array) with at least two values")
    if not np.all(np.isfinite(arr)):
        raise ValueError("dist must contain only finite values")
    peak = arr.max()
    if peak <= 0.0:
        raise ValueError("dist must contain a positive intensity")
    return arr / peak, peak


def calc_heavy_percent(dist, formula, heavy_atom, weights="none"):
    """Estimate the percent enrichment of a heavy isotope.

    This inverts ``isodist(formula, type="ATOM", heavy_isotopes={heavy_atom: f})``.
    The model is ``(1 - f) * natural + f * heavy`` for every atom of the heavy
    atom's element; the intensity scale is fitted and removed, so ``dist`` may be
    on any scale. All other elements keep their natural abundance.

    Args:
        dist: Measured intensities starting at the light-isotope (monoisotopic)
            peak and spaced by one nominal mass unit. A 1-D vector or an
            ``(m, 2)`` array of (mass, intensity) rows (only the intensity
            column is used). The length of the vector is the fitted window.
        formula: Elemental formula such as ``"C6H12O6"``.
        heavy_atom: Heavy isotope to test, such as ``"2H"``, ``"D"``, ``"13C"``,
            ``"15N"`` or ``"18O"``. A bare symbol (``"H"``, ``"C"``, ``"N"``,
            ``"O"``, ``"S"``) selects 2H, 13C, 15N, 18O or 34S.
        weights: ``"none"`` for ordinary least squares or ``"poisson"`` to weight
            each peak by the inverse of its intensity.

    Returns:
        A :class:`HeavyPercentResult` that unpacks as ``(percent, fit_dist)``.
        ``percent`` is the percentage (0-100) of the element's atoms that are
        the heavy isotope, and ``fit_dist`` is the base-peak-normalized
        distribution from ``isodist`` at that enrichment, with the same length
        as ``dist``.

    Raises:
        TypeError: If an argument has the wrong type.
        ValueError: If the distribution, formula, heavy atom or weights are
            invalid, if the element is not in the formula, or if the heavy
            isotope falls outside the distribution window.
    """
    if weights not in ("none", "poisson"):
        raise ValueError("weights must be 'none' or 'poisson'")
    element, mass_number = _parse_heavy_atom(heavy_atom)
    observed, peak = _prepare_observed(dist)
    length = observed.size
    counts = atom_formula_to_vector(formula)
    atomic_number = _ELEMENTS.index(element) + 1
    n_atoms = int(counts[atomic_number - 1])
    if n_atoms <= 0:
        raise ValueError("Heavy atom {!r} refers to an element not in the formula".format(heavy_atom))

    shift = int(round(isotope_mass(element, mass_number) - atom_masses_monoisotopic[atomic_number - 1]))
    if shift <= 0:
        raise ValueError("{!r} is not heavier than the lightest isotope of {}".format(heavy_atom, element))
    if shift >= length:
        raise ValueError("The mass shift of {!r} lies outside the supplied distribution".format(heavy_atom))

    natural = _atom_pmf(element, length)
    delta = -natural.copy()
    delta[shift] += 1.0  # d(pmf)/df = heavy - natural

    rest = np.zeros(length)
    rest[0] = 1.0
    for z in np.nonzero(counts)[0]:
        if z + 1 != atomic_number:
            rest = _convolve(rest, _power(_atom_pmf(_ELEMENTS[z], length), int(counts[z]), length), length)

    sqrt_w = np.ones(length) if weights == "none" else 1.0 / np.sqrt(np.maximum(observed, 1e-3))

    def model(f):
        pmf = (1.0 - f) * natural
        pmf[shift] += f
        q = _power(pmf, n_atoms - 1, length)
        return _convolve(rest, _convolve(q, pmf, length), length), q

    def projected_scale(m):
        wm = sqrt_w * m
        denominator = wm @ wm
        return float(wm @ (sqrt_w * observed) / denominator) if denominator > 0.0 else 0.0

    def residual(params):
        return sqrt_w * (observed - params[0] * model(params[1])[0])

    def jacobian(params):
        s, f = params
        m, q = model(f)
        dm = n_atoms * _convolve(rest, _convolve(q, delta, length), length)
        return -sqrt_w[:, None] * np.column_stack((m, s * dm))

    # Coarse scan of the 1-D variable-projection objective guards against local minima.
    grid = np.append(np.linspace(0.0, 1.0, _GRID_POINTS), _moment_estimate(observed, natural, rest, shift, n_atoms))
    costs = []
    for f in grid:
        m = model(f)[0]
        costs.append(np.sum((sqrt_w * (observed - projected_scale(m) * m)) ** 2))
    f0 = float(grid[int(np.argmin(costs))])
    s0 = projected_scale(model(f0)[0])

    fit = least_squares(residual, [s0, f0], jac=jacobian, bounds=([0.0, 0.0], [np.inf, 1.0]),
                        xtol=1e-14, ftol=1e-14, gtol=1e-14)
    scale, f = float(fit.x[0]), float(fit.x[1])
    m, _ = model(f)
    rmse = float(np.sqrt(np.mean((observed - scale * m) ** 2)))
    stderr = _stderr(fit.jac, fit.fun, length)

    label = "{}{}".format(mass_number, element)
    fit_dist = np.asarray(
        isodist(formula, type="ATOM", isolen=length, dist_only=True, heavy_isotopes={label: f}), dtype=float)
    return HeavyPercentResult(
        100.0 * f, fit_dist, stderr=100.0 * stderr, scale=scale * peak * float(m.max()), rmse=rmse,
        heavy_atom=label, n_atoms=n_atoms)


def _moment_estimate(observed, natural, rest, shift, n_atoms):
    """Closed-form mean-shift estimate of f; the mean is linear in f when untruncated."""
    index = np.arange(observed.size)
    observed_mean = observed @ index / observed.sum()
    atom_mean = natural @ index
    natural_mean = rest @ index / rest.sum() + n_atoms * atom_mean
    return float(np.clip((observed_mean - natural_mean) / (n_atoms * (shift - atom_mean)), 0.0, 1.0))


def _stderr(jac, residual, length):
    """Standard error of f from the Gauss-Newton covariance of (scale, f)."""
    if length <= 2:
        return float("nan")
    try:
        cov = np.linalg.inv(jac.T @ jac) * (residual @ residual) / (length - 2)
    except np.linalg.LinAlgError:
        return float("nan")
    return float(np.sqrt(max(cov[1, 1], 0.0)))
