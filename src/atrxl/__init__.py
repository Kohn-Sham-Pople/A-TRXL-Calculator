"""A-TRXL neutral-atom pair coefficient calculator. q: Å⁻¹; energy: keV."""
from .factors import atomic_factor, factor_components
from .weights import pair_weight, delta_pair_weight, pair_weight_curve, delta_weight_curve
__version__ = "1.0.0"
