#!/usr/bin/env python3
"""
Ultra-fast redshift distribution forecaster for strong lenses
Generates deflector and source redshift distributions in minutes

Usage:
    python FastRedshiftForecast.py <survey> <N_target_lenses> [output_file]

Example:
    python FastRedshiftForecast.py Euclid 15000 euclid_redshift_forecast.pkl
"""

import numpy as np
import pickle
import time
import sys
from scipy import interpolate
import matplotlib.pyplot as plt

# Import lenspop components
from PopulationFunctions import LensPopulation_, SourcePopulation_, EinsteinRadiusTools
import distances

class FastLensForecast:
    """
    Vectorized lens forecasting for redshift distributions only
    Skips all image generation and disk I/O
    """

    def __init__(self, survey='Euclid', cosmo=[0.319, 0.681, 0.67],
                 sigfloor=100, zlmax=2, sourcepop='lsst'):
        """
        Initialize forecast

        Parameters:
        -----------
        survey : str
            Survey name (e.g., 'Euclid', 'LSST', 'DES')
        cosmo : list
            [Omega_m, Omega_Lambda, h]
        sigfloor : float
            Minimum velocity dispersion (km/s)
        zlmax : float
            Maximum lens redshift
        sourcepop : str
            Source population ('lsst' or 'jaguar')
        """
        print("="*70)
        print("FAST LENS REDSHIFT FORECAST")
        print("="*70)
        print(f"Survey: {survey}")
        print(f"Cosmology: Omega_m={cosmo[0]}, Omega_L={cosmo[1]}, h={cosmo[2]}")
        print(f"Source population: {sourcepop}")

        self.survey = survey
        self.cosmo = cosmo
        self.sigfloor = sigfloor
        self.zlmax = zlmax
        self.sourcepop = sourcepop

        # Initialize cosmology
        self.D = distances.Distance(cosmo=cosmo)

        # Initialize lens population
        print("\nInitializing lens population...")
        self.L = LensPopulation_(zlmax=zlmax, sigfloor=sigfloor, D=self.D, cosmo=cosmo)

        # Initialize source population
        print(f"Loading {sourcepop} source catalog...")
        self.S = SourcePopulation_(D=self.D, population=sourcepop, cosmo=cosmo)

        # Initialize Einstein radius calculator
        self.E = EinsteinRadiusTools(D=self.D)

        # Get survey depth parameters
        self._set_survey_params()

        print(f"\nSetup complete!")
        print("="*70)

    def _set_survey_params(self):
        """Set survey-specific selection parameters"""
        # These are approximate depth limits for different surveys
        # Used for quick magnitude-based filtering
        survey_params = {
            'Euclid': {'mag_limit': 24.5, 'area_deg2': 14000},
            'LSST': {'mag_limit': 27.0, 'area_deg2': 18000},
            'DES': {'mag_limit': 24.0, 'area_deg2': 5000},
            'COSMOS-Web': {'mag_limit': 27.5, 'area_deg2': 0.54},
        }

        if self.survey in survey_params:
            self.mag_limit = survey_params[self.survey]['mag_limit']
            self.area = survey_params[self.survey]['area_deg2']
        else:
            print(f"Warning: Unknown survey {self.survey}, using Euclid defaults")
            self.mag_limit = 24.5
            self.area = 14000

        print(f"  Survey area: {self.area} deg²")
        print(f"  Approximate magnitude limit: {self.mag_limit}")

    def generate_forecast(self, N_target_lenses=15999,
                          oversample_factor=10,
                          mag_cut=3.0,
                          snr_proxy_cut=5.0,
                          sourceplaneoverdensity=1):
        """
        Generate lens forecast with vectorized operations

        Parameters:
        -----------
        N_target_lenses : int
            Target number of detected lenses
        oversample_factor : float
            Oversampling factor (generate more pairs than needed)
        mag_cut : float
            Magnification threshold
        snr_proxy_cut : float
            Approximate SNR threshold (simplified)
        sourceplaneoverdensity : float
            Source plane overdensity factor

        Returns:
        --------
        dict : Forecast results including redshift distributions
        """
        t_start = time.perf_counter()

        print(f"\nGenerating forecast for {N_target_lenses} lenses...")
        print(f"Oversampling factor: {oversample_factor}x")

        # Estimate how many lens-source pairs to generate
        # Typical lensing probability is ~0.1-1%, so we need lots of pairs
        N_pairs = int(N_target_lenses * oversample_factor * 100)

        print(f"\nStep 1: Drawing {N_pairs:,} lens-source pairs...")
        t0 = time.perf_counter()

        # Draw lens population (vectorized)
        zl, sigl, ml, rl, ql = self.L.drawLensPopulation(N_pairs)

        # Draw source population (vectorized)
        zs, ms, xs, ys, qs, ps, rs = self.S.drawSourcePopulation(
            N_pairs,
            sourceplaneoverdensity=sourceplaneoverdensity,
            returnmasses=False
        )

        t1 = time.perf_counter()
        print(f"  Generated in {t1-t0:.2f}s")

        # Step 2: Geometric lensing criterion (vectorized)
        print(f"\nStep 2: Applying geometric lensing criterion...")
        t0 = time.perf_counter()

        # Calculate Einstein radii (vectorized)
        b = self.E.sie_rein(sigl, zl, zs)

        # Geometric criterion: b² > x² + y²
        r_src_sq = xs**2 + ys**2
        b_sq = b**2
        is_lensed = b_sq > r_src_sq

        N_geometrically_lensed = np.sum(is_lensed)
        t1 = time.perf_counter()
        print(f"  {N_geometrically_lensed:,} / {N_pairs:,} satisfy geometric criterion ({100*N_geometrically_lensed/N_pairs:.2f}%)")
        print(f"  Applied in {t1-t0:.2f}s")

        if N_geometrically_lensed == 0:
            raise RuntimeError("No lenses satisfy geometric criterion! Increase N_pairs or adjust parameters.")

        # Step 3: Apply selection cuts (vectorized)
        print(f"\nStep 3: Applying selection cuts...")
        t0 = time.perf_counter()

        # Filter to geometrically lensed systems only
        zl_lensed = zl[is_lensed]
        zs_lensed = zs[is_lensed]
        sigl_lensed = sigl[is_lensed]
        b_lensed = b[is_lensed]
        xs_lensed = xs[is_lensed]
        ys_lensed = ys[is_lensed]
        rs_lensed = rs[is_lensed]

        # Get source magnitudes (use first available band)
        band = list(ms.keys())[0]
        ms_lensed = ms[band][is_lensed]

        # Approximate magnification (SIS approximation)
        # μ ≈ b / sqrt(x² + y²) for SIS
        r_src = np.sqrt(xs_lensed**2 + ys_lensed**2)
        r_src[r_src < 0.01] = 0.01  # Avoid division by zero
        mu_approx = b_lensed / r_src

        # Magnification cut
        mag_cut_mask = mu_approx > mag_cut

        # Magnitude cut (source must be detectable when magnified)
        # Approximate: magnified source brightness
        ms_magnified = ms_lensed - 2.5*np.log10(mu_approx)
        mag_limit_mask = ms_magnified < self.mag_limit

        # SNR proxy: brighter sources and larger sources easier to detect
        # This is a simplified proxy without full ray-tracing
        # Assume SNR ∝ flux / sqrt(area) ∝ 10^(-0.4*mag) / r_src
        snr_proxy = 10**(-0.4 * ms_magnified) / (rs_lensed + 0.1)
        snr_proxy_mask = snr_proxy > (snr_proxy_cut * np.median(snr_proxy))

        # Combined selection
        detected = mag_cut_mask & mag_limit_mask & snr_proxy_mask

        N_detected = np.sum(detected)
        t1 = time.perf_counter()
        print(f"  Magnification > {mag_cut}: {np.sum(mag_cut_mask):,}")
        print(f"  Magnitude limit < {self.mag_limit}: {np.sum(mag_limit_mask):,}")
        print(f"  SNR proxy cut: {np.sum(snr_proxy_mask):,}")
        print(f"  Combined: {N_detected:,} detected lenses")
        print(f"  Applied in {t1-t0:.2f}s")

        if N_detected == 0:
            raise RuntimeError("No lenses pass selection cuts! Relax cuts or increase oversampling.")

        # Extract final sample
        if N_detected > N_target_lenses:
            # Randomly subsample to target
            detected_indices = np.where(detected)[0]
            selected_indices = np.random.choice(detected_indices, N_target_lenses, replace=False)
            final_mask = np.zeros(len(detected), dtype=bool)
            final_mask[selected_indices] = True
        else:
            final_mask = detected
            print(f"  Warning: Only found {N_detected} lenses, less than target {N_target_lenses}")

        # Final redshifts
        zl_final = zl_lensed[final_mask]
        zs_final = zs_lensed[final_mask]
        sigl_final = sigl_lensed[final_mask]
        b_final = b_lensed[final_mask]
        mu_final = mu_approx[final_mask]

        t_end = time.perf_counter()

        # Summary
        print(f"\n" + "="*70)
        print("FORECAST COMPLETE")
        print("="*70)
        print(f"Generated {len(zl_final):,} detected lenses in {t_end-t_start:.2f} seconds")
        print(f"\nLens redshifts:")
        print(f"  Range: {zl_final.min():.3f} - {zl_final.max():.3f}")
        print(f"  Median: {np.median(zl_final):.3f}")
        print(f"  Mean: {np.mean(zl_final):.3f}")
        print(f"\nSource redshifts:")
        print(f"  Range: {zs_final.min():.3f} - {zs_final.max():.3f}")
        print(f"  Median: {np.median(zs_final):.3f}")
        print(f"  Mean: {np.mean(zs_final):.3f}")
        print(f"\nEinstein radii:")
        print(f"  Range: {b_final.min():.3f} - {b_final.max():.3f} arcsec")
        print(f"  Median: {np.median(b_final):.3f} arcsec")
        print(f"\nMagnifications:")
        print(f"  Range: {mu_final.min():.1f} - {mu_final.max():.1f}")
        print(f"  Median: {np.median(mu_final):.1f}")
        print("="*70)

        # Package results
        results = {
            'survey': self.survey,
            'cosmo': self.cosmo,
            'sourcepop': self.sourcepop,
            'N_lenses': len(zl_final),
            'N_target': N_target_lenses,
            'runtime_seconds': t_end - t_start,
            'zl': zl_final,
            'zs': zs_final,
            'sigl': sigl_final,
            'b': b_final,
            'mu_approx': mu_final,
            'mag_cut': mag_cut,
            'mag_limit': self.mag_limit,
        }

        return results

    def plot_results(self, results, output_file='forecast_redshifts.png'):
        """
        Plot redshift distributions
        """
        try:
            plt.style.use('sanglier')
        except:
            print("Warning: 'sanglier' style not found, using default")

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Lens redshifts
        ax = axes[0, 0]
        ax.hist(results['zl'], bins=30, alpha=0.7, edgecolor='black')
        ax.set_xlabel('Lens Redshift $z_l$', fontsize=12)
        ax.set_ylabel('Number of Lenses', fontsize=12)
        ax.set_title(f'Lens Redshift Distribution (N={results["N_lenses"]:,})', fontsize=14)
        ax.axvline(np.median(results['zl']), color='red', linestyle='--',
                   label=f'Median = {np.median(results["zl"]):.3f}')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Source redshifts
        ax = axes[0, 1]
        ax.hist(results['zs'], bins=30, alpha=0.7, edgecolor='black', color='orange')
        ax.set_xlabel('Source Redshift $z_s$', fontsize=12)
        ax.set_ylabel('Number of Sources', fontsize=12)
        ax.set_title(f'Source Redshift Distribution', fontsize=14)
        ax.axvline(np.median(results['zs']), color='red', linestyle='--',
                   label=f'Median = {np.median(results["zs"]):.3f}')
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Lens-source redshift plane
        ax = axes[1, 0]
        h = ax.hist2d(results['zl'], results['zs'], bins=30, cmap='viridis')
        plt.colorbar(h[3], ax=ax, label='Number of lenses')
        ax.set_xlabel('Lens Redshift $z_l$', fontsize=12)
        ax.set_ylabel('Source Redshift $z_s$', fontsize=12)
        ax.set_title('Lens-Source Redshift Plane', fontsize=14)
        # Add diagonal (zs = zl) - physically impossible
        zmax = max(results['zl'].max(), results['zs'].max())
        ax.plot([0, zmax], [0, zmax], 'r--', alpha=0.5, label='$z_s = z_l$')
        ax.legend()

        # Einstein radius distribution
        ax = axes[1, 1]
        ax.hist(results['b'], bins=30, alpha=0.7, edgecolor='black', color='green')
        ax.set_xlabel('Einstein Radius (arcsec)', fontsize=12)
        ax.set_ylabel('Number of Lenses', fontsize=12)
        ax.set_title('Einstein Radius Distribution', fontsize=14)
        ax.axvline(np.median(results['b']), color='red', linestyle='--',
                   label=f'Median = {np.median(results["b"]):.3f}"')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.suptitle(f'{results["survey"]} Strong Lens Forecast\n' +
                     f'Source catalog: {results["sourcepop"]} | ' +
                     f'Runtime: {results["runtime_seconds"]:.1f}s',
                     fontsize=16, y=0.995)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved: {output_file}")

        return fig


def main():
    """Command-line interface"""
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    survey = sys.argv[1]
    N_target = int(sys.argv[2])
    output_file = sys.argv[3] if len(sys.argv) > 3 else f'{survey.lower()}_forecast.pkl'

    # Initialize forecast
    forecast = FastLensForecast(survey=survey, sourcepop='lsst')

    # Generate forecast
    results = forecast.generate_forecast(N_target_lenses=N_target)

    # Save results
    with open(output_file, 'wb') as f:
        pickle.dump(results, f, protocol=5)
    print(f"\nResults saved: {output_file}")

    # Plot
    plot_file = output_file.replace('.pkl', '.png')
    forecast.plot_results(results, plot_file)

    print("\nDone!")


if __name__ == '__main__':
    main()
