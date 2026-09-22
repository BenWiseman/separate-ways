# External data for the geometry-only fits

The optional Pantheon+ distance table and covariance are acquired separately.
The publisher's [data repository](https://github.com/PantheonPlusSH0ES/DataRelease)
identifies the products under `Pantheon+_Data/4_DISTANCES_AND_COVAR/`; its
[format description](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/README)
defines the columns and covariance ordering. Cite the Pantheon+ analysis
[arXiv:2202.04077](https://arxiv.org/abs/2202.04077) when using these inputs.

The inputs used by the retained calculations match these upstream files byte for
byte at commit `c447f0fea703fcd0fff57de5000947b5ca81286b` (checked 18 September
2026). The release uses `Plus` in local filenames where the upstream names use
`+`; the contents are unchanged.

| Local filename | Bytes | SHA256 |
|---|---:|---|
| `PantheonPlusSH0ES.dat` | 579283 | `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8` |
| `PantheonPlusSH0ES_STAT+SYS.cov` | 33284960 | `abf806d966485e64afdb359c87bffc0ecc00d05eff0a31ced66f247385df0fdc` |

From the release root, acquire and verify the two files:

```bash
mkdir -p data/pantheon_plus
curl --fail --location 'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat' -o data/pantheon_plus/PantheonPlusSH0ES.dat
curl --fail --location 'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov' -o 'data/pantheon_plus/PantheonPlusSH0ES_STAT+SYS.cov'
python3 data/check_inputs.py
```

An existing copy can stay outside the release. Set `PAPER2_PANTHEON_DIR` to its
directory; both fit scripts and `check_inputs.py` respect that setting. The
default is this release's `data/pantheon_plus`, resolved from the script location.

The upstream repository did not expose a license through the GitHub API, and no
license file was found in its recursive tree during this check. This release
therefore supplies acquisition instructions and hashes, without redistributing
the external data or asserting new distribution rights over them.

Both fits select `zHD > 0.01` and `IS_CALIBRATOR == 0`: 1580 of the 1701 input
rows. They use `m_b_corr` and the full selected STAT+SYS covariance, with the
supernova magnitude offset marginalised analytically. They do not use the
diagonal plotting uncertainties as a likelihood covariance.

The DESI DR2 BAO table and block covariance builder are in
`calc/desi_dr2_geometry.py`, with source [arXiv:2503.14738](https://arxiv.org/abs/2503.14738).
They yield a 13-element distance vector. No external DESI file, CAMB installation,
or Planck likelihood package is needed for these particular geometry-only fits.
The compressed CMB priors are numerical assumptions in the fit sources, not a
substitute for the complete Planck likelihood.
