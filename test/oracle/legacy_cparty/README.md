# External CParty inside oracle

The semantic reference is upstream CParty inside at a fixed commit plus the
ordered local patches in [config.json](config.json):

```text
upstream  https://github.com/TheCOBRALab/CParty.git
commit    800e95af9d1a97555233e2c7b37d73dbaa1a92b4
```

Compare with explicit `WmbClosurePolicy::LegacyCompatible`. The production
LinearKnot default requires completed-WMB weak closure; the official oracle
has no completed-WMB closure patch. Legacy MFE, BPP and sampling are not
inside semantic oracles.

## Build

From the repository root:

```bash
test/oracle/legacy_cparty/build.sh
test/oracle/legacy_cparty/build.sh /path/to/output
```

Building the optional oracle requires network access. The script clones
upstream, resets its disposable output checkout to the pin, applies the
explicit patch list, and builds with -j8. Rebuilding replaces local changes
in that generated checkout; use a dedicated output directory.
The default output is `build/legacy_cparty_oracle/`.

The `cparty-oracle` wrapper invokes CParty with sampling disabled. The
manifest records the pin, ordered patch hashes, patched source and binary
identity. Probe helpers check that manifest before linking.

## Patch inventory

| Patch | Scope |
| --- | --- |
| 0001 | Include algorithm in centroid.cc for std::sort; build only |
| 0002 | Read the input-structure option before parser cleanup; CLI only |
| 0006 | VPR right-run guard includes the full unpaired run through j |
| 0007 | Explicit PF extrapolation and capped asymmetry for scaffold loops larger than the parameter table |
| 0008 | Include the minimal-hairpin endpoint in WIP forward and traceback scans |
| 0009 | Supply the two inner neighbors to the multiloop PF closing helper at sequence endpoints |

Former patches 0003 (WMB endpoint), 0004 (VP internal-loop border), and 0005
(VPR exponent) are incorporated in the pinned upstream source. Local patch
numbers are not upstream issue numbers. The patch files retain original
CParty context and attribution; see the [third-party notices](../../../THIRD_PARTY_NOTICES.md).

Frozen inside and local-weight witnesses run without this download. For
renewed random comparisons, use `test/tools/audit_legacy_inside_parity.py`.
See [validation](../../../documents/validation.md) and
[model policy](../../../documents/model.md).
