# Example Data

This directory is intended for example GROMACS trajectory files used for
testing and demonstration purposes.

## Files expected

| File          | Description                                           |
|---------------|-------------------------------------------------------|
| `sample.tpr`  | GROMACS binary run input (topology) file              |
| `traj.xtc`    | Trajectory for analysis                               |
| `broken.xtc`  | Trajectory with intentional PBC/drift issues (demo)  |
| `fixed.xtc`   | The same trajectory after `gmx trjconv` preprocessing |

## Generating test data

If you do not have GROMACS data handy, you can generate a minimal test
system using the GROMACS built-in alanine dipeptide tutorial:

```bash
# Download the GROMACS benchmark dataset
wget https://www.gromacs.org/benchmark/gromacs_benchmark.tar.gz
tar xvf gromacs_benchmark.tar.gz
```

Or use the MDAnalysis test files available at:
```python
import MDAnalysis.tests.datafiles as mda_data
print(mda_data.TPR)   # path to a bundled .tpr file
print(mda_data.XTC)   # path to a bundled .xtc file
```

See the repository README for full usage instructions.
