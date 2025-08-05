import numpy as np
import subprocess
from pathlib import Path

def save_field(path, arr):
    path = Path(path)
    arr.ravel(order='F').astype(np.float64).tofile(path)

def read_snapshot(path, nlon, nlat):
    path = Path(path)
    data = np.fromfile(path, dtype=np.float32, count=nlon*nlat)
    return data.reshape((nlon, nlat), order='F').astype(np.float64)

def main():
    build_dir = Path(__file__).resolve().parents[1] / 'build'
    exe_fwd = build_dir / 'shallow_water_test1_forward.out'
    exe_rev = build_dir / 'shallow_water_test1_reverse.out'

    nlon, nlat = 128, 64
    rng = np.random.default_rng(0)
    x = rng.standard_normal((nlon, nlat))
    u = rng.standard_normal((nlon, nlat))
    v = rng.standard_normal()

    x_file = build_dir / 'x.bin'
    u_file = build_dir / 'u.bin'
    save_field(x_file, x)
    save_field(u_file, u)

    res = subprocess.run(
        [str(exe_fwd), '0', '0', str(x_file), str(u_file)],
        check=True, cwd=build_dir, capture_output=True, text=True
    )
    Ju = float(res.stdout.strip().split()[0])

    subprocess.run(
        [str(exe_rev), '0', '1', str(x_file)],
        check=True, cwd=build_dir, capture_output=True, text=True
    )
    g = read_snapshot(build_dir / 'snapshot_0000.bin', nlon, nlat)

    JT_v = v * g
    lhs = v * Ju
    rhs = np.vdot(u.ravel(order='F'), JT_v.ravel(order='F'))
    diff = abs(lhs - rhs)
    tol = 1e-6 * max(abs(lhs), abs(rhs))
    print(f"vTJu={lhs:.6e} uTJTv={rhs:.6e} diff={diff:.6e}")
    assert diff < tol

if __name__ == '__main__':
    main()
