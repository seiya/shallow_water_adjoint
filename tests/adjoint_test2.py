import numpy as np
import subprocess
from pathlib import Path


def save_field(path, arr):
    path = Path(path)
    arr.ravel(order='F').astype(np.float64).tofile(path)


def read_snapshot(path, nx, ny):
    path = Path(path)
    data = np.fromfile(path, dtype=np.float32, count=nx*ny)
    return data.reshape((nx, ny), order='F').astype(np.float64)


def geostrophic_height(nx, ny):
    pi = np.pi
    radius = 6371220.0
    g = 9.80616
    day = 86400.0
    omega = 2.0 * pi / (12.0 * day)
    h0 = 10000.0
    u0 = 20.0
    dlat = pi / ny
    lat = -pi/2 + (np.arange(ny) + 0.5) * dlat
    coeff = radius * omega * u0 / g
    hlat = h0 + coeff * np.sin(lat) ** 2
    return np.repeat(hlat[np.newaxis, :], nx, axis=0)


def main():
    build_dir = Path(__file__).resolve().parents[1] / 'build'
    exe_fwd = build_dir / 'shallow_water_test2_forward.out'
    exe_rev = build_dir / 'shallow_water_test2_reverse.out'

    nx, ny = 128, 64
    rng = np.random.default_rng(0)
    x = geostrophic_height(nx, ny)
    u = rng.standard_normal((nx, ny))
    v = rng.standard_normal()

    x_file = build_dir / 'x2.bin'
    u_file = build_dir / 'u2.bin'
    save_field(x_file, x)
    save_field(u_file, u)

    res = subprocess.run(
        [str(exe_fwd), '-1', str(x_file), str(u_file)],
        check=True, cwd=build_dir, capture_output=True, text=True
    )
    Ju = float(res.stdout.strip().split()[0])

    subprocess.run(
        [str(exe_rev), '0', str(x_file)],
        check=True, cwd=build_dir, capture_output=True, text=True
    )
    g = read_snapshot(build_dir / 'snapshot_0000.bin', nx, ny)

    JT_v = v * g
    lhs = v * Ju
    rhs = np.vdot(u.ravel(order='F'), JT_v.ravel(order='F'))
    diff = abs(lhs - rhs)
    tol = max(1e-12, 1e-6 * max(abs(lhs), abs(rhs)))
    print(f"vTJu={lhs:.6e} uTJTv={rhs:.6e} diff={diff:.6e}")
    assert diff < tol


if __name__ == '__main__':
    main()
