import struct
import random
import sys
import math
import numpy as np

def _sample_truncated_exponential_radius(r_min, r_max, r_scale):
    """
    Sample r in [r_min, r_max] from an exponential disk profile:
      p(r) ∝ exp(-(r-r_min)/r_scale)

    Uses inverse CDF of the truncated distribution for speed (no rejection).
    """
    if r_max <= 0:
        return 0.0
    if r_min < 0:
        r_min = 0.0
    if r_max <= r_min:
        return float(r_min)
    if r_scale <= 0:
        # Degenerate: everything at the center.
        return float(r_min)
    u = random.random()
    span = r_max - r_min
    z = 1.0 - math.exp(-span / r_scale)
    # t = -r_scale * ln(1 - u * z), r = r_min + t
    return float(r_min) + (-r_scale * math.log(1.0 - u * z))


def generate_disk(
    N,
    filename,
    center_x=0.0,
    center_y=0.0,
    radius=10.0,
    drift_vx=0.0,
    drift_vy=0.0,
    *,
    g_sim=1.0,
    central_mass=100.0,
    r_scale=None,
    r_min=0.05,
    r_soft=0.05,
    mass_min=1e-3,
    mass_max=1e-2,
    v_jitter=0.01,
):
    """
    Stable benchmark disk: 1 heavy central body + (N-1) light orbiters.

    - Positions follow a center-dense exponential disk (truncated to `radius`).
    - Velocities approximate circular orbits around the central mass.
    - `g_sim` is constant (does NOT depend on N).
    """
    if N <= 0:
        raise ValueError("N must be positive")

    if r_scale is None:
        # Heuristic: visible center density without extreme clustering.
        r_scale = radius / 5.0 if radius > 0 else 1.0

    print(
        f"Generating {N} particles in a center-dense DISK at ({center_x}, {center_y}) "
        f"R=[{r_min}, {radius}] (r_scale={r_scale}) with drift ({drift_vx}, {drift_vy}) to {filename}..."
    )
    with open(filename, 'wb') as f:
        # Particle 0: central mass at the origin of the disk (dominant potential).
        brightness = random.uniform(0.5, 1.5)
        data = struct.pack('6d', center_x, center_y, float(central_mass), drift_vx, drift_vy, brightness)
        f.write(data)

        # Remaining particles: light orbiters.
        for _ in range(N - 1):
            theta = random.uniform(0.0, 2.0 * math.pi)
            r = _sample_truncated_exponential_radius(r_min, radius, r_scale)

            pos_x = center_x + r * math.cos(theta)
            pos_y = center_y + r * math.sin(theta)

            # Circular orbit speed around the central mass (softened near r=0).
            r_eff = r + r_soft
            v_speed = math.sqrt(g_sim * central_mass / r_eff) if r_eff > 0 else 0.0

            # Tangential velocity (perpendicular to radial vector).
            # Using theta avoids an extra normalization step.
            vx = -v_speed * math.sin(theta)
            vy = v_speed * math.cos(theta)

            # Small noise keeps the disk from being perfectly uniform.
            if v_jitter and v_jitter > 0:
                vx *= 1.0 + random.uniform(-v_jitter, v_jitter)
                vy *= 1.0 + random.uniform(-v_jitter, v_jitter)

            vx += drift_vx
            vy += drift_vy

            mass = random.uniform(mass_min, mass_max)
            brightness = random.uniform(0.5, 1.5)
            data = struct.pack('6d', pos_x, pos_y, mass, vx, vy, brightness)
            f.write(data)
    print("Done.")

def generate_random(N, filename):
    print(f"Generating {N} random particles to {filename}...")
    with open(filename, 'wb') as f:
        for _ in range(N):
            pos_x = random.uniform(-1, 1)
            pos_y = random.uniform(-1, 1)
            mass = random.uniform(0.1, 10.0)
            vx = random.uniform(-0.05, 0.05)
            vy = random.uniform(-0.05, 0.05)
            brightness = random.uniform(0.5, 1.5)
            data = struct.pack('6d', pos_x, pos_y, mass, vx, vy, brightness)
            f.write(data)
    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 generate_data.py <N> <filename> [mode: random|disk] [center_x center_y radius drift_vx drift_vy]")
        sys.exit(1)

    N = int(sys.argv[1])
    filename = sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "disk"

    if mode == "disk":
        center_x = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
        center_y = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
        radius = float(sys.argv[6]) if len(sys.argv) > 6 else 10.0
        drift_vx = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
        drift_vy = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0
        generate_disk(N, filename, center_x=center_x, center_y=center_y, radius=radius, drift_vx=drift_vx, drift_vy=drift_vy)
    elif mode == "random":
        generate_random(N, filename)
    else:
        print(f"Unknown mode: {mode}")
        print("Usage:")
        print("  Random: python3 generate_data.py <N> <filename> random")
        print("  Disk:   python3 generate_data.py <N> <filename> disk [center_x center_y radius drift_vx drift_vy]")
        sys.exit(1)
