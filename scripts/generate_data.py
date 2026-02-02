import struct
import random
import sys
import math
import numpy as np

def generate_disk(N, filename, center_x=0, center_y=0, radius=0.5, mass_val=1.0, drift_vx=0.0, drift_vy=0.0):
    print(f"Generating {N} particles in a DISK at ({center_x}, {center_y}) with drift ({drift_vx}, {drift_vy}) to {filename}...")
    with open(filename, 'wb') as f:
        for i in range(N):
            # Distance from center (exponential or power law distribution for density)
            r = radius * math.sqrt(random.random())
            theta = random.uniform(0, 2 * math.pi)
            
            pos_x = center_x + r * math.cos(theta)
            pos_y = center_y + r * math.sin(theta)
            
            # Orbital velocity: v = sqrt(G * M / r)
            # Heuristic G=100/N as used in the C code
            G = 100.0 / N
            # Simple orbital speed heuristic
            v_speed = math.sqrt(G * (N/2) / (r + 0.05)) * 0.4 
            
            vx = -v_speed * math.sin(theta) + drift_vx
            vy = v_speed * math.cos(theta) + drift_vy
            
            mass = random.uniform(0.1, 2.0)
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
        print("Usage: python3 generate_data.py <N> <filename> [mode: random|disk] [center_x center_y radius]")
        sys.exit(1)

    N = int(sys.argv[1])
    filename = sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "disk"

    if mode == "disk":
        center_x = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
        center_y = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
        radius = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5
        drift_vx = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
        drift_vy = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0
        generate_disk(N, filename, center_x=center_x, center_y=center_y, radius=radius, drift_vx=drift_vx, drift_vy=drift_vy)
    elif mode == "random":
        generate_random(N, filename)
    else:
        print(f"Unknown mode: {mode}")
        print("Usage:")
        print("  Random: python3 generate_data.py <N> <filename> random")
        print("  Disk:   python3 generate_data.py <N> <filename> disk [center_x center_y radius]")
        sys.exit(1)
