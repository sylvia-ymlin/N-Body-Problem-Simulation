import struct
import random
import sys
import os

def generate_file(filename, N):
    print(f"Generating {N} particles to {filename}...")
    with open(filename, 'wb') as f:
        for _ in range(N):
            # pos_x, pos_y in range [-1000, 1000]
            pos_x = random.uniform(0, 1)
            pos_y = random.uniform(0, 1)
            mass = 1.0
            vx = 0.0
            vy = 0.0
            brightness = 1.0
            
            data = struct.pack('6d', pos_x, pos_y, mass, vx, vy, brightness)
            f.write(data)
    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 generate_data.py <N> <filename>")
        sys.exit(1)
    
    N = int(sys.argv[1])
    filename = sys.argv[2]
    generate_file(filename, N)
