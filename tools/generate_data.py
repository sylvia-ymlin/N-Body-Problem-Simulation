import argparse
import numpy as np
import struct

def generate_random(N):
    pos_x = np.random.rand(N)
    pos_y = np.random.rand(N)
    mass = np.ones(N) * (1.0 / N)
    vx = np.random.rand(N) * 0.1
    vy = np.random.rand(N) * 0.1
    brightness = np.random.rand(N)
    return pos_x, pos_y, mass, vx, vy, brightness

def generate_cluster(N):
    # One main dense cluster
    N_cluster = int(N * 0.8)
    N_bg = N - N_cluster
    
    # Cluster center (0.5, 0.5), small spread
    pos_x_c = np.random.normal(0.5, 0.05, N_cluster)
    pos_y_c = np.random.normal(0.5, 0.05, N_cluster)
    
    # Background uniform
    pos_x_bg = np.random.rand(N_bg)
    pos_y_bg = np.random.rand(N_bg)
    
    pos_x = np.concatenate([pos_x_c, pos_x_bg])
    pos_y = np.concatenate([pos_y_c, pos_y_bg])
    
    # Clip to [0,1]
    pos_x = np.clip(pos_x, 0.01, 0.99)
    pos_y = np.clip(pos_y, 0.01, 0.99)
    
    mass = np.ones(N) * (1.0 / N)
    vx = np.zeros(N)
    vy = np.zeros(N)
    brightness = np.ones(N)
    
    return pos_x, pos_y, mass, vx, vy, brightness

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('N', type=int)
    parser.add_argument('output')
    parser.add_argument('--mode', default='random', choices=['random', 'cluster'])
    args = parser.parse_args()
    
    if args.mode == 'random':
        pos_x, pos_y, mass, vx, vy, bright = generate_random(args.N)
    else:
        pos_x, pos_y, mass, vx, vy, bright = generate_cluster(args.N)
        
    with open(args.output, 'wb') as f:
        for i in range(args.N):
            data = struct.pack('dddddd', pos_x[i], pos_y[i], mass[i], vx[i], vy[i], bright[i])
            f.write(data)

if __name__ == '__main__':
    main()
