import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button
import struct
import sys
import os

def read_movie_gal(filename, n_particles):
    """Reads the binary movie.gal file and returns a list of frames."""
    frames = []
    frame_size = n_particles * 3 * 8 
    
    if not os.path.exists(filename):
        print(f"Error: {filename} not found.")
        return None

    with open(filename, 'rb') as f:
        while True:
            data = f.read(frame_size)
            if not data or len(data) < frame_size:
                break
            unpacked = struct.unpack('d' * (n_particles * 3), data)
            frame = np.array(unpacked).reshape((n_particles, 3))
            frames.append(frame)
    return frames

def int_to_rgb(val):
    """Converts a float containing packed RGB back to normalized RGB tuple."""
    ival = int(val)
    r = ((ival >> 16) & 0xFF) / 255.0
    g = ((ival >> 8) & 0xFF) / 255.0
    b = (ival & 0xFF) / 255.0
    return (r, g, b)

def animate_simulation(frames, n_particles, save_path=None, speedup=1, fps=30):
    if not frames:
        print("No frames to animate.")
        return

    # Apply speedup
    if speedup > 1:
        frames = frames[::speedup]
        print(f"Speedup applied: Using every {speedup}th frame ({len(frames)} total).")

    num_frames = len(frames)
    
    # Setup Figure
    fig, ax = plt.subplots(figsize=(10, 10))
    if save_path is None:
        plt.subplots_adjust(bottom=0.15)
    
    # Ultra-dark cinematic theme
    fig.patch.set_facecolor('#050505')
    ax.set_facecolor('#050505')
    ax.axis('off') # Hide axes for a cleaner look

    # Calculate speeds for color mapping
    # Note: movie.gal only stores x, y, mass. We'll approximate speed from displacement
    speeds = []
    for i in range(num_frames):
        if i == 0:
            s = np.zeros(n_particles)
        else:
            # dx^2 + dy^2
            dist = np.sqrt(np.sum((frames[i][:, :2] - frames[i-1][:, :2])**2, axis=1))
            s = dist
        speeds.append(s)
    
    # Normalize speeds for color mapping
    all_speeds = np.concatenate(speeds)
    max_speed = np.percentile(all_speeds, 99) if len(all_speeds) > 0 else 1.0
    if max_speed == 0: max_speed = 1.0

    # Initial plot with "glow" effect
    # We use two scatter layers: a core and a glow
    glow = ax.scatter(frames[0][:, 0], frames[0][:, 1], s=12, c='#5e81ac', alpha=0.15, edgecolors='none')
    core = ax.scatter(frames[0][:, 0], frames[0][:, 1], s=1.5, c='#88c0d0', alpha=0.9, edgecolors='none')
    
    # Smart Zoom: Focus on the center of the action
    def get_frame_limits(frame, padding=1.2):
        # Use percentiles to ignore outliers (particles ejected at high speed)
        x_min, x_max = np.percentile(frame[:, 0], [2, 98])
        y_min, y_max = np.percentile(frame[:, 1], [2, 98])
        
        center_x, center_y = (x_min + x_max)/2, (y_min + y_max)/2
        range_val = max(x_max - x_min, y_max - y_min) * padding
        if range_val < 0.1: range_val = 0.1 # Prevent zero division
        
        return [center_x - range_val/2, center_x + range_val/2], \
               [center_y - range_val/2, center_y + range_val/2]

    # Pre-calculate limits for all frames to smooth them
    all_limits = [get_frame_limits(f) for f in frames]
    
    # Smooth limits using a moving average window
    window = 10
    smoothed_xlim = []
    smoothed_ylim = []
    for i in range(num_frames):
        start = max(0, i - window)
        end = min(num_frames, i + window)
        avg_x = np.mean([all_limits[j][0] for j in range(start, end)], axis=0)
        avg_y = np.mean([all_limits[j][1] for j in range(start, end)], axis=0)
        smoothed_xlim.append(avg_x)
        smoothed_ylim.append(avg_y)

    ax.set_xlim(smoothed_xlim[0])
    ax.set_ylim(smoothed_ylim[0])
    ax.set_aspect('equal')
    
    # Cinematic title
    title = ax.text(0.5, 0.95, f"N-Body Cinematic Simulation (N={n_particles})", 
                    transform=ax.transAxes, color='#d8dee9', ha='center', fontsize=14, fontweight='light')

    def update_scatter(i):
        frame = frames[i]
        speed = speeds[i]
        norm_speed = np.clip(speed / max_speed, 0, 1)
        
        # Update positions
        core.set_offsets(frame[:, :2])
        glow.set_offsets(frame[:, :2])
        
        # Update limits (Tracking Camera)
        ax.set_xlim(smoothed_xlim[i])
        ax.set_ylim(smoothed_ylim[i])
        
        # Dynamic color mapping
        colors = plt.cm.magma(norm_speed)
        core.set_color(colors)
        glow.set_color(colors)
        return core, glow

    if save_path:
        print(f"Exporting cinematic animation to {save_path}...")
        def animate_frame(i):
            if i % 20 == 0:
                print(f"Processing frame {i}/{num_frames}...")
            return update_scatter(i)
        
        ani = animation.FuncAnimation(fig, animate_frame, frames=num_frames, blit=True)
        
        if save_path.endswith('.gif'):
            ani.save(save_path, writer='pillow', fps=fps, savefig_kwargs={'facecolor': '#050505'})
        else:
            ani.save(save_path, writer='ffmpeg', fps=fps, savefig_kwargs={'facecolor': '#050505'})
        print(f"Export complete: {save_path}")
        plt.close()
        return

    # Interactive mode
    ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03], facecolor='#2e3440')
    slider = Slider(ax_slider, 'Time', 0, num_frames - 1, valinit=0, valstep=1, color='#81a1c1')
    slider.label.set_color('#d8dee9')

    def update_slider(val):
        update_scatter(int(slider.val))
        fig.canvas.draw_idle()

    slider.on_changed(update_slider)

    is_paused = False
    def animate(i):
        if not is_paused:
            new_val = (slider.val + 1) % num_frames
            slider.set_val(new_val)
        return core, glow

    ani = animation.FuncAnimation(fig, animate, interval=30, blit=True, cache_frame_data=False)
    plt.show()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Visualize N-Body simulation results.")
    parser.add_argument("filename", help="Path to movie.gal")
    parser.add_argument("n_particles", type=int, help="Number of particles")
    parser.add_argument("--save", help="Save animation to file (e.g., output.gif or output.mp4)")
    parser.add_argument("--speedup", type=int, default=1, help="Speedup factor (skip frames)")
    parser.add_argument("--fps", type=int, default=30, help="Frames per second for export")
    
    args = parser.parse_args()

    print(f"Reading {args.filename}...")
    frames = read_movie_gal(args.filename, args.n_particles)
    if frames:
        animate_simulation(frames, args.n_particles, save_path=args.save, speedup=args.speedup, fps=args.fps)
