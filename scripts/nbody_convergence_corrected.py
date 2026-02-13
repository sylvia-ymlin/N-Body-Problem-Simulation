#!/usr/bin/env python3
"""
N-Body Convergence Order Validation (Corrected)
在实际的N-body问题上验证Symplectic Euler、Velocity Verlet和RK4的收敛阶数
只计算position误差，reference line与对应方法对齐
"""

import numpy as np
import matplotlib.pyplot as plt
import time

# Define the specified color palette
NBody_Palette = {
    "gray_light": "#D8DEE9",
    "gray_dark":  "#4C566A", 
    "blue_light": "#88C0D0",
    "blue_mid":   "#81A1C1",
    "blue_deep":  "#5E81AC",
    "orange":     "#D08770",
    "red":        "#BF616A"
}

# Gravitational constant
G = 1.0

def compute_force(pos, mass):
    """Compute gravitational forces between all particles"""
    n = len(pos)
    force = np.zeros_like(pos)
    
    for i in range(n):
        for j in range(i + 1, n):
            r = pos[j] - pos[i]
            r_mag = np.linalg.norm(r)
            
            # Softening to avoid singularities
            softening = 0.01
            r_mag_soft = np.sqrt(r_mag**2 + softening**2)
            
            force_ij = G * mass[i] * mass[j] * r / (r_mag_soft**3)
            force[i] += force_ij
            force[j] -= force_ij
    
    return force

# Integration methods
def symplectic_euler_step(pos, vel, mass, dt):
    """Symplectic Euler (1st order)"""
    force = compute_force(pos, mass)
    vel += force / mass[:, None] * dt
    pos += vel * dt
    return pos, vel

def velocity_verlet_step(pos, vel, mass, dt):
    """Velocity Verlet (2nd order)"""
    force = compute_force(pos, mass)
    pos += vel * dt + 0.5 * force / mass[:, None] * dt**2
    
    force_new = compute_force(pos, mass)
    vel += 0.5 * (force + force_new) / mass[:, None] * dt
    
    return pos, vel

def rk4_step(pos, vel, mass, dt):
    """RK4 (4th order)"""
    def derivative(state, t):
        pos_flat, vel_flat = state[:len(pos)*2], state[len(pos)*2:]
        pos_reshaped = pos_flat.reshape(-1, 2)
        vel_reshaped = vel_flat.reshape(-1, 2)
        
        force = compute_force(pos_reshaped, mass)
        dpos = vel_reshaped.flatten()
        dvel = (force / mass[:, None]).flatten()
        
        return np.concatenate([dpos, dvel])
    
    state = np.concatenate([pos.flatten(), vel.flatten()])
    
    k1 = derivative(state, 0)
    k2 = derivative(state + 0.5 * dt * k1, 0.5 * dt)
    k3 = derivative(state + 0.5 * dt * k2, 0.5 * dt)
    k4 = derivative(state + dt * k3, dt)
    
    new_state = state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    n_particles = len(pos)
    new_pos = new_state[:n_particles*2].reshape(-1, 2)
    new_vel = new_state[n_particles*2:].reshape(-1, 2)
    
    return new_pos, new_vel

def integrate(integrator, pos, vel, mass, dt, total_time):
    """Integrate using specified method"""
    steps = int(total_time / dt)
    
    for _ in range(steps):
        pos, vel = integrator(pos, vel, mass, dt)
    
    return pos, vel

def create_binary_system():
    """Create a stable binary star system for testing"""
    # Two stars in circular orbit
    pos = np.array([
        [0.5, 0.0],   # Star 1
        [-0.5, 0.0]   # Star 2
    ])
    
    mass = np.array([1.0, 1.0])
    
    # Circular orbit velocities
    r = np.linalg.norm(pos[1] - pos[0])
    v_circ = np.sqrt(G * mass[0] / r)
    
    vel = np.array([
        [0.0, v_circ],    # Star 1
        [0.0, -v_circ]    # Star 2
    ])
    
    return pos, vel, mass

def compute_reference_solution(dt_ref, total_time):
    """Compute high-accuracy reference solution using RK4 with very small dt"""
    pos, vel, mass = create_binary_system()
    
    # Use RK4 with very small time step for reference
    pos_ref, vel_ref = integrate(rk4_step, pos.copy(), vel.copy(), mass, dt_ref, total_time)
    
    return pos_ref, vel_ref

def compute_convergence_rate(errors, dt_values):
    """Compute convergence rate using linear regression"""
    log_dt = np.log(dt_values)
    log_error = np.log(errors)
    p, _ = np.polyfit(log_dt, log_error, 1)
    return p

def main():
    print("N-BODY CONVERGENCE ORDER VALIDATION")
    print("=" * 50)
    print("Testing on actual binary star system")
    print("Only position errors, reference lines aligned with methods")
    print()
    
    # Create binary system
    pos0, vel0, mass = create_binary_system()
    total_time = 1.0  # Simulation time
    
    # Reference solution (high accuracy)
    dt_ref = 0.0001
    pos_ref, vel_ref = compute_reference_solution(dt_ref, total_time)
    
    # Time step values for convergence testing
    dt_values = np.array([0.1, 0.05, 0.025, 0.0125, 0.00625])
    
    methods = {
        "Symplectic Euler": symplectic_euler_step,
        "Velocity Verlet": velocity_verlet_step,
        "RK4": rk4_step
    }
    
    theoretical_orders = {
        "Symplectic Euler": 1,
        "Velocity Verlet": 2,
        "RK4": 4
    }
    
    results = {}
    
    print("Computing convergence rates...")
    print("-" * 40)
    
    for name, method in methods.items():
        position_errors = []
        
        for dt in dt_values:
            # Integrate with current method
            pos, vel = integrate(method, pos0.copy(), vel0.copy(), mass, dt, total_time)
            
            # Compute position error compared to reference
            pos_error = np.linalg.norm(pos - pos_ref)
            position_errors.append(pos_error)
        
        # Compute convergence rate
        pos_order = compute_convergence_rate(np.array(position_errors), dt_values)
        
        results[name] = {
            'position_errors': position_errors,
            'position_order': pos_order,
            'theoretical_order': theoretical_orders[name]
        }
        
        print(f"{name}:")
        print(f"  Position order: {pos_order:.3f} (theoretical: {theoretical_orders[name]})")
        print()
    
    # Create convergence plot
    plt.figure(figsize=(10, 8))
    
    colors = {
        "Symplectic Euler": NBody_Palette["blue_light"],
        "Velocity Verlet": NBody_Palette["blue_deep"],
        "RK4": NBody_Palette["orange"]
    }
    
    # Plot actual errors
    for name, data in results.items():
        color = colors[name]
        
        plt.loglog(dt_values, data['position_errors'], 'o-', 
                  color=color, linewidth=3, markersize=10, 
                  label=f'{name} (order: {data["position_order"]:.2f})')
    
    # Add reference lines aligned with each method
    for name, data in results.items():
        color = colors[name]
        theoretical_order = data['theoretical_order']
        
        # Create reference line starting from the first data point
        ref_dt = np.array([dt_values[0], dt_values[-1]])
        ref_error = data['position_errors'][0] * (ref_dt / dt_values[0])**theoretical_order
        
        plt.loglog(ref_dt, ref_error, '--', 
                  color=color, alpha=0.7, linewidth=2,
                  label=f'{theoretical_order} order reference')
    
    plt.xlabel('Time step Δt', fontsize=14)
    plt.ylabel('Position error', fontsize=14)
    plt.title('N-Body Convergence Order Validation', fontsize=16, fontweight='bold')
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    
    # Add theoretical order annotations
    for name, data in results.items():
        color = colors[name]
        theoretical_order = data['theoretical_order']
        
        # Position annotation near the middle of the line
        mid_dt = np.sqrt(dt_values[0] * dt_values[-1])
        mid_error = data['position_errors'][0] * (mid_dt / dt_values[0])**theoretical_order
        
        plt.annotate(f'order {theoretical_order}', 
                    xy=(mid_dt, mid_error),
                    xytext=(10, 10), textcoords='offset points',
                    fontsize=10, color=color,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=color, alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('nbody_convergence_corrected.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print final results
    print("\n" + "=" * 80)
    print("FINAL CONVERGENCE RESULTS")
    print("=" * 80)
    print(f"{'Method':<20} {'Theoretical':<12} {'Empirical':<12} {'Error':<12} {'Status':<10}")
    print("-" * 80)
    
    for name, data in results.items():
        error = abs(data['position_order'] - data['theoretical_order'])
        status = "✅" if error < 0.3 else "⚠️" if error < 0.8 else "❌"
        
        print(f"{name:<20} {data['theoretical_order']:<12} {data['position_order']:<12.3f} {error:<12.3f} {status:<10}")
    
    print("\n" + "=" * 100)
    print("VALIDATION METHODOLOGY")
    print("=" * 100)
    print("1. Uses actual N-body binary system (not simplified oscillator)")
    print("2. Reference solution: RK4 with Δt = 0.0001 (high accuracy)")
    print("3. Only position errors (as requested)")
    print("4. Reference lines aligned with corresponding methods")
    print("5. Linear regression on log(error) vs log(Δt) for empirical order")

if __name__ == "__main__":
    main()