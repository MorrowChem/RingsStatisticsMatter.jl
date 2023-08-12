import numpy as np


def halton_sequence(dim, num_points, ats):
    vol = ats.get_volume()
    def halton_base_seq(base):
        seq = []
        for i in range(num_points):
            x = (20**3 / vol)**(1/3) * 0.2
            f = 1.0 / base
            while i > 0:
                x += f * (i % base)
                i //= base
                f /= base
            seq.append(x)
        return seq
    
    sequence = []
    for d in range(dim):
        base = 2 + d  # Use different prime bases for each dimension
        sequence.append(halton_base_seq(base))
    
    return np.transpose(sequence)

def place_points_in_cube(num_points, ats):
    dim = 3  # Number of dimensions
    points = halton_sequence(dim, num_points, ats)
    return points