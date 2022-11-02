import matplotlib.pyplot as plt
import numpy as np

class Beam:

    def __init__(self, height=1, base=1, length=1, resolution=1e-2):
        """
        Using metric system
        """
        self.height = height
        self.base = base
        self.length = length
        self.stress_invariant = self.compute_stress_invariant()
        self.resolution = resolution
        
        self.len_matrix_length = int(self.length/self.resolution)
        self.len_matrix_height = int(self.height/self.resolution)

        self.sigma_x_matrix = np.zeros((self.len_matrix_height, self.len_matrix_length))
        self.tau_xy_matrix = np.zeros((self.len_matrix_height, self.len_matrix_length))

    def compute_stress_invariant(self):
        return (1/12)*self.base*(self.height**3)


    def compute_sigma_x(self, x, y, moment_function):
        return -(moment_function(x)*y)/(self.stress_invariant)

    def compute_tau_xy(self, V, Q, t):
        return -(V*Q(length=self.length, thickness=t, pos_above_y=self.height/2))/(self.stress_invariant*self.height) 

    def compute_full_beam_sigma_x(self):
        for x_i in range(self.len_matrix_length):
            for y_i in range(self.len_matrix_height):
                x_position = x_i*(self.length/self.resolution)
                y_position = y_i*(self.height/self.resolution)
                self.sigma_x_matrix[y_i][x_i] = self.compute_sigma_x(x_position, y_position, self.sample_moment_function)

    def compute_full_beam_tau_xy(self):
        for x_i in range(self.len_matrix_length):
            for y_i in range(self.len_matrix_height):
                x_position = x_i*(self.length/self.resolution)
                y_position = y_i*(self.height/self.resolution)
                self.tau_xy_matrix[y_i][x_i] = self.compute_tau_xy(V=10, Q=self.compute_Q, t=self.height)                

    def export_sigma_x_heatmap(self, filepath="sample-sigma-x.png"):
        plt.pcolormesh(self.sigma_x_matrix)
        plt.colorbar()
        plt.savefig(fname=filepath)
        plt.close()

    def export_tau_xy_heatmap(self, filepath="sample-tau-xy.png"):
        plt.pcolormesh(self.tau_xy_matrix)
        plt.colorbar()
        plt.savefig(fname=filepath)
        plt.close()

    def compute_Q(self, length, thickness, pos_above_y):
        return length*thickness*pos_above_y

    def sample_moment_function(self, x):
        """
        Constant force of 10 N applied over a range of 0 meters to n-distance meters.
        """
        return 0.5*(x**2)*10

if __name__ == "__main__":
    beam = Beam(height=0.2, length=2)
    beam.compute_full_beam_sigma_x()
    beam.compute_full_beam_tau_xy()
    beam.export_sigma_x_heatmap()
    beam.export_tau_xy_heatmap()
