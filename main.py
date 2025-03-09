import matplotlib.pyplot as plt
import math
import os
import sys

class CyclicVoltammetryToolbox:
    def __init__(self, file_name):
        self.file_name = file_name
        self.data_points = []
        self.oxidation_peak = []
        self.reduction_peak = []
        self.concentration = None
        self.variables = {}

# method to get the file path into the data_files folder
    def get_data_file_path(self, file_name):
        if hasattr(sys, '_MEIPASS'):
            base_dir = os.path.join(sys._MEIPASS, 'data_files')
        else:
            base_dir = os.path.join(os.getcwd(), 'data_files')
        return os.path.join(base_dir, file_name)
# method to load the data from the file
    def load_data(self):
        file_path = self.get_data_file_path(self.file_name)
        try:
            with open(file_path, "r") as file:
                return file.read()
        except Exception as e:
            print(f"Error loading data from {file_path}: {e}")
            return None
# method to parse the data from the file
    def parse_data(self, content):
        data_points = []
        for line in content.split('\n'):
            if line.strip():
                parts = line.split()
                if len(parts) == 2:
                    try:
                        data_points.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        continue
        self.data_points = data_points

    def find_peaks(self):
        oxidation_peak = []
        reduction_peak = []
        increasing = None  # Start with no assumption

        for i in range(1, len(self.data_points) - 1):
            prev_point, curr_point, next_point = self.data_points[i - 1], self.data_points[i], self.data_points[i + 1]

            if increasing is None:  # Determine initial trend
                if curr_point[0] < next_point[0]:
                    increasing = True
                else:
                    increasing = False

            if increasing:
                if curr_point[1] > prev_point[1] and curr_point[1] > next_point[1]:
                    oxidation_peak.append(curr_point)
                elif curr_point[0] > next_point[0]:  # Switch to reduction detection
                    increasing = False

            if not increasing:
                if curr_point[1] < prev_point[1] and curr_point[1] < next_point[1]:
                    reduction_peak.append(curr_point)
                elif curr_point[0] < next_point[0]:  # Switch back to oxidation detection
                    increasing = True

        self.oxidation_peak = oxidation_peak
        self.reduction_peak = reduction_peak

    def generate_graph(self, output_file):
        try:
            plt.figure()
            potentials, currents = zip(*self.data_points)
            plt.plot(potentials, currents, label='Cyclic Voltammogram')
            if self.oxidation_peak:
                max_potentials, max_currents = zip(*self.oxidation_peak)
                plt.scatter(max_potentials, max_currents, color='red', label='Oxidation Peak')
                for peak in self.oxidation_peak:
                    plt.annotate(f"({peak[0]:.2f}, {peak[1]:.2f})", (peak[0], peak[1]), textcoords="offset points", xytext=(0,10), ha='center')
            if self.reduction_peak:
                min_potentials, min_currents = zip(*self.reduction_peak)
                plt.scatter(min_potentials, min_currents, color='blue', label='Reduction Peak')
                for peak in self.reduction_peak:
                    plt.annotate(f"({peak[0]:.2f}, {peak[1]:.2f})", (peak[0], peak[1]), textcoords="offset points", xytext=(0,10), ha='center')
            plt.xlabel('Potential (V)')
            plt.ylabel('Current ($\mu$A)')
            plt.legend()
            plt.show()

            plt.savefig(output_file)

            print(f"Graph successfully saved at: {output_file}")
            plt.close()
        except Exception as e:
            raise

    def calculate_concentration(self):
        n = int(input("Enter the number of electrons transferred (n): "))
        radius_mm = float(input("Enter the radius of the electrode (mm): "))
        radius_cm = radius_mm / 10
        area = math.pi * radius_cm**2
        print(f"Calculated electrode area (A): {area:.4f} cm²")
        diffusion_file_name = input("Enter the name of the file containing the diffusion coefficient: ")
        print(f"Reading diffusion coefficient file: {diffusion_file_name}")
        ion_name = input("Enter the name of the ion: ")
        
        try:
            diffusion_file_path = self.get_data_file_path(diffusion_file_name)
            with open(diffusion_file_path, "r") as file:
                diffusion_content = file.read()
                diffusion_coefficient = self.parse_diffusion_coefficient(diffusion_content, ion_name)
                diffusion_coefficient_sf = diffusion_coefficient * 1e-6
                print(f"Diffusion coefficient (D) for {ion_name}: {diffusion_coefficient_sf:.4e} cm²/s")
        except Exception as e:
            print(f"Error reading diffusion coefficient file: {e}")
            return
        
        scan_rate = float(input("Enter the scan rate (V/s): "))
        if self.oxidation_peak:
            peak_current = self.oxidation_peak[0][1] * 1e-6
            print(f"Peak current (I_p): {peak_current:.4e} A")
            concentration = peak_current / (2.69e5 * n**(3/2) * area * diffusion_coefficient_sf**(1/2) * scan_rate**(1/2))
            print(f"Calculated concentration (C) of {ion_name} in the solution: {concentration:.4e} mol/cm³")
            self.concentration = concentration
            self.variables = {
                "n": n,
                "radius_mm": radius_mm,
                "radius_cm": radius_cm,
                "area": area,
                "diffusion_coefficient_sf": diffusion_coefficient_sf,
                "scan_rate": scan_rate,
                "peak_current": peak_current,
                "ion_name": ion_name
            }
        else:
            print("No maxima found for peak current calculation.")

    def parse_diffusion_coefficient(self, content, ion_name):
        for line in content.split('\n'):
            if ion_name in line:
                for part in line.split():
                    try:
                        return float(part)
                    except ValueError:
                        continue
        raise ValueError(f"Diffusion coefficient for {ion_name} not found.")

def main():
    file_path = input("Enter the absolute path to the Cyclic Voltammogram (.txt) file: ")
    cv_toolbox = CyclicVoltammetryToolbox(file_path)
    
    content = cv_toolbox.load_data()
    if content:
        cv_toolbox.parse_data(content)
        cv_toolbox.find_peaks()

        output_graph_path = os.path.join(
            os.path.dirname(cv_toolbox.get_data_file_path(file_path)),
            f"{os.path.basename(file_path).split('.')[0]}_graph.png")
        print(f"Attempting to save graph at: {output_graph_path}")
        cv_toolbox.generate_graph(output_graph_path)
        cv_toolbox.calculate_concentration()

if __name__ == '__main__':
    main()

