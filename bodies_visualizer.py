import sys
import matplotlib.pyplot as plt

def visualize(input_filename):
    input_file = open(input_filename, "r")
    bodies_string = input_file.read()
    input_file_arr = bodies_string.splitlines()

    x_coords = []
    y_coords = []
    num_bodies = int(input_file_arr[1])

    for i in range(num_bodies):
        x, y = input_file_arr[i + 2].split(" ")
        x_coords.append(int(x))
        y_coords.append(int(y))
    
    plt.plot(x_coords, y_coords, 'ro')
    plt.show()

if __name__ == '__main__':
    bodies_input = sys.argv[1]

    visualize(bodies_input)
