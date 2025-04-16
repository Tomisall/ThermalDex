import matplotlib.pyplot as plt

def scatterSelection(x_values, y_values, x_name, y_name, current_x, current_y):
    plt.ion()
    plt.figure(num='ThermalDex Database Scatter Plot')
    #plt.close()
    plt.clf()
    plt.plot(x_values, y_values, 'bo', current_x, current_y, 'ro')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.show()