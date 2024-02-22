from matplotlib import pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.animation as animation

def plot3d(name, xStep, yStep):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    x = np.arange(0, 1 + xStep / 2, xStep)
    y = np.arange(0, 1 + yStep / 2, yStep)
    x, y = np.meshgrid(x, y)
    u = np.genfromtxt('simulation/' + name + ".csv", delimiter=',')

    surf = ax.plot_surface(x, y, u, cmap=cm.coolwarm, linewidth=0, antialiased=False)

    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title(name)
    plt.show()


def plotAnimate(name, xStep, yStep, duration=5):
    fig, ax = plt.subplots()

    x = np.arange(0, 1 + xStep / 2, xStep)
    y = np.arange(0, 1 + yStep / 2, yStep)
    u = np.genfromtxt('simulation/' + name + ".csv", delimiter=',')

    minVal = min(u[0])
    maxVal = max(u[0])
    
    def doFrame(snapshot):
        ax.clear()
        ax.set_ylim(minVal * 1.1, maxVal * 1.1)
        ax.plot(x, snapshot, color='blue')
        ax.scatter(x, [0] * len(x), c=cm.coolwarm((snapshot - minVal) / (maxVal - minVal)))

    plt.title(name)

    ani = animation.FuncAnimation(fig, doFrame, frames=u, interval=1000*duration/len(y))
    writer = animation.PillowWriter(fps=len(y) / duration,
                                    metadata={},
                                    bitrate=1800)
    ani.save('simulation/' + name + '.gif', writer=writer)
    plt.show()


#plot3d('heatflow', 0.01, 0.01)
#plotAnimate('heatflow', 0.01, 0.01)

#plot3d('heatflowCN', 0.01, 0.01)
#plotAnimate('heatflowCN', 0.01, 0.01)

#plot3d('poisson5', 1.0 / 6, 1.0 / 6)

#plot3d('poisson9', 1.0 / 6, 1.0 / 6)

#plot3d('transport', 0.01, 0.01)
plotAnimate('transport', 0.01, 0.01)