'''
from matplotlib import pyplot as plt
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

fig = plt.figure()
ax = p3.Axes3D(fig)

# create the parametric curve
t=np.arange(0, 2*np.pi, 2*np.pi/100)
x=np.cos(t)
y=np.sin(t)
z=t/(2.*np.pi)

# create the first plot
point, = ax.plot([x[0]], [y[0]], [z[0]], 'o')
line, = ax.plot(x, y, z, label='parametric curve')
ax.legend()
ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])

# second option - move the point position at every frame
def update_point(n, x, y, z, point):
    point.set_data(np.array([x[n], y[n]]))
    point.set_3d_properties(z[n], 'z')
    return point

ani=animation.FuncAnimation(fig, update_point, 99, fargs=(x, y, z, point))

plt.show()
'''
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
import pandas as pd

a = np.random.rand(2000, 3) * 10
t = np.array([np.ones(100) * i for i in range(20)]).flatten()
df = pd.DataFrame({"time": t, "x": a[:, 0], "y": a[:, 1], "z": a[:, 2]})


def update_graph(num):
    data = df[df['time'] == num]
    graph.set_data(data.x, data.y)
    graph.set_3d_properties(data.z)
    title.set_text('3D Test, time={}'.format(num))
    return title, graph,


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('3D Test')

data = df[df['time'] == 0]
graph, = ax.plot(data.x, data.y, data.z, linestyle="", marker="o")

ani = matplotlib.animation.FuncAnimation(fig, update_graph, 19,
                                         interval=40, blit=True)

plt.show()
