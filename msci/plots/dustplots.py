import matplotlib

import matplotlib.pyplot as plt
import matplotlib.animation
import mpl_toolkits.mplot3d.axes3d as p3

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation


def pplot(dustanalysis, name='Dust_rotation', save=False, jn=False):
    def _update_graph(n_iter):
        data = dustanalysis.positions_df[dustanalysis.positions_df['time'] == n_iter]
        point.set_data(data.x, data.y)
        point.set_3d_properties(data.z)
        return point,

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=90., azim=90)
    if min(dustanalysis.position_array[:, 0]) == max(dustanalysis.position_array[:, 0]):
        ax.set_xlim([-10 * dustanalysis.const.lambdaD, 10 * dustanalysis.const.lambdaD])
    else:
        ax.set_xlim([min(dustanalysis.position_array[:, 0]), max(dustanalysis.position_array[:, 0])])
    if min(dustanalysis.position_array[:, 1]) == max(dustanalysis.position_array[:, 1]):
        ax.set_ylim([-10 * dustanalysis.const.lambdaD, 10 * dustanalysis.const.lambdaD])
    else:
        ax.set_ylim([min(dustanalysis.position_array[:, 1]), max(dustanalysis.position_array[:, 1])])
    if min(dustanalysis.position_array[:, 2]) == max(dustanalysis.position_array[:, 2]):
        ax.set_zlim([-10 * dustanalysis.const.lambdaD, 10 * dustanalysis.const.lambdaD])
    else:
        ax.set_zlim(
            [min(dustanalysis.position_array[:, 2]) - 0.5 * min(dustanalysis.position_array[:, 2]),
             max(dustanalysis.position_array[:, 2])])

    ax.set_xlim([-dustanalysis.modified_b_field['rmax'] * 1.5, dustanalysis.modified_b_field['rmax'] * 1.5])
    ax.set_ylim([-dustanalysis.modified_b_field['rmax'] * 1.5, dustanalysis.modified_b_field['rmax'] * 1.5])
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")

    data = dustanalysis.positions_df[dustanalysis.positions_df['time'] == 0]
    point, = ax.plot(data.x, data.y, data.z, linestyle="", marker=".")

    if not jn:
        plt.cla()

    ani = matplotlib.animation.FuncAnimation(
        fig,
        _update_graph,
        frames=dustanalysis.iterationsB + dustanalysis.init_iterations,
        interval=1,
        blit=True
    )

    if save:
        ani.save('{}.mp4'.format(name), fps=30, extra_args=['-vcodec', 'libx264'])

    if jn:
        return ani
    else:
        fig.show()


def plotgibson2d(dustanalysis, save):
    font = {'family': 'normal',
            'weight': 'bold',
            'size': 25}

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    data = dustanalysis.positions_df[dustanalysis.positions_df['time'] == max(dustanalysis.positions_df['time'])]
    ax2.scatter(data.x, data.y, color='r', marker=".")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.title("Crystal without Gibson's modified E field")
    plt.xlim([-max(data.x) * 1.1, max(data.x) * 1.1])
    plt.ylim([-max(data.y) * 1.1, max(data.y) * 1.1])
    if save:
        plt.savefig('fig', format='eps')
    plt.show()
