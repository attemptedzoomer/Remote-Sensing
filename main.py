import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scene import generate_scene

plt.style.use('dark_background')

fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.add_subplot(1, 1, 1)
k = 0
lxb = -0.5
rxb = 1.5
domain_density = 1000
nx = int((rxb - lxb) * domain_density)
x0 = np.linspace(lxb, rxb, nx)
a = generate_scene(lxb, rxb, nx)


def animate(i):
    global k
    x = a[k]
    k += 1
    ax1.clear()
    plt.plot(x0, x, color='cyan')
    plt.ylim([-2, 2])
    dll, dlr = 0, 1
    plt.xlim([0, 1])

    maximum, minimum = 0, 0
    for i in range(int((dll - lxb) * domain_density), int((dlr - lxb) * domain_density)):
        try:
            if a[k][i] > maximum:
                maximum = a[k][i]
            if a[k][i] < minimum:
                minimum = a[k][i]
        except:
            maximum, minimum = 0, 0

    plt.annotate(f'Max amplitude: {round(maximum, 4).__format__(".04f")}\n '
                 f'Min Amplitude: {round(minimum, 4).__format__(".04f")}', (0.25, -0.75))


with open("Output/attemptdata.txt", "r") as f:
    contents = f.readlines()
    attempt = int(contents[0])

with open("Output/attemptdata.txt", "w") as f:
    f.write(f"{attempt + 1}")

anim = animation.FuncAnimation(fig, animate, frames=len(a) - 1, interval=20)
anim.save(f"Output/Wave{format(attempt, '04d')}.gif", writer='ffmpeg', fps=50)
