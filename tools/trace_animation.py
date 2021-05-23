import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tqdm import trange
from pathlib import Path

status_color = {
    'NOT_EXPOSED': 'g',
    'EXPOSED': 'y',
    'INFECTED': 'r',
    'IMMUNE': 'b',
}


def load_data(res_dir: Path, countries):
    # Parse and merge csv files
    df = pd.DataFrame()
    for c in range(countries):
        filepath = res_dir.joinpath(f'trace_{c}.csv')
        df_c = pd.read_csv(filepath, index_col=['t', 'id'])
        df = df.append(df_c)

    return df


def main(res_dir: Path, countries, world_w, world_l, country_w, country_l, t_step, t_target=None):
    print("Loading data...")
    df = load_data(res_dir, countries)

    # Compute additional parameters from given data
    cols, rows = int(world_w // country_w), int(world_l // country_l)
    if not t_target:
        t_target = int(df.index.get_level_values('t').max())

    # Draw the initial situation
    u = 3
    fig, ax = plt.subplots(figsize=(u * cols, u * rows))
    ax.set(xlim=(0, world_w), ylim=(0, world_l))
    ax.set_xticks(np.linspace(0, world_w, cols + 1))
    ax.set_yticks(np.linspace(0, world_l, rows + 1))
    ax.grid()
    del u

    print("Initial plot...")
    # Define context for the animate clojure
    df_t = df.loc[0]
    scat = ax.scatter(
        'pos_x',
        'pos_y',
        c=df_t['status'].map(status_color),
        data=df_t
    )

    print("Animating...")

    # Define the animation function
    def animate(t):
        df_t = df.loc[t]
        scat.set_offsets(df_t[['pos_x', 'pos_y']])
        scat.set_color(df_t['status'].map(status_color))
        if t in []:
            print(f'Snapshot {t}')
            plt.savefig(f'anim_f_{t}.pdf')
        return scat

    # Set up the animation
    anim = FuncAnimation(
        fig,
        animate,
        frames=trange(t_step, t_target, t_step),
        init_func=lambda: scat
    )

    # Render and save the animation
    anim.save('anim.mp4', dpi=200)


if __name__ == '__main__':
    # Set up the parameters
    res_dir = Path.cwd().joinpath('../src/results')
    countries = 4
    world_w, world_l = 2e3, 2e3
    country_w, country_l = 1e3, 1e3
    t_step = 1
    t_target = 600

    # Call the main function
    main(res_dir, countries, world_w, world_l,
         country_w, country_l, t_step, t_target)

# Suggested running parameters (for make run)
# @mpirun -np 4 --oversubscribe $(exec) \
# 	-N 1000 \
# 	-I 100 \
# 	-W 2000 \
# 	-L 2000 \
# 	-w 1000 \
# 	-l 1000 \
# 	-v 1.4 \
# 	-d 30 \
# 	--t-infection=$$((1 * 10)) \
# 	--t-recovery=$$((2 * 60)) \
# 	--t-immunity=$$((2 * 60)) \
# 	--sim-step=1 \
# 	--sim-length=1 \
# 	--write-trace \
# 	--log-level INFO
