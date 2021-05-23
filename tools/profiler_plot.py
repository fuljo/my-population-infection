import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def latex_time(seconds):
    h = seconds // 3600
    m = (seconds % 3600) // 60
    s = seconds % 60

    h_str = f'\SI{{{int(h):2d}}}{{h}}~' if h > 0 else ''
    m_str = f'\SI{{{int(m):2d}}}{{m}}~' if m > 0 or h > 0 else ''
    s_str = f'\SI{{{int(s):2d}}}{{s}}'

    return f'${h_str}{m_str}{s_str}$'


def plot(csvpath: Path, var: str, fit=True, report_dir=None):
    print(f'{var}...')

    df = pd.read_csv(csvpath, header=0)

    if report_dir:
        # Convert elapsed time to latex hour-minute-second format
        df['elapsed_latex'] = df['elapsed'].apply(latex_time)
        df.to_csv(Path(report_dir).joinpath(csvpath.name), index=False)

    if fit:
        from sklearn.linear_model import LinearRegression

        # Fit to an inverse-quadratic
        reg = LinearRegression()
        x = df[var].to_numpy()
        if var == 'countries':
            x = 1 / (x ** 2)
        elif var == 'individuals':
            x = x ** 2
        x = np.expand_dims(x, -1)
        y = df['elapsed'].to_numpy()
        reg.fit(x, y)
        print(f'Intercept: {reg.intercept_}, Coefficients: {reg.coef_}')

        df['predicted'] = reg.predict(x)

    # Draw total time vs. num of individuals
    fig, ax = plt.subplots()
    ax.plot(var, 'elapsed', data=df)
    if fit:
        ax.plot(var, 'predicted', '--', c='tab:gray', data=df)
    ax.set_xlabel(var)
    ax.set_ylabel('exec. time (s)')
    ax.legend()

    fig.savefig(f'{csvpath.stem}.png', dpi=200)


if __name__ == '__main__':
    report_dir = None  # Change this to './report' to rewrite latex results
    for path in Path.cwd().glob('profile_countries*.csv'):
        plot(path, 'countries', report_dir=report_dir)
    for path in Path.cwd().glob('profile_individuals*.csv'):
        plot(path, 'individuals', report_dir=report_dir)
