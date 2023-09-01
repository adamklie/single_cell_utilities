import os
import subprocess
import loompy as lp
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def show_values(axs, orient="v", space=.01):
    """Show values on top of bars in barplot
    
    Parameters
    ----------
    axs : matplotlib.axes.Axes
        Axes object
    orient : str
        Orientation of barplot, "v" for vertical, "h" for horizontal
    space : float
        Space between bar and value
    
    Returns
    -------
    None
    """
    def _single(ax):
        if orient == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height() + (p.get_height()*0.01)
                value = '{:.1f}'.format(p.get_height())
                ax.text(_x, _y, value, ha="center")
        elif orient == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height() - (p.get_height()*0.5)
                value = '{:.1f}'.format(p.get_width())
                ax.text(_x, _y, value, ha="left")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _single(ax)
    else:
        _single(axs)


def make_dirs(path):
    """Make directory if it doesn't exist
    
    Parameters
    ----------
    path : str
        Path to directory
        
    Returns
    -------
    str
        Path to directory
    """
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],  # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid)  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    t0 = get_ticks()
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    t1 = get_ticks()
    err_str = (
        'PID={pid}, PGID={pgid}, RC={rc}, DURATION_SEC={dur:.1f}\n'
        'STDERR={stde}\nSTDOUT={stdo}'
    ).format(
        pid=pid, pgid=pgid, rc=rc, dur=t1 - t0, stde=stderr.strip(), stdo=stdout.strip()
    )
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


def is_counts(data):
    if np.all(data >= 0) and np.all(data.astype(int) == data):
        print("The matrix contains count data.")
    else:
        print("The matrix does not contain count data.")
        
        
def is_mostly_counts(data, percent=0.9):
    """Want to check if some percent of the data is counts

    Args:
        data (_type_): _description_
    """
    if np.all(data >= 0) and np.all(data.astype(int) == data):
        print("The matrix contains all count data.")
    elif np.sum(data >= 0) / data.size >= percent and np.sum(data.astype(int) == data) / data.size >= percent:
        greater_than_0 = (np.sum(data >= 0) / data.size)*100
        int_equals = (np.sum(data.astype(int) == data) / data.size)*100
        print(f"The matrix contains mostly count data. {greater_than_0}% of the data is greater than 0 and {int_equals}% of the data is equal to its integer value.")
    else:
        print("The matrix does not contain mostly count data.")