import os
import re
import subprocess
import signal


def strip_ext_gz(f):
    return re.sub(r'\.gz$', '', str(f))


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
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid, pgid, rc, stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    return stdout.strip('\n')


def gunzip(f, suffix, out_dir):
    if not f.endswith('.gz'):
        raise Exception('Cannot gunzip a file without .gz extension.')
    gunzipped = os.path.join(out_dir,
                             os.path.basename(strip_ext_gz(f)))
    if suffix:
        gunzipped += '.{}'.format(suffix)
    # cmd = 'gzip -cd {} > {}'.format(f, gunzipped)
    cmd = 'zcat -f {} > {}'.format(f, gunzipped)
    run_shell_cmd(cmd)
    return gunzipped


def rm_f(files):
    if files:
        if type(files) == list:
            run_shell_cmd('rm -f {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -f {}'.format(files))