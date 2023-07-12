import subprocess

def run_shell(cmd):
    status = 0
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        status = e.returncode
        output = e.output
    output = output.decode()
    return (status, output)

def is_in_path(path):
    try:
        if subprocess.run(["which",path],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).returncode != 0:
            return False
    except:
        return False
    return True