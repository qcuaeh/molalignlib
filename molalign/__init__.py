try:
    from .wrapper import Alignment
except ModuleNotFoundError:
    from os import path
    from shutil import copyfile
    from subprocess import Popen, PIPE, STDOUT
    rootdir = path.dirname(__file__)
    source = path.join(rootdir, 'config', 'gnu.cfg')
    destination = path.join(rootdir, 'build.cfg')
    command = [path.join(rootdir, 'build.sh'), '-lpy']
    try:
        copyfile(source, destination)
    except FileExistsError:
        pass
    with Popen(command, stdout=PIPE, stderr=STDOUT, bufsize=0) as p:
        for line in p.stdout:
            print(line.decode('utf-8').rstrip())
