#!/usr/bin/env python3

'''
Usage:
  docker_run.py <file>
'''

from docker import Client
import sys
import os
from warnings import warn

# Parse the configuration file that contains a dictionary
with open(sys.argv[1]) as fin:
    config = eval(fin.read())

# Instantiate the client that will communicate with the Docker daemon
cli = Client(base_url='unix://var/run/docker.sock')

# Pull the lastest image
for line in cli.pull(repository=config['image'], tag=config['tag'], stream=True):
    sys.stdout.write(line.decode(sys.stdout.encoding))

# Get the list of mountpoints and declare volume mappings
volumes = []
for volume in config['volumes']:
    volumes.append(volume.split(':')[1])
host_config = cli.create_host_config(binds=config['volumes'])

# Create a container and start it
container = cli.create_container(image=config['image'] + ':' + config['tag'],
                                 command='tail -f /dev/null',
                                 detach=True,
                                 stdin_open=True,
                                 tty=True,
                                 environment=config['environment'],
                                 volumes=volumes,
                                 name=config['name'],
                                 host_config=host_config)
# Forward warning messages to stderr if any
if container.get('Warnings') is not None:
    warn(container.get('Warnings'), RuntimeWarning)
# Start the container
cli.start(container=container.get('Id'))

# Execute the commands
for cmd in config['cmd']:
    print('[+] ' + cmd)
    execute = cli.exec_create(container['Id'], cmd=cmd, stdout=True, stderr=True)
    for char in cli.exec_start(execute['Id'], tty=True, stream=True):
        sys.stdout.write(char.decode(sys.stdout.encoding))
    status = cli.exec_inspect(execute['Id'])['ExitCode']
    if status != 0:
        break

# Stop the container and remove it
cli.stop(container=container.get('Id'))
cli.remove_container(container=container['Id'])

sys.exit(status)
