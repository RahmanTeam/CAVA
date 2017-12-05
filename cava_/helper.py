
import os

# Reads default config file path from the default_config_path file
def defaultConfigPath():
    d = os.path.dirname(os.path.realpath(__file__))
    dir = d[:d.rfind('env/lib')]
    for line in open(dir + '/default_config_path'):
        line = line.strip()
        if line != '': return line
    return None
