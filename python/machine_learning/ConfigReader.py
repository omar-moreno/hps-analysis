
import yaml

class ConfigReader: 

    """ Reader used to parse and store hps-learn configuration files  """

    def __init__(self, file_name):
        self.parse_config(file_name)

    def parse_config(self, file_name):
        print "Loading configuration from " + file_name
        config_file = open(file_name, 'r')
        yaml.load(config_file)
        

