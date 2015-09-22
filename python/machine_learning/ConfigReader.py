
import yaml

class ConfigReader: 

    """ Reader used to parse and store hps-learn configuration files  """

    def __init__(self, file_name):
        self.config = self.parse_config(str(file_name))

    def parse_config(self, file_name):
        print "Loading configuration from " + str(file_name)
        config_file = open(file_name, 'r')
        return yaml.load(config_file)
    
    def get_signal_files(self):
        return self.config["Signal Files"]

    def get_background_files(self):
        return self.config["Background Files"]

    def get_training_variables(self):
        return self.config["Variables"]
