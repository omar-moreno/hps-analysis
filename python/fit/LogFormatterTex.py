
import re
from matplotlib.ticker import LogFormatter

class LogFormatterTex(LogFormatter, object) : 

    def __init__(self, *args, **kwargs) : 
        super(LogFormatterTex, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs) :
        label = super(LogFormatterTex, self).__call__(*args, **kwargs)
        label = re.sub(r'e(\S)0?(\d+)',
                       r'\\times 10^{\1\2}',
                       str(label))
        label = "$" + label + "$"
        return label


