from config.types import *
from config.method import *


def get_info_types(config):
    info_types = [x.value for x in InfoType]
    return info_types