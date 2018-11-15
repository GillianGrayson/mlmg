from config.types.attributes.attribute import Attribute
from infrastructure.load.attributes import get_attributes


def get_attributes_dict(config):

    attributes = get_attributes(config)

    min_atr = min(attributes)
    max_atr = max(attributes)

    min_atr = int(min_atr / config.shift) * config.shift
    max_atr = (int(max_atr / config.shift) + 1) * config.shift
    attributes_dict = {}
    for age_id in range(0, len(attributes)):
        age = attributes[age_id]
        key = int((age - min_atr) / config.shift)
        if key in attributes_dict:
            attributes_dict[key].append(age_id)
        else:
            attributes_dict[key] = [age_id]

    return attributes_dict

def get_attributes_group(config, attribute=Attribute.age):

    attributes = get_attributes(config, attribute)

    min_atr = min(attributes)
    max_atr = max(attributes)

    min_atr = int(min_atr / config.shift) * config.shift
    max_atr = (int(max_atr / config.shift) + 1) * config.shift
    group = []
    for age_id in range(0, len(attributes)):
        age = attributes[age_id]
        key = int((age - min_atr) / config.shift)
        group.append(key)

    return group