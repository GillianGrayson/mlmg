from config.types.attributes.attribute import Attribute


def age_less(config, limit):
    indexes = config.indexes
    attributes = config.attributes[Attribute.age.value]

    indexes_new = []
    for index in indexes:
        if int(attributes[index]) < limit:
            indexes_new.append(index)

    config.indexes = indexes_new


def age_more(config, limit):
    indexes = config.indexes
    attributes = config.attributes[Attribute.age.value]

    indexes_new = []
    for index in indexes:
        if int(attributes[index]) >= limit:
            indexes_new.append(index)

    config.indexes = indexes_new