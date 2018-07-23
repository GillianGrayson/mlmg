from infrastructure.file_system import *

def get_indexes(config):
    indexes = []
    if config.db is DataBaseType.GSE40279:
        genders = config.attributes['gender']
        if config.gt is Gender.any:
            indexes = list(range(0, len(genders)))
        else:
            for p_id in range(0, len(genders)):
                if config.gt.value == genders[p_id]:
                    indexes.append(p_id)
    elif config.db is DataBaseType.GSE52588:
        print('TODO')
    return indexes