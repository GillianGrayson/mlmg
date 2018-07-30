from infrastructure.file_system import *

def get_indexes(config):
    geo_accessions = config.attributes[Attribute.geo_accession.value]
    indexes = list(range(0, len(geo_accessions)))

    if Attribute.gender.value in config.attributes:
        gender_indexes = []
        genders = config.attributes[Attribute.gender.value]
        if config.gender is Gender.any:
            gender_indexes = list(range(0, len(genders)))
        else:
            for p_id in range(0, len(genders)):
                if config.gender.value == genders[p_id]:
                    gender_indexes.append(p_id)
        indexes = list(set(indexes).intersection(gender_indexes))

    if Attribute.disease in config.attributes:
        disease_indexes = []
        diseases = config.attributes[Attribute.disease.value]
        if config.disease is Disease.any:
            disease_indexes = list(range(0, len(diseases)))
        else:
            for p_id in range(0, len(diseases)):
                if config.disease.value == diseases[p_id]:
                    disease_indexes.append(p_id)
        indexes = list(set(indexes).intersection(disease_indexes))

    return indexes