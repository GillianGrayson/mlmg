from enum import Enum


class Attribute(Enum):
    title = 'title'
    geo_accession = 'geo_accession'
    source_name_ch1 = 'source_name_ch1'
    tissue = 'tissue'
    age = 'age'
    gender = 'gender'
    disease = 'disease'
    group = 'Group'
    batch = 'Batch'
    dnam_age = 'dnam_age'
    age_acceleration_ds_vs_controls = 'age_acceleration_ds_vs_controls'
    age_acceleration_ds_vs_controls_in_crbm = 'age_acceleration_ds_vs_controls_in_crbm'
    age_acceleration_ds_vs_controls_in_frontal_cortex = 'age_acceleration_ds_vs_controls_in_frontal_cortex'
    age_acceleration_ds_vs_controls_other_regions = 'age_acceleration_ds_vs_controls_other_regions'
    age_acceleration_residual = 'age_acceleration_residual'
    post_mortem_interval = 'post_mortem_interval'
    supplementary_file = 'supplementary_file'