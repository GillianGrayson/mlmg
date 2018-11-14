from enum import Enum


class Disease(Enum):
    any = 'any'
    healthy = 'healthy'
    down_syndrome = 'down_syndrome'
    versus = 'versus'


class Gender(Enum):
    M = 'M'
    F = 'F'
    any = 'any'
    versus = 'versus'