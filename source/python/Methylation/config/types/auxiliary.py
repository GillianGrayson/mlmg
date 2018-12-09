from enum import Enum


class MANOVATest(Enum):
    wilks = 'wilks'
    pillai_bartlett = 'pillai_bartlett'
    lawley_hotelling = 'lawley_hotelling'
    roy = 'roy'


class ClockExogType(Enum):
    all = 'all'
    single = 'single'
    slide = 'slide'