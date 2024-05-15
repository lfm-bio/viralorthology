class Args:
    def __init__(self):
        self.proteinortho = ''
        self.orffinder = ''
        self.blastp = ''
        self.HMMsearch = ''

    def get_software_params(self, software_params):
        return software_params if software_params != '' else 'default'

    def get_all_params(self):
        return f'proteinortho: {self.get_software_params(self.proteinortho)}\norfinder: {self.get_software_params(self.orffinder)}\nblastp: {self.get_software_params(self.blastp)}\nHMMsearch: {self.get_software_params(self.HMMsearch)}'

    def __repr__(self):
        return f'{self.proteinortho} {self.orffinder} {self.blastp} {self.HMMsearch}'