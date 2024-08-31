class DivisionError(Exception):
    def __init__(self, message):
        self.message = message

class FortranBinaryError(Exception):
    def __init__(self, message):
        self.message = message
