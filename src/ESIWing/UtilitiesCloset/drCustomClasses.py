## BASIC PYTHON LIBRARIES
import os

class DirectoryPath:
    def __init__(self, path: str):
        if not os.path.isdir(path):
            raise ValueError(f"{path} is not a valid directory path")
        self.path = os.path.abspath(path)

    def __str__(self):
        return self.path

class FilePath:
    def __init__(self, path: str):
        if not os.path.isfile(path):
            raise ValueError(f"{path} is not a valid file path")
        self.path = os.path.abspath(path)

    def __str__(self):
        return self.path