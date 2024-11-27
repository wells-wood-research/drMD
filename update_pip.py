import argpass
import os
from os import path as p
import subprocess
import re
import shutil


class FilePath:
    pass
class DirectoryPath:
    pass

def main():
    newVersion: str = get_input_version()
    # find files
    setupDotPy, pyprojectDotToml = find_setup_files()
    ## delete previous builds
    clean_previous_builds()
    ## edit setup files
    update_version(setupDotPy, newVersion)
    update_version(pyprojectDotToml, newVersion)
    build_package()
    upload_package()



def get_input_version():
    ## get new version as command line arg
    parser = argpass.ArgumentParser()
    parser.add_argument(f"--newVersion")
    args = parser.parse_args()
    newVersion: str = args.newVersion


    if newVersion == None:
        raise ValueError("ENTER A VERSION NUMBER USING --newVersion INPUT!")
    return newVersion

def find_setup_files() -> tuple[FilePath, FilePath]:
    print("Current working directory:", os.getcwd())

    cwd: DirectoryPath = os.getcwd()
    print("setup.py:", p.join(cwd, "setup.py"))
    setupDotPy: FilePath = p.join(cwd, "setup.py")
    print("pyproject.toml:", p.join(cwd, "pyproject.toml"))
    pyprojectDotToml: FilePath = p.join(cwd, "pyproject.toml")
    return setupDotPy, pyprojectDotToml

def build_package():
    """Build the package using setup.py."""
    print("Building the package...")
    try:
        result = subprocess.run(['python', 'setup.py', 'sdist', 'bdist_wheel'], check=True, capture_output=True, text=True)
        print(result.stdout)
        print("Package built successfully.")
    except subprocess.CalledProcessError as e:
        print("An error occurred while building the package:", e)
        print(e.stdout)
        print(e.stderr)

def upload_package():
    """Upload the package to PyPI using twine."""
    print("Uploading the package to PyPI using twine...")
    try:
        subprocess.run(['twine', 'upload', 'dist/*'], check=True)
        print("Package uploaded successfully.")
    except subprocess.CalledProcessError as e:
        print("An error occurred while uploading the package:", e)


def clean_previous_builds():
    """Remove previous build directories."""
    for directory in ['build', 'dist']:
        if os.path.exists(directory):
            shutil.rmtree(directory)
    for file in os.listdir('.'):
        if file.endswith('.egg-info'):
            shutil.rmtree(file)

def update_version(inputFile: FilePath, newVersion: str):
    """Update the version in pip pyproject.toml and setup.py"""
    with open(inputFile, 'r') as file:
        content = file.read()

    updated_content = re.sub(
        r'version\s*=\s*["\'][^"\']*["\']',
        f"version = '{newVersion}'",
        content
    )
    with open(inputFile, 'w') as file:
        file.write(updated_content)


if __name__ == "__main__":
    main()