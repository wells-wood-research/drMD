import argpass
import os
from os import path as p
import subprocess
import re
import shutil

def main():
    newVersion = get_input_version()
    # find files
    setupDotPy, pyprojectDotToml = find_setup_files()
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
        raise ValueError("ENTER A VERISON NUMBER!")
    return newVersion

def find_setup_files():
    print("Current working directory:", os.getcwd())

    cwd = os.getcwd()
    print("setup.py:", p.join(cwd, "setup.py"))
    setupDotPy = p.join(cwd, "setup.py")
    print("pyproject.toml:", p.join(cwd, "pyproject.toml"))
    pyprojectDotToml = p.join(cwd, "pyproject.toml")
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

def update_version(file_path, new_version):
    print(file_path)
    with open(file_path, 'r') as file:
        content = file.read()

    updated_content = re.sub(
        r'version\s*=\s*["\'][^"\']*["\']',
        f"version = '{new_version}'",
        content
    )
    
    with open(file_path, 'w') as file:
        file.write(updated_content)


if __name__ == "__main__":
    main()