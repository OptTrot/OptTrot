from setuptools import setup, Extension, find_packages
from pathlib import Path

FILE_PATH = Path(__file__).parent
src_name = "opttrot"
def get_path(st:str):
    return str((FILE_PATH/st).absolute())

sources_list = [
    f'{src_name}/c_src/pauli_bn.c',
    f'{src_name}/c_src/bn/bn.c',
    f'{src_name}/c_src/bn/bn_ext.c',
    f'{src_name}/c_src/bn/bn_python.c'
]
dirs_list = [
    f'{src_name}',
    f'{src_name}/c_src',
    f'{src_name}/c_src/bn'
]

# Define the extension module
pauli_module = Extension(
    'pauli_c',
    sources=[get_path(s) for s in sources_list],
    include_dirs=[get_path(s) for s in dirs_list],  # Include directories for header files
    define_macros=[("BIG_NUM_BYTES", "512")]
)

# Run the setup
setup(
    name='opttrot',
    version='0.1',
    description='Pauli element manipulation library',
    packages=find_packages(where=src_name),
    package_dir={'': src_name},
    ext_modules=[pauli_module],
    py_modules=[
        'opttrot'
        #'pauli',
        #'pauli_utils',
        #'utils'
    ]
)
