from setuptools import setup, Extension, find_packages
from pathlib import Path
import numpy

sep_path = __file__.split("..")
FILE_PATH = Path(sep_path[0]).parent
print("First:", FILE_PATH)
src_name = "opttrot"
def get_path(st:str):
    return str((FILE_PATH/st).absolute())

print(get_path(src_name))
sources_list = [
    f'{src_name}/c_src/pauli_c.c',
    f'{src_name}/c_src/pauli_bn.c',
    f'{src_name}/c_src/pauli_bn_methods.c',
    f'{src_name}/c_src/pauli_bn_utils.c',
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
    include_dirs=[get_path(s) for s in dirs_list]+[numpy.get_include()],  # Include directories for header files
    define_macros=[("BIG_NUM_BYTES", "512")]
)

# Run the setup
setup(
    name='opttrot',
    version='0.1',
    description='Pauli element manipulation library',
    packages=find_packages(where=get_path(src_name)),
    package_dir={'': get_path(src_name)},
    ext_modules=[pauli_module],
    py_modules=[
        'opttrot'
    ]
)
