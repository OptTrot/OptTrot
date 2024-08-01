from setuptools import setup, Extension
from os import getcwd

print(getcwd())

# Define the extension module
pauli_module = Extension(
    'pauli_c',
    sources=[
        'src/c_src/pauli_bn.c',
        'src/c_src/bn/bn.c',
        'src/c_src/bn/bn_ext.c',
        'src/c_src/bn/bn_python.c'
    ],
    include_dirs=['src', 'src/c_src', 'src/c_src/bn'],  # Include directories for header files
    define_macros=[("BIG_NUM_BYTES", "512")]
)

# Run the setup
setup(
    name='pauli_c',
    version='0.1',
    description='Pauli element manipulation library',
    ext_modules=[pauli_module],
)
