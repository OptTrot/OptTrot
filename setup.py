from setuptools import setup, Extension, find_packages
from pathlib import Path
import numpy

sep_path = __file__.split("..")
FILE_PATH = Path(sep_path[0]).parent
print("File path of setup.py:", FILE_PATH)
src_name = "opttrot"
def get_path(st:str):
    return str((FILE_PATH/st).absolute())
print(get_path(src_name))

opt_source_path = FILE_PATH/src_name

sources_list1_path =[
    opt_source_path/"c_src"/"pauli_c.c",
    opt_source_path/"c_src"/"pauli_bn"/"pauli_bn.c",
    opt_source_path/"c_src"/"pauli_bn"/"pauli_bn_methods.c",
    opt_source_path/"c_src"/"pauli_bn"/"pauli_bn_utils.c",
    opt_source_path/"c_src"/"bn"/"bn.c",
    opt_source_path/"c_src"/"bn"/"bn_ext.c",
    opt_source_path/"c_src"/"bn"/"bn_python.c",
]
sources_list2_path = [opt_source_path/"c_src"/"c_utils.c"]
dirs_list_path = [
    opt_source_path,
    opt_source_path/"c_src",
    opt_source_path/"c_src"/"bn",
    opt_source_path/"c_src"/"pauli_bn",
]

#sources_list1 = [
#    f'{src_name}/c_src/pauli_c.c',
#    f'{src_name}/c_src/pauli_bn/pauli_bn.c',
#    f'{src_name}/c_src/pauli_bn/pauli_bn_methods.c',
#    f'{src_name}/c_src/pauli_bn/pauli_bn_utils.c',
#    f'{src_name}/c_src/bn/bn.c',
#    f'{src_name}/c_src/bn/bn_ext.c',
#    f'{src_name}/c_src/bn/bn_python.c'
#]
#sources_list2 = [
#    f'{src_name}/c_src/c_utils.c',
#]
#dirs_list = [
#    f'{src_name}',
#    f'{src_name}/c_src',
#    f'{src_name}/c_src/bn',
#    f'{src_name}/c_src/pauli_bn',
#    numpy.get_include()
#]

sources_list1 =[str(p.absolute()) for p in sources_list1_path]
sources_list2 = [str(p.absolute()) for p in sources_list2_path]
sources_list2 = sources_list1 + sources_list2
dirs_list = [str(p.absolute()) for p in dirs_list_path]
dirs_list.append(numpy.get_include())

print(sources_list1)
print(sources_list2)
print(dirs_list)


# Define the extension module
pauli_module = Extension(
    'pauli_c',
    sources=sources_list1, #[get_path(s) for s in sources_list1],
    include_dirs=dirs_list, #get_path(s) for s in dirs_list],  # Include directories for header files
    define_macros=[("BIG_NUM_BYTES", "512")]
)
c_utils_module = Extension(
    'c_utils',
    sources=sources_list2, #[get_path(s) for s in sources_list2],
    include_dirs=dirs_list #[get_path(s) for s in dirs_list],  # Include directories for header files
)

# Run the setup
setup(
    name='opttrot',
    version='0.1',
    description='Pauli element manipulation library',
    packages=find_packages(where=get_path(src_name)),
    package_dir={'': get_path(src_name)},
    ext_modules=[
        pauli_module, 
        c_utils_module
        ],
    py_modules=[
        'opttrot'
    ]
)
