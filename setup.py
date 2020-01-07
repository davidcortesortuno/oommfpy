import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

with open('README.md') as f:
    long_description = f.read()

com_args = ['-std=c99', '-O3', '-Wno-cpp', '-Wno-unused-function']
extensions = [
    Extension("oommfpy.tools.clib",
              ["oommfpy/tools/clib/clib.pyx",
               "oommfpy/tools/clib/vtk_writer.c"],
              extra_compile_args=com_args
              ),
]

setuptools.setup(
    name='oommfpy',
    version='0.9',
    description=('Python lib to read and process OOMMF data files'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='D. Cortes',
    author_email='d.cortes@soton.ac.uk',
    packages=setuptools.find_packages(),
    ext_modules=cythonize(extensions),
    install_requires=['matplotlib',
                      'numpy',
                      'click',
                      'cython'
                      ],
    classifiers=['License :: BSD2 License',
                 'Programming Language :: Python :: 3 :: Only',
                 ],
    entry_points='''
        [console_scripts]
        plot_omf=oommfpy.tools.plot_slices:plot_omf_slices_cli
        omf2vtk=oommfpy.tools.omf2vtk:omf2vtk_cli
    ''',
)
