import setuptools
from setuptools.extension import Extension
import sys

# Check that Cython is installed:
# https://github.com/pypa/setuptools/issues/1317
# TODO: Implement pyproject.toml file with dependencies
# See: https://pip.pypa.io/en/stable/reference/pip/#pep-517-and-518-support
try:
    from Cython.Build import cythonize
except ImportError:
    import subprocess

    errno = subprocess.call([sys.executable, "-m", "pip", "install", "cython"])
    if errno:
        print("Please install Cython")
        raise SystemExit(errno)
    else:
        from Cython.Build import cythonize

with open('README.md') as f:
    long_description = f.read()

com_args = ['-std=c99', '-O3']
extensions = [
    Extension("oommfpy.tools.clib",
              ["oommfpy/tools/clib/clib.pyx",
               "oommfpy/tools/clib/vtk_writer.c"],
              extra_compile_args=com_args
              ),
]

setuptools.setup(
    # setup_requires=['cython'],  # not working (see the link at top)
    name='oommfpy',
    version='1.0a',
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
